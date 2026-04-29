# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ["ComputeRoughPsfShapeletsTask", "ComputeRoughPsfShapeletsConfig"]

import itertools
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import scipy.signal
import scipy.stats
from sklearn.covariance import MinCovDet
from sklearn.neighbors import KernelDensity

from lsst.afw.geom import ellipses
from lsst.afw.image import ImageD, ImageF, MaskedImageF
from lsst.afw.table import Point2DKey, QuadrupoleKey, Schema, SourceCatalog
from lsst.geom import Box2D
from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import AlgorithmError, Struct, Task
from lsst.shapelet import LAGUERRE, ShapeletFunction, computeOffset

from ._algorithmsLib import SpanSetMoments

if TYPE_CHECKING:
    import matplotlib.axes
    import matplotlib.figure
    import matplotlib.image
    import matplotlib.patches


class NoStarsForShapeletsError(AlgorithmError):
    """Exception raised when any of the methods involving selection fail to
    find any usable stars (e.g. compute_raw_moments,  _threshold_with_bounds,
    _find_first_radius_mode).
    """

    @property
    def metadata(self) -> dict[str, Any]:
        return {}


class ComputeRoughPsfShapeletsConfig(Config):
    bad_mask_planes = ListField(
        "Mask planes to identify pixels to drop from the calculation.",
        dtype=str,
        default=["SAT", "SUSPECT", "INTRP"],
    )
    bad_pixel_max_fraction = Field(
        "Maximum fraction of a footprint's pixels that can be bad (according "
        "to bad_mask_planes) before that footprint is fully dropped.",
        dtype=float,
        default=0.25,
    )
    bad_pixel_exclusion_radius = Field(
        "If a bad pixel (according to bad_mask_planes) falls within this "
        "radius of the unweighted centroid of a footprint, drop that footprint "
        "entirely.",
        dtype=float,
        default=2.0,
    )
    max_footprint_area = Field(
        "Footprints with a pixel count larger than this threshold are dropped before computing moments.",
        dtype=int,
        default=10000,
    )
    min_snr = Field(
        "Mininum flux S/N for inclusion in the star sample.",
        dtype=float,
        default=50.0,
    )
    max_radius_factor = Field(
        "Maximium multiple of the mode of the radius distribution for inclusion in the star sample.",
        dtype=float,
        default=2.0,
    )
    max_shape_distance = Field(
        "Maximum Mahalanobis distance (distance from the center of the shape "
        "distribution in elliptical sigma units) to select a star from the candidate "
        "sample, comparing the shape of that source vs. a robust estimate of the "
        "distribution of the shapes of all sources.",
        dtype=float,
        default=3.0,
    )
    min_n_stars = Field(
        "Minimum number of stars to select.  The S/N, radius, and shape distance thresholds are "
        "relaxed as needed to meet this target.",
        dtype=int,
        default=10,
    )
    max_n_stars = Field(
        "Maximum number of stars to select.  High shape-distance sources "
        "are dropped as needed to meet this target.",
        dtype=int,
        default=20,
    )
    logarithmic_shapes = Field(
        "If True, transform the (xx, yy, xy) moments to conformal shear "
        "and log trace radius when selecting stars in order to map the ellipse "
        "parameter to a space with (-inf, inf) bounds on all quantities (but "
        "a less-linear relationship to the pixel data).",
        dtype=bool,
        default=False,
    )
    radius_mode_min = Field(
        "Minimum radius in pixels at which to start searching for the first mode. "
        "This just needs be large enough to avoid unphysically small garbage (e.g. CRs).",
        dtype=float,
        default=1.0,
    )
    radius_kde_bandwidth = Field(
        "Bandwidth of the Gaussian kernel density estimator used to find the "
        "first mode of the radius distribution.",
        dtype=float,
        default=1.0,
    )
    shapelet_order = Field(
        "Order of the shapelet expansion fit to the stars.",
        dtype=int,
        default=4,
    )
    shapelet_scale_factor = Field(
        "Scale factor to apply to the moments ellipse when computing the ellipse for the shapelet basis.",
        dtype=float,
        default=1.0,
    )
    shapelet_circular_basis = Field(
        "Whether to use a circular shapelet basis with the same moments trace instead of an elliptical one.",
        dtype=bool,
        default=True,
    )

    def validate(self) -> None:
        if self.min_n_stars > self.max_n_stars:
            raise ValueError(
                f"min_n_stars={self.min_n_stars} is greater than max_n_stars={self.max_n_stars}."
            )
        if self.shapelet_order < 0:
            raise ValueError(f"shapelet order {self.shapelet_order} must be nonnegative.")
        if self.shapelet_scale_factor <= 0.0:
            raise ValueError(f"shapelet scale factor {self.shapelet_scale_factor} must be positive.")


class ComputeRoughPsfShapeletsTask(Task):
    """A task that computes a rough shapelet expansion of the PSF from a set
    of high S/N detections.

    Notes
    -----
    This task is expected to be run early in single-epoch processing - just
    after background subtraction and an initial high S/N detection phase, and
    before any deblending or measurement - in order to identify out-of-focus
    or otherwise bad PSFs.

    Given a background-subtracted `lsst.afw.image.MaskedImage`, an
    `lsst.afw.table.SourceCatalog` with footprints attached, and a random
    number generator seed, the `run` method will:

    - Compute the *unweighted* 0th-2nd moments of every non-child source over
      the footprint (except certain configurable masked pixels).  This is
      delegated to the `compute_raw_moments` method (which uses the C++
      `SpanSetMoments` class for the pixel-level processing).  Unweighted
      moments are used to avoid "latching onto" a small piece of PSF
      substructure, but can be much noiser than the Gaussian-weighed moments
      we usually use.

    - Select a "candidate" sample of sources with successfully measured
      moments that satisfy a S/N cut and a radius cut (determined from the
      first mode of the radius distribution, via kernel density estimation),
      and then use a robust covariance estimator (`scikit_learn.MinCovDet`)
      to select presumed isolated stars that are close to the center of that
      distribution, in 3-parameter shape space. This is delegated to the
      `select_stars` method.

    - Fit a single shapelet expansion to the selected stars.  This is
      mostly delegated to the `SpanSetMoments.fit_shapelets` method.

    The radial shapelet terms at 0th, 2nd, and 4th order are expected to form
    a space in which donut-shaped PSFs are well-separated from those with
    monotonic profiles.  Other terms *may* be useful in identifying other
    kinds of undesirable PSF structure.
    """

    ConfigClass = ComputeRoughPsfShapeletsConfig
    _DefaultName = "computeRoughPsfShapelets"
    config: ComputeRoughPsfShapeletsConfig

    def __init__(
        self,
        config: ComputeRoughPsfShapeletsConfig | None = None,
        *,
        schema: Schema,
        **kwargs: Any,
    ):
        super().__init__(config=config, **kwargs)
        self.schema = schema
        self._flux_key = schema.addField(
            "RoughPsfShapelets_flux",
            type=float,
            doc="Unweighted zeroth moment.",
        )
        self._flux_err_key = schema.addField(
            "RoughPsfShapelets_fluxErr",
            type=float,
            doc="Uncertainty on the unweighted zeroth moment.",
        )
        self._center_key = Point2DKey.addFields(
            schema, "RoughPsfShapelets", "Center from unweighted first moments.", "pixels"
        )
        self._shape_key = QuadrupoleKey.addFields(
            schema, "RoughPsfShapelets", "Shape from unweighted second moments."
        )
        self._flag_key = schema.addField(
            "RoughPsfShapelets_flag",
            type="Flag",
            doc="Flag set if the raw PSF moments were not computed.",
        )
        self._candidate_key = schema.addField(
            "RoughPsfShapelets_candidate",
            type="Flag",
            doc="Flag set if this source passed the radius_fraction cut (see configuration).",
        )
        self._used_key = schema.addField(
            "RoughPsfShapelets_used",
            type="Flag",
            doc=(
                "Flag set if this source passed the radius_fraction and shape_distance cuts "
                "(see configuration) and was used to fit the shapelet expansion."
            ),
        )

    def run(self, *, masked_image: MaskedImageF, catalog: SourceCatalog, seed: int) -> Struct:
        """Compute raw moments, select stars, and fit a shapelet expansion to
        them.

        Parameters
        ----------
        masked_image
            Masked image to measure on.  Must be background-subtracted.
        catalog
            Catalog of detections to extract footprints from and fill output
            columns of.  Its schema must be a superset of ``self.schema``.
        seed
            A random-number generator seed, used for the robust covariance
            estimator.

        Returns
        -------
        `lsst.pipe.base.Struct`
            A struct of results containing:

            - ``shapelet`` (`lsst.shapelet.ShapeletFunction`): A
              Gauss-Laguerre (polar shaplet) expansion of the PSF.
            - ``radial`` (`list` [`float`]): the purely radial coefficients
              of the shapelet expansion.
            - all attributes returned by the `select_stars` method.
        """
        moments = self.compute_raw_moments(masked_image=masked_image, catalog=catalog)
        result = self.select_stars(catalog, seed=seed)
        star_moments = [moments[star_id] for star_id in result.star_ids]
        result.shapelet = SpanSetMoments.fit_shapelets(
            masked_image,
            star_moments,
            self.config.shapelet_order,
            self.config.shapelet_scale_factor,
            self.config.shapelet_circular_basis,
        )
        result.shapelet.getEllipse().setCore(result.mean_shape)
        result.shapelet.changeBasisType(LAGUERRE)
        result.radial = result.shapelet.getCoefficients()[
            [computeOffset(i) for i in range(0, self.config.shapelet_order + 1, 2)]
        ]
        self.log.info("Rough PSF shapelet radial terms: %s.", result.radial)
        return result

    def compute_raw_moments(
        self, *, masked_image: MaskedImageF, catalog: SourceCatalog
    ) -> dict[int, SpanSetMoments]:
        """Compute the unweighted moments of the footprints in a catalog.

        Parameters
        ----------
        masked_image
            Masked image to measure on.  Must be background-subtracted.
        catalog
            Catalog of detections to extract footprints from and fill output
            columns of.  Its schema must be a superset of ``self.schema``.

        Returns
        -------
        `dict` [`int`, `SpanSetMoments`]
            Objects used to construct and hold the unweighted moments and the
            pixel region used to computed them, keyed by source ID.
        """
        bitmask = masked_image.mask.getPlaneBitMask(self.config.bad_mask_planes)
        all_moments: dict[int, SpanSetMoments] = {}
        for record in catalog:
            if record.getParent() != 0:
                record.set(self._flag_key, True)
                self.log.debug("Skipping child source %s", record.getId())
                continue
            footprint_spans = record.getFootprint().getSpans()
            if footprint_spans.getArea() > self.config.max_footprint_area:
                record.set(self._flag_key, True)
                self.log.debug(
                    "Skipping source %s with footprint area %d > %d.",
                    record.getId(),
                    footprint_spans.getArea(),
                    self.config.max_footprint_area,
                )
                continue
            moments = SpanSetMoments.compute(
                record.getFootprint().getSpans(),
                masked_image=masked_image,
                bad_bitmask=bitmask,
                bad_pixel_max_fraction=self.config.bad_pixel_max_fraction,
                bad_pixel_exclusion_radius=self.config.bad_pixel_exclusion_radius,
            )
            record.set(self._flux_key, moments.flux)
            record.set(self._flux_err_key, moments.variance**0.5)
            record.set(self._center_key, moments.center)
            record.set(self._shape_key, moments.shape)
            record.set(self._flag_key, moments.any_flags_set)
            if not moments.any_flags_set:
                all_moments[record.getId()] = moments
        if all_moments:
            self.log.verbose(
                "Successfully measured raw moments for %d of %d sources.",
                len(all_moments),
                len(catalog),
            )
        else:
            raise NoStarsForShapeletsError("No raw moments could be measured.")
        return all_moments

    def select_stars(self, catalog: SourceCatalog, seed: int) -> Struct:
        """Select probable stars from the distribution of second moments.

        Parameters
        ----------
        catalog
            Catalog of detections to extract footprints from and fill output
            columns of.  Its schema must be a superset of ``self.schema``.
        seed
            A random-number generator seed, used for the robust covariance
            estimator.

        Returns
        -------
        `lsst.pipe.base.Struct`
            A struct of results containing:

            - ``star_ids`` (`numpy.ndarray`): the source IDs that are expected
              to be stars.
            - ``mean_shape`` (`lsst.afw.geom.ellipses.BaseCore`): the mean of
              the shape distribution.
            - ``shape_covariance`` (`numpy.ndarray`): the covariance of the
              distribution of shapes; a 3x3 matrix.  This uses the same
              parameterization of the shapes as ``mean_shape``.
            - ``radius_cut`` (`float`): the indended radius cut (i.e. the
              mode of the radius distribution multipled by the
              ``radius_factor`` configuration option).
            - ``radius_kde`` (`sklearn.neighbors.KernelDensity`): kernel
              density estimator on the radius distribution, used to determine
              the radius cut.
        """
        # Cut on flags and SNR first.
        indices = np.arange(len(catalog), dtype=int)[np.logical_not(catalog[self._flag_key])]
        signalToNoise = []
        for index in indices:
            if (np.isfinite(catalog[self._flux_err_key][index])
               and np.isfinite(catalog[self._flux_err_key][index])
               and catalog[self._flux_err_key][index] != 0.0):
                signalToNoise.append(catalog[self._flux_key][index] / catalog[self._flux_err_key][index])
            else:
                signalToNoise.append(np.nan)
        signalToNoise = np.array(signalToNoise)
        indices = indices[
            self._threshold_with_bounds(
                signalToNoise,
                threshold=self.config.min_snr,
                min_count=self.config.min_n_stars,
                max_count=len(catalog),
                name="S/N",
                kind=">",
            )
        ]
        # Cut on radius next.
        radii = np.zeros(indices.size, dtype=np.float64)
        for n, index in enumerate(indices):
            record = catalog[index]
            shape = record.get(self._shape_key)
            radii[n] = shape.getTraceRadius()
        radius_mode, radius_kde = self._find_first_radius_mode(radii)
        radius_cut = self.config.max_radius_factor * radius_mode
        indices = indices[
            self._threshold_with_bounds(
                radii,
                threshold=radius_cut,
                min_count=self.config.min_n_stars,
                max_count=len(catalog),
                name="radius",
                kind="<",
            )
        ]
        shape_data = np.zeros((len(indices), 3), dtype=np.float64)
        for n, index in enumerate(indices):
            record = catalog[index]
            record.set(self._candidate_key, True)
            shape = record.get(self._shape_key)
            if self.config.logarithmic_shapes:
                shape = ellipses.SeparableConformalShearLogTraceRadius(shape)
            shape_data[n, :] = shape.getParameterVector()
        shape_dist = MinCovDet(random_state=seed).fit(shape_data)
        m_distances = shape_dist.mahalanobis(shape_data)
        indices = indices[
            self._threshold_with_bounds(
                m_distances,
                threshold=self.config.max_shape_distance,
                min_count=self.config.min_n_stars,
                max_count=self.config.max_n_stars,
                name="Mahalanobis distance",
                kind="<",
            )
        ]
        for index in indices:
            catalog[index].set(self._used_key, True)
        star_ids = catalog["id"][indices]
        if self.config.logarithmic_shapes:
            mean_shape = ellipses.SeparableConformalShearLogTraceRadius(shape_dist.location_)
        else:
            mean_shape = ellipses.Quadrupole(shape_dist.location_)
        return Struct(
            star_ids=star_ids,
            mean_shape=mean_shape,
            shape_covariance=shape_dist.covariance_,
            radius_cut=radius_cut,
            radius_kde=radius_kde,
            catalog=catalog,
        )

    def plot_selection(
        self, figure: matplotlib.figure.Figure, *, catalog: SourceCatalog, results: Struct
    ) -> None:
        """Create plots of the shape distribution space used to select stars.

        Parameters
        ----------
        figure
            Matplotlib figure to plot to.
        catalog
            Catalog of sources with columns populated by the `run` method (at
            least through the `select_stars` step).
        results
            Result struct returned by `run` or `select_stars`.
        """
        from matplotlib.lines import Line2D

        shape_data = np.zeros((len(catalog), 3), dtype=np.float64)
        radii = np.zeros(len(catalog), dtype=np.float64)
        for n, record in enumerate(catalog):
            if record[self._flag_key]:
                continue
            shape = record[self._shape_key]
            if self.config.logarithmic_shapes:
                shape = ellipses.SeparableConformalShearLogTraceRadius(shape)
            shape_data[n, :] = shape.getParameterVector()
            radii[n] = shape.getTraceRadius()
        used_mask = catalog[self._used_key]
        candidate_mask = np.logical_and(catalog[self._candidate_key], np.logical_not(used_mask))
        measured_mask = np.logical_and(
            np.logical_not(catalog[self._flag_key]), np.logical_not(catalog[self._candidate_key])
        )
        # Set up the axes.
        axes = figure.subplot_mosaic(
            [
                ["radius", "radius", "radius"],
                ["hist0", ".", "."],
                ["scatter01", "hist1", "."],
                ["scatter02", "scatter12", "hist2"],
            ],
            gridspec_kw=dict(bottom=0.2, top=1.0),
        )
        axes["scatter01"].sharex(axes["hist0"])
        axes["scatter02"].sharex(axes["hist0"])
        axes["scatter12"].sharex(axes["hist1"])
        axes["scatter12"].sharey(axes["scatter02"])
        for tk in itertools.chain(
            axes["hist0"].get_xticklabels(),
            axes["scatter01"].get_xticklabels(),
            axes["hist1"].get_xticklabels(),
            axes["scatter12"].get_yticklabels(),
        ):
            tk.set_visible(False)
        # Move hist y axes to the outside.
        for ax in [axes["hist0"], axes["hist1"], axes["hist2"]]:
            ax.yaxis.set_label_position("right")
            ax.yaxis.tick_right()
        # Add labels to outer axes.
        names = ["Ixx", "Iyy", "Ixy"] if not self.config.logarithmic_shapes else ["𝜂1", "𝜂2", "ln(r)"]
        axes["scatter02"].set_xlabel(names[0])
        axes["scatter12"].set_xlabel(names[1])
        axes["hist2"].set_xlabel(names[2])
        axes["scatter01"].set_ylabel(names[1])
        axes["scatter02"].set_ylabel(names[2])
        # Make the plots.
        mu = results.mean_shape.getParameterVector()
        sigma = np.diagonal(results.shape_covariance) ** 0.5
        lower_bounds = [max(mu[i] - 3 * sigma[i], min(shape_data[:, i])) for i in range(3)]
        upper_bounds = [min(mu[i] + 3 * sigma[i], max(shape_data[:, i])) for i in range(3)]
        grids = [np.linspace(lower_bounds[i], upper_bounds[i], 50) for i in range(3)]
        axes["radius"].axvline(results.radius_cut, color="k")
        for color, mask, alpha in [
            ("grey", measured_mask, 0.5),
            ("blue", candidate_mask, 0.75),
            ("green", used_mask, 1.0),
        ]:
            axes["radius"].hist(
                radii[mask],
                color=color,
                alpha=alpha,
                bins=16,
                histtype="step",
                range=(radii.min(), 15.0),
                density=True,
            )
            for i in range(3):
                axes[f"hist{i}"].hist(
                    shape_data[mask, i],
                    bins=16,
                    range=(lower_bounds[i], upper_bounds[i]),
                    density=True,
                    color=color,
                    histtype="step",
                    alpha=alpha,
                )
                for j in range(3):
                    if (ax := axes.get(f"scatter{i}{j}")) is not None:
                        ax.scatter(
                            shape_data[mask, i],
                            shape_data[mask, j],
                            c=color,
                            s=4,
                            alpha=alpha,
                            edgecolors=None,
                        )
        for i in range(3):
            axes[f"hist{i}"].plot(
                grids[i], scipy.stats.norm.pdf(grids[i], loc=mu[i], scale=sigma[i]), "k", alpha=0.5
            )
            axes[f"hist{i}"].set_xlim(lower_bounds[i], upper_bounds[i])
            for j in range(3):
                if (ax := axes.get(f"scatter{i}{j}")) is not None:
                    sigma_ellipse = ellipses.Quadrupole(
                        results.shape_covariance[i, i],
                        results.shape_covariance[j, j],
                        results.shape_covariance[i, j],
                    )
                    for factor in [1, 2, 3]:
                        self._draw_ellipse(
                            ax,
                            sigma_ellipse,
                            x=mu[i],
                            y=mu[j],
                            scale=factor,
                            fill=False,
                            edgecolor="k",
                            alpha=0.5,
                        )
                        ax.set_xlim(lower_bounds[i], upper_bounds[i])
                        ax.set_ylim(lower_bounds[j], upper_bounds[j])
        figure.legend(
            [
                Line2D([], [], color="green", alpha=1.0),
                Line2D([], [], color="blue", alpha=0.75),
                Line2D([], [], color="gray", alpha=0.5),
            ],
            [
                f"RoughPsfShapelets_used ({used_mask.sum()})",
                f"RoughPsfShapelets_candidate & ~RoughPsfShapelets_used ({candidate_mask.sum()})",
                f"~RoughPsfShapelets_flag & ~RoughPsfShapelets_candidate ({measured_mask.sum()})",
            ],
            loc="lower center",
        )
        return figure

    def plot_shapelets(
        self,
        figure: matplotlib.figure.Figure,
        *,
        image: ImageF,
        catalog: SourceCatalog,
        results: Struct,
        n_stars: int = 3,
        stamp_size: float = 2.0,
    ) -> None:
        """Create data/model/residual plots of stars and the shapelet model.

        Parameters
        ----------
        figure
            Matplotlib figure to plot to.
        image
            The image the stars were measured on.
        catalog
            Catalog of sources with columns populated by the `run` method .
        results
            Result struct returned by `run`.
        n_stars
            Number of stars to include.
        stamp_size
            Stamp size in inches.
        """
        from matplotlib.colors import Normalize

        width = stamp_size * 3 + 1.5
        figure.set_size_inches(w=width, h=stamp_size * n_stars)
        axes = figure.subplot_mosaic(
            [
                ["image_cbar", f"d{star_id}", f"m{star_id}", f"r{star_id}", "res_cbar"]
                for star_id in results.star_ids[:n_stars]
            ],
            gridspec_kw=dict(
                wspace=0.01, hspace=0.01, left=0.5 / width, right=1.0 - 0.5 / width, bottom=0.01, top=0.99
            ),
            width_ratios=[0.25, stamp_size, stamp_size, stamp_size, 0.25],
        )
        for name, ax in axes.items():
            if not name.endswith("cbar"):
                ax.axis("off")
        norm: Normalize | None = None
        res_norm: Normalize | None = None
        for star_id in results.star_ids[:n_stars]:
            record = catalog.find(star_id)
            star_bbox = record.getFootprint().getBBox()
            star_model = ImageD(star_bbox)
            star_center = record[self._center_key]
            star_ellipse = ellipses.Ellipse(ellipses.Axes(record[self._shape_key]), star_center)
            star_shapelet = ShapeletFunction(results.shapelet)
            star_shapelet.setEllipse(star_ellipse)
            star_shapelet.evaluate().addToImage(star_model)
            if norm is None:
                norm = Normalize(vmin=star_model.array.min(), vmax=star_model.array.max())
            star_image = image[star_bbox].clone()
            star_image /= record[self._flux_key]
            self._draw_image(axes[f"d{star_id}"], star_image, norm=norm, cmap="YlGnBu")
            self._draw_ellipse(axes[f"d{star_id}"], star_ellipse, fill=False, edgecolor="blue", alpha=0.5)
            image_plot = self._draw_image(axes[f"m{star_id}"], star_model, norm=norm, cmap="YlGnBu")
            self._draw_ellipse(axes[f"m{star_id}"], star_ellipse, fill=False, edgecolor="blue", alpha=0.5)
            star_image -= star_model.convertF()
            amax = np.abs(star_image.array).max()
            if res_norm is None:
                res_norm = Normalize(vmin=-amax, vmax=amax)
            res_plot = self._draw_image(axes[f"r{star_id}"], star_image, norm=res_norm, cmap="RdBu")
            self._draw_ellipse(axes[f"r{star_id}"], star_ellipse, fill=False, edgecolor="blue", alpha=0.5)
        figure.colorbar(image_plot, cax=axes["image_cbar"], location="left")
        figure.colorbar(res_plot, cax=axes["res_cbar"], location="right")
        return figure

    def _threshold_with_bounds(
        self,
        values: np.ndarray,
        threshold: float,
        min_count: int,
        max_count: int,
        name: str,
        kind: Literal["<", ">"],
    ) -> np.ndarray:
        """Return the indices of an array that satisfy an inequality
        and/or lower and upper bounds on the number of indices returned.

        Parameters
        ----------
        values
            Array of values to threshold on.
        threshold
            Threshold value that selected elements must be above or below.
        min_count
            The minimum number of indices returned.  When thresholding would
            yield fewer than this number, the threshold is ignored.  Note that
            the number of indices may still be less than this if the size of
            ``values`` is less than this.
        max_count
            The maximum number of indices returned.
        name
            Name of the quantity being thresholded, for log messages.
        kind
            Whether the threshold is a upper bound (``<``) or lower bound
            (``>``).  This also sets how values are ranked when they are added
            or dropped to satisfy the count constraints.

        Returns
        -------
        indices
            Indices into ``values``.
        """
        if min_count > len(values):
            raise NoStarsForShapeletsError(
                f"Not enough sources ({len(values)}) for {name} cut that must yield at least {min_count}."
            )
        sorter = values.argsort()
        n = np.searchsorted(values[sorter], threshold)
        if kind == ">":
            sorter = sorter[::-1]
            n = len(sorter) - n
        if n < min_count:
            self.log.verbose(
                "Applying a %s %s %f cut yields only %d sources; keeping the top %d (%s %s %f) instead.",
                name,
                kind,
                threshold,
                n,
                min_count,
                name,
                kind,
                values[sorter[min_count - 1]],
            )
            n = min_count
        elif n > max_count:
            self.log.verbose(
                "%d sources have %s %s %f; keeping only the top %d (%s %s %f) instead.",
                n,
                name,
                kind,
                threshold,
                max_count,
                name,
                kind,
                values[sorter[max_count - 1]],
            )
            n = max_count
        else:
            self.log.verbose("Keeping %d sources with %s %s %f.", n, name, kind, threshold)
        return sorter[:n]

    def _find_first_radius_mode(self, radii: np.ndarray) -> tuple[float, KernelDensity]:
        """Find the first peak in a 1-d distribution of radii."""
        kde = KernelDensity(bandwidth=self.config.radius_kde_bandwidth).fit(radii.reshape(-1, 1))
        sorted_radii = radii.copy()
        sorted_radii.sort()
        sorted_radii = sorted_radii[sorted_radii.searchsorted(self.config.radius_mode_min):]
        scores = kde.score_samples(sorted_radii.reshape(-1, 1))
        peaks, _ = scipy.signal.find_peaks(scores)
        if not peaks.size:
            raise NoStarsForShapeletsError("Radius distribute has no mode.")
        return sorted_radii[peaks.min()], kde

    @staticmethod
    def _draw_ellipse(
        axes: matplotlib.axes.Axes,
        ellipse: ellipses.BaseCore | ellipses.Ellipse,
        *,
        x: float | None = None,
        y: float | None = None,
        scale: float = 1.0,
        **kwargs: Any,
    ) -> matplotlib.patches.Ellipse:
        from matplotlib.patches import Ellipse as EllipsePatch

        if isinstance(ellipse, ellipses.Ellipse):
            if x is None:
                x = ellipse.getCenter().getX()
            if y is None:
                y = ellipse.getCenter().getY()
            ellipse = ellipse.getCore()
        else:
            if x is None:
                x = 0.0
            if y is None:
                y = 0.0
        ellipse = ellipses.Axes(ellipse)
        patch = EllipsePatch(
            (x, y),
            ellipse.getA() * 2 * scale,  # factor of 2 for radius->diameter
            ellipse.getB() * 2 * scale,
            angle=ellipse.getTheta() * 180.0 / np.pi,
            **kwargs,
        )
        axes.add_patch(patch)
        return patch

    @staticmethod
    def _draw_image(axes: matplotlib.axes.Axes, image: ImageF, **kwargs: Any) -> matplotlib.image.AxesImage:
        fp_bbox = Box2D(image.getBBox())
        return axes.imshow(
            image.array,
            interpolation="nearest",
            origin="lower",
            aspect="equal",
            extent=(fp_bbox.x.min, fp_bbox.x.max, fp_bbox.y.min, fp_bbox.y.max),
            **kwargs,
        )
