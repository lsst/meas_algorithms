#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
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
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#

__all__ = ("MeasureApCorrConfig", "MeasureApCorrTask", "MeasureApCorrError")

import numpy as np
from scipy.stats import median_abs_deviation

import lsst.pex.config
from lsst.afw.image import ApCorrMap
from lsst.afw.math import ChebyshevBoundedField, ChebyshevBoundedFieldConfig
from lsst.pipe.base import Task, Struct
from lsst.meas.base.apCorrRegistry import getApCorrNameSet

from .sourceSelector import sourceSelectorRegistry


class MeasureApCorrError(RuntimeError):
    pass


class _FluxNames:
    """A collection of flux-related names for a given flux measurement algorithm.

    Parameters
    ----------
    name : `str`
        Name of flux measurement algorithm, e.g. ``base_PsfFlux``.
    schema : `lsst.afw.table.Schema`
        Catalog schema containing the flux field.  The ``{name}_instFlux``,
        ``{name}_instFluxErr``, ``{name}_flag`` fields are checked for
        existence, and the ``apcorr_{name}_used`` field is added.

    Raises
    ------
    KeyError if any of instFlux, instFluxErr, or flag fields is missing.
    """
    def __init__(self, name, schema):
        self.fluxName = name + "_instFlux"
        if self.fluxName not in schema:
            raise KeyError("Could not find " + self.fluxName)
        self.errName = name + "_instFluxErr"
        if self.errName not in schema:
            raise KeyError("Could not find " + self.errName)
        self.flagName = name + "_flag"
        if self.flagName not in schema:
            raise KeyError("Cound not find " + self.flagName)
        self.usedName = "apcorr_" + name + "_used"
        schema.addField(self.usedName, type="Flag",
                        doc="Set if source was used in measuring aperture correction.")


class MeasureApCorrConfig(lsst.pex.config.Config):
    """Configuration for MeasureApCorrTask.
    """
    refFluxName = lsst.pex.config.Field(
        doc="Field name prefix for the flux other measurements should be aperture corrected to match",
        dtype=str,
        default="slot_CalibFlux",
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="Selector that sets the stars that aperture corrections will be measured from.",
        default="science",
    )
    minDegreesOfFreedom = lsst.pex.config.RangeField(
        doc="Minimum number of degrees of freedom (# of valid data points - # of parameters);"
        " if this is exceeded, the order of the fit is decreased (in both dimensions), and"
        " if we can't decrease it enough, we'll raise ValueError.",
        dtype=int,
        default=1,
        min=1,
    )
    fitConfig = lsst.pex.config.ConfigField(
        doc="Configuration used in fitting the aperture correction fields.",
        dtype=ChebyshevBoundedFieldConfig,
    )
    numIter = lsst.pex.config.Field(
        doc="Number of iterations for robust MAD sigma clipping.",
        dtype=int,
        default=4,
    )
    numSigmaClip = lsst.pex.config.Field(
        doc="Number of robust MAD sigma to do clipping.",
        dtype=float,
        default=4.0,
    )
    allowFailure = lsst.pex.config.ListField(
        doc="Allow these measurement algorithms to fail without an exception.",
        dtype=str,
        default=[],
    )

    def setDefaults(self):
        selector = self.sourceSelector["science"]

        selector.doFlags = True
        selector.doUnresolved = True
        selector.doSignalToNoise = True
        selector.doIsolated = False
        selector.flags.good = []
        selector.flags.bad = [
            "base_PixelFlags_flag_edge",
            "base_PixelFlags_flag_interpolatedCenter",
            "base_PixelFlags_flag_saturatedCenter",
            "base_PixelFlags_flag_crCenter",
            "base_PixelFlags_flag_bad",
            "base_PixelFlags_flag_interpolated",
            "base_PixelFlags_flag_saturated",
        ]
        selector.signalToNoise.minimum = 200.0
        selector.signalToNoise.maximum = None
        selector.signalToNoise.fluxField = "base_PsfFlux_instFlux"
        selector.signalToNoise.errField = "base_PsfFlux_instFluxErr"

    def validate(self):
        lsst.pex.config.Config.validate(self)
        if self.sourceSelector.target.usesMatches:
            raise lsst.pex.config.FieldValidationError(
                MeasureApCorrConfig.sourceSelector,
                self,
                "Star selectors that require matches are not permitted.")


class MeasureApCorrTask(Task):
    """Task to measure aperture correction.
    """
    ConfigClass = MeasureApCorrConfig
    _DefaultName = "measureApCorr"

    def __init__(self, schema, **kwargs):
        """Construct a MeasureApCorrTask

        For every name in lsst.meas.base.getApCorrNameSet():
        - If the corresponding flux fields exist in the schema:
            - Add a new field apcorr_{name}_used
            - Add an entry to the self.toCorrect dict
        - Otherwise silently skip the name
        """
        Task.__init__(self, **kwargs)
        self.refFluxNames = _FluxNames(self.config.refFluxName, schema)
        self.toCorrect = {}  # dict of flux field name prefix: FluxKeys instance
        for name in sorted(getApCorrNameSet()):
            try:
                self.toCorrect[name] = _FluxNames(name, schema)
            except KeyError:
                # if a field in the registry is missing, just ignore it.
                pass
        self.makeSubtask("sourceSelector")

    def run(self, exposure, catalog):
        """Measure aperture correction

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure aperture corrections are being measured on. The
            bounding box is retrieved from it, and it is passed to the
            sourceSelector. The output aperture correction map is *not*
            added to the exposure; this is left to the caller.
        catalog : `lsst.afw.table.SourceCatalog`
            SourceCatalog containing measurements to be used to
            compute aperture corrections.

        Returns
        -------
        Struct : `lsst.pipe.base.Struct`
            Contains the following:

            ``apCorrMap``
                aperture correction map (`lsst.afw.image.ApCorrMap`)
                that contains two entries for each flux field:
                - flux field (e.g. base_PsfFlux_instFlux): 2d model
                - flux sigma field (e.g. base_PsfFlux_instFluxErr): 2d model of error
        """
        bbox = exposure.getBBox()
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        doPause = lsstDebug.Info(__name__).doPause

        self.log.info("Measuring aperture corrections for %d flux fields", len(self.toCorrect))

        # First, create a subset of the catalog that contains only selected stars
        # with non-flagged reference fluxes.
        selected = self.sourceSelector.run(catalog, exposure=exposure)

        use = (
            ~selected.sourceCat[self.refFluxNames.flagName]
            & (np.isfinite(selected.sourceCat[self.refFluxNames.fluxName]))
        )
        goodRefCat = selected.sourceCat[use].copy()

        apCorrMap = ApCorrMap()

        # Outer loop over the fields we want to correct
        for name, fluxNames in self.toCorrect.items():
            # Create a more restricted subset with only the objects where the to-be-correct flux
            # is not flagged.
            fluxes = goodRefCat[fluxNames.fluxName]
            with np.errstate(invalid="ignore"):  # suppress NaN warnings.
                isGood = (
                    (~goodRefCat[fluxNames.flagName])
                    & (np.isfinite(fluxes))
                    & (fluxes > 0.0)
                )

            # The 1 is the minimum number of ctrl.computeSize() when the order
            # drops to 0 in both x and y.
            if (isGood.sum() - 1) < self.config.minDegreesOfFreedom:
                if name in self.config.allowFailure:
                    self.log.warning("Unable to measure aperture correction for '%s': "
                                     "only %d sources, but require at least %d.",
                                     name, isGood.sum(), self.config.minDegreesOfFreedom + 1)
                    continue
                else:
                    msg = ("Unable to measure aperture correction for required algorithm '%s': "
                           "only %d sources, but require at least %d." %
                           (name, isGood.sum(), self.config.minDegreesOfFreedom + 1))
                    self.log.warning(msg)
                    raise MeasureApCorrError("Aperture correction failed on required algorithm.")

            goodCat = goodRefCat[isGood].copy()

            x = goodCat['slot_Centroid_x']
            y = goodCat['slot_Centroid_y']
            z = goodCat[self.refFluxNames.fluxName]/goodCat[fluxNames.fluxName]
            ids = goodCat['id']

            # We start with an initial fit that is the median offset; this
            # works well in practice.
            fitValues = np.median(z)

            ctrl = self.config.fitConfig.makeControl()

            allBad = False
            for iteration in range(self.config.numIter):
                resid = z - fitValues
                # We add a small (epsilon) amount of floating-point slop because
                # the median_abs_deviation may give a value that is just larger than 0
                # even if given a completely flat residual field (as in tests).
                apCorrErr = median_abs_deviation(resid, scale="normal") + 1e-7
                keep = np.abs(resid) <= self.config.numSigmaClip * apCorrErr

                self.log.debug("Removing %d sources as outliers.", len(resid) - keep.sum())

                x = x[keep]
                y = y[keep]
                z = z[keep]
                ids = ids[keep]

                while (len(x) - ctrl.computeSize()) < self.config.minDegreesOfFreedom:
                    if ctrl.orderX > 0:
                        ctrl.orderX -= 1
                    else:
                        allBad = True
                        break
                    if ctrl.orderY > 0:
                        ctrl.orderY -= 1
                    else:
                        allBad = True
                        break

                if allBad:
                    if name in self.config.allowFailure:
                        self.log.warning("Unable to measure aperture correction for '%s': "
                                         "only %d sources remain, but require at least %d." %
                                         (name, keep.sum(), self.config.minDegreesOfFreedom + 1))
                        break
                    else:
                        msg = ("Unable to measure aperture correction for required algorithm "
                               "'%s': only %d sources remain, but require at least %d." %
                               (name, keep.sum(), self.config.minDegreesOfFreedom + 1))
                        self.log.warning(msg)
                        raise MeasureApCorrError("Aperture correction failed on required algorithm.")

                apCorrField = ChebyshevBoundedField.fit(bbox, x, y, z, ctrl)
                fitValues = apCorrField.evaluate(x, y)

            if allBad:
                continue

            self.log.info(
                "Aperture correction for %s from %d stars: MAD %f, RMS %f",
                name,
                median_abs_deviation(fitValues - z, scale="normal"),
                np.mean((fitValues - z)**2.)**0.5,
                len(x),
            )

            if display:
                plotApCorr(bbox, x, y, z, apCorrField, "%s, final" % (name,), doPause)

            # Record which sources were used.
            used = np.zeros(len(catalog), dtype=bool)
            used[np.searchsorted(catalog['id'], ids)] = True
            catalog[fluxNames.usedName] = used

            # Save the result in the output map
            # The error is constant spatially (we could imagine being
            # more clever, but we're not yet sure if it's worth the effort).
            # We save the errors as a 0th-order ChebyshevBoundedField
            apCorrMap[fluxNames.fluxName] = apCorrField
            apCorrMap[fluxNames.errName] = ChebyshevBoundedField(
                bbox,
                np.array([[apCorrErr]]),
            )

        return Struct(
            apCorrMap=apCorrMap,
        )


def plotApCorr(bbox, xx, yy, zzMeasure, field, title, doPause):
    """Plot aperture correction fit residuals

    There are two subplots: residuals against x and y.

    Intended for debugging.

    Parameters
    ----------
    bbox : `lsst.geom.Box2I`
        Bounding box (for bounds)
    xx, yy : `numpy.ndarray`, (N)
        x and y coordinates
    zzMeasure : `float`
        Measured value of the aperture correction
    field : `lsst.afw.math.ChebyshevBoundedField`
        Fit aperture correction field
    title : 'str'
        Title for plot
    doPause : `bool`
        Pause to inspect the residuals plot? If
        False, there will be a 4 second delay to
        allow for inspection of the plot before
        closing it and moving on.
    """
    import matplotlib.pyplot as plt

    zzFit = field.evaluate(xx, yy)
    residuals = zzMeasure - zzFit

    fig, axes = plt.subplots(2, 1)

    axes[0].scatter(xx, residuals, s=3, marker='o', lw=0, alpha=0.7)
    axes[1].scatter(yy, residuals, s=3, marker='o', lw=0, alpha=0.7)
    for ax in axes:
        ax.set_ylabel("ApCorr Fit Residual")
        ax.set_ylim(0.9*residuals.min(), 1.1*residuals.max())
    axes[0].set_xlabel("x")
    axes[0].set_xlim(bbox.getMinX(), bbox.getMaxX())
    axes[1].set_xlabel("y")
    axes[1].set_xlim(bbox.getMinY(), bbox.getMaxY())
    plt.suptitle(title)

    if not doPause:
        try:
            plt.pause(4)
            plt.close()
        except Exception:
            print("%s: plt.pause() failed. Please close plots when done." % __name__)
            plt.show()
    else:
        print("%s: Please close plots when done." % __name__)
        plt.show()
