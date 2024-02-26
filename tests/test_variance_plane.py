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

import re
import unittest
from contextlib import nullcontext

import galsim
import lsst.utils.tests
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from lsst.ip.isr.isrMock import IsrMock
from lsst.meas.algorithms import remove_signal_from_variance
from lsst.utils.tests import methodParametersProduct
from matplotlib.legend_handler import HandlerTuple
from matplotlib.patheffects import withStroke
from matplotlib.ticker import FixedLocator, FuncFormatter

# Set to True to save the plot of the variance plane before and after
# correction for a representative test case.
SAVE_PLOT = True


def outline_effect(lw, alpha=0.8):
    """Generate a path effect for enhanced text visibility.

    Parameters
    ----------
    lw : `float`
        Line width of the outline.
    alpha : `float`, optional
        Transparency of the outline.

    Returns
    -------
    `list`
        A list containing the withStroke path effect.
    """
    return [withStroke(linewidth=lw, foreground="white", alpha=alpha)]


class CustomHandler(HandlerTuple):
    """Custom handler for handling grouped items in the legend."""

    def create_artists(self, *args):
        artists = super().create_artists(*args)
        for a in artists:
            a.set_transform(args[-1])
        return artists


def get_valid_color(handle):
    """Extracts a valid color from a Matplotlib handle.

    Parameters
    ----------
    handle : `matplotlib.artist.Artist`
        The handle from which to extract the color.

    Returns
    -------
    color : `str` or `tuple`
        The color extracted from the handle, or "default" if no valid color is
        found.
    """
    for attr in ["get_facecolor", "get_edgecolor", "get_color"]:
        if hasattr(handle, attr):
            color = getattr(handle, attr)()
            # If the handle is a collection, use the first color.
            if isinstance(color, np.ndarray) and color.shape[0] > 0:
                color = color[0]
            # If the color is RGBA with alpha = 0, continue the search.
            if len(color) == 4 and color[3] == 0:
                continue
            return color
    return "default"  # If no valid color is found


def get_emptier_side(ax):
    """Analyze a matplotlib Axes object to determine which side (left or right)
    has more whitespace, considering cases where artists' bounding boxes span
    both sides of the midpoint.

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The Axes object to analyze.

    Returns
    -------
    more_whitespace_side : `str`
        'left' if the left side has more whitespace, or 'right' if the right
        side does.
    """
    # Get the total plotting area's midpoint on the x-axis.
    xlim = ax.get_xlim()
    midpoint = sum(xlim) / 2

    # Initialize areas as zero
    left_area, right_area = 0, 0

    # Loop through all children (artists) in the Axes
    for artist in ax.get_children():
        # Skip if artist is invisible or lacks a bounding box.
        if not artist.get_visible() or not hasattr(artist, "get_window_extent"):
            continue
        bbox = artist.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())
        # Check if the artist's bounding box spans the midpoint.
        if bbox.x0 < midpoint < bbox.x1:
            # Calculate the proportion of the bbox on each side of the
            # midpoint.
            left_proportion = (midpoint - bbox.x0) / bbox.width
            right_proportion = 1 - left_proportion
            # Adjust area calculations for both sides.
            left_area += bbox.width * bbox.height * left_proportion
            right_area += bbox.width * bbox.height * right_proportion
        elif bbox.x0 + bbox.width / 2 < midpoint:
            # Entirely on the left.
            left_area += bbox.width * bbox.height
        else:
            # Entirely on the right.
            right_area += bbox.width * bbox.height

    # Determine which side has more whitespace by comparing occupied areas.
    return "left" if left_area <= right_area else "right"


def adjust_legend_with_groups(ax, combine_groups, colors="default", yloc="upper", **kwargs):
    """Adjusts the legend of a given Axes object by combining specified handles
    based on provided groups, setting the marker location and text alignment
    based on the inferable emptier side, and optionally setting the text color
    of legend entries to a provided list or inferring colors from the handles.
    Additionally, allows specifying the vertical location of the legend within
    the plot.

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The Axes object for which to adjust the legend.
    combine_groups : `list` of `int` or iterable of `int`
        A list that can contain a mix of individual integers and/or iterables
        (lists, tuples, or sets) of integers. An individual integer specifies
        the index of a single legend entry. An iterable of integers specifies
        a group of indices to be combined into a single legend entry.
    colors : `list` of `str` or `tuple`, or `str`, optional
        Specifies the colors for the legend entries. This parameter can be:
        - A list of color specifications, where each element is a string (for
          named colors or hex values) or a tuple (for RGB or RGBA values). This
          list explicitly assigns colors to each legend entry post-combination.
        - A single string value:
            - "match": Colors are inferred from the properties of the first
              handle in each group that corresponds to a non-white-space label.
              This aims to match the legend text color with the color of the
              plotted data.
            - "default": The function does not alter the default colors
              assigned by Matplotlib, preserving the automatic color assignment
              for all legend entries.
    yloc : `str`, optional
        The vertical location of the legend within the Axes. Valid options are
        'upper', 'lower', or 'middle'. This parameter is combined with the
        inferable emptier side ('left' or 'right') to determine the legend's
        placement. For example, 'upper right' or 'lower left'.
    **kwargs :
        Keyword arguments forwarded to the `ax.legend` function.
    """

    handles, labels = ax.get_legend_handles_labels()
    new_handles = []
    new_labels = []

    if colors == "match":
        colors = []
        infer_colors = True
    else:
        infer_colors = False

    for group in combine_groups:
        # Assume the first non-white-space label represents the group. If no
        # such label is found, just use the first label in the group which is
        # in fact empty.
        if isinstance(group, (list, tuple, set)):
            group = list(group)  # Just in case
            label_index = next((i for i in group if labels[i].strip()), group[0])
            combined_handle = tuple(handles[i] for i in group)
            combined_label = labels[label_index]
        elif isinstance(group, int):
            label_index = group
            combined_handle = handles[group]
            combined_label = labels[group]
        else:
            raise ValueError("Invalid value in 'combine_groups'")
        new_handles.append(combined_handle)
        new_labels.append(combined_label)
        if infer_colors:
            # Attempt to infer color from the representative handle in the
            # group.
            handle = handles[label_index]
            color = get_valid_color(handle)
            colors.append(color)

    # Determine the emptier side to decide legend and text alignment.
    emptier_side = get_emptier_side(ax)
    markerfirst = emptier_side != "right"

    # Create the legend with custom adjustments.
    legend = ax.legend(
        new_handles,
        new_labels,
        handler_map={tuple: CustomHandler()},
        loc=f"{yloc} {emptier_side}",
        fontsize=8,
        frameon=False,
        markerfirst=markerfirst,
        **kwargs,
    )

    # Right- or left-align the legend text based on the emptier side.
    for text in legend.get_texts():
        text.set_ha(emptier_side)

    # Set legend text colors if necessary.
    if colors != "default":
        for text, color in zip(legend.get_texts(), colors):
            if not (isinstance(color, str) and color == "default"):
                text.set_color(color)


def adjust_tick_scale(ax, axis_label_templates):
    """Scales down tick labels to make them more readable and updates axis
    labels accordingly.

    Calculates a power of 10 scale factor (common divisor) to reduce the
    magnitude of tick labels. It automatically determines which axes to adjust
    based on the provided axis label templates, which should include `{scale}`
    for inserting the scale factor dynamically.

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The Axes object to modify.
    axis_label_templates : `dict`
        Templates for axis labels, including `{scale}` for scale factor
        insertion. Keys should be one or more of axes names ("x", "y", "z")
        and values should be the corresponding label templates.
    """

    def trailing_zeros(n):
        """Determines the number of trailing zeros in a number."""
        return len(n := str(int(float(n))) if float(n).is_integer() else n) - len(n.rstrip("0"))

    def format_tick(val, pos, divisor):
        """Formats tick labels using the determined divisor."""
        return str(int(val / divisor)) if (val / divisor).is_integer() else str(val / divisor)

    # Iterate through the specified axes and adjust their tick labels and axis
    # labels.
    for axis in axis_label_templates.keys():
        # Gather current tick labels.
        labels = [label.get_text() for label in getattr(ax, f"get_{axis}ticklabels")()]

        # Calculate the power of 10 divisor based on the minimum number of
        # trailing zeros in the tick labels.
        divisor = 10 ** min(trailing_zeros(label) for label in labels if float(label) != 0)

        # Set a formatter for the axis ticks that scales them according to the
        # common divisor.
        getattr(ax, f"{axis}axis").set_major_formatter(
            FuncFormatter(lambda val, pos: format_tick(val, pos, divisor))
        )

        # Ensure the tick positions remain unchanged despite the new
        # formatting.
        getattr(ax, f"{axis}axis").set_major_locator(FixedLocator(getattr(ax, f"get_{axis}ticks")()))

        # Prepare 'scale', empty if divisor <= 1.
        scale = f"{int(divisor)}" if divisor > 1 else ""

        # Fetch the corresponding label template for the axis.
        label_template = axis_label_templates[axis]

        # If 'scale' is empty, remove whitespace around "{scale}" in the
        # template. Also remove any trailing "/{scale}".
        if scale == "":
            label_template = re.sub(r"\s*{\s*scale\s*}\s*", "{scale}", label_template)
            label_template = label_template.replace("/{scale}", "")

        # Always strip remaining whitespace from the template.
        label_template = label_template.strip()

        # Set the formatted axis label.
        label_text = label_template.format(scale=scale)
        getattr(ax, f"set_{axis}label")(label_text, labelpad=8)


class VariancePlaneTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        # Testing with a single detector that has 8 amplifiers in a 4x2
        # configuration. Each amplifier measures 100x51 in dimensions.
        config = IsrMock.ConfigClass()
        config.isLsstLike = True
        config.doAddBias = False
        config.doAddDark = False
        config.doAddFlat = False
        config.doAddFringe = False
        config.doGenerateImage = True
        config.doGenerateData = True
        config.doGenerateAmpDict = True
        self.mock = IsrMock(config=config)

    def tearDown(self):
        del self.mock

    def buildExposure(
        self,
        average_gain,
        gain_sigma_factor,
        sky_level,
    ):
        """Build and return an exposure with different types of simulated
        source profiles and a background sky level. It's intended for testing
        and analysis, providing a way to generate exposures with controlled
        conditions.

        Parameters
        ----------
        average_gain : `float`
            The average gain value of amplifiers in e-/ADU.
        gain_sigma_factor : float
            The standard deviation of the gain values as a factor of the
            'average_gain'.
        sky_level : `float`
            The background sky level in e-/arcsec^2.

        Returns
        -------
        exposure : `~lsst.afw.image.Exposure`
            An exposure object with simulated sources and background. The units
            are in detector counts (ADU).
        """

        # Set the random seed for reproducibility.
        random_seed = galsim.BaseDeviate(1905).raw() + 1
        np.random.seed(random_seed)
        rng = galsim.BaseDeviate(random_seed)

        # Get the exposure, detector, and amps from the mock.
        exposure = self.mock.getExposure()
        detector = exposure.getDetector()
        amps = detector.getAmplifiers()
        num_amps = len(amps)
        table = str.maketrans("", "", ":,")  # Remove ':' and ',' from names
        self.amp_names = [amp.getName().translate(table) for amp in amps]

        # Adjust instrument and observation parameters to some nominal values.
        pixel_scale = 0.2  # arcsec/pixel
        self.background = sky_level * pixel_scale**2  # e-/pixel

        # Get the bounding boxes for the exposure and amplifiers and convert
        # them to galsim bounds.
        exp_bbox = exposure.getBBox()
        image_bounds = galsim.BoundsI(exp_bbox.minX, exp_bbox.maxX, exp_bbox.minY, exp_bbox.maxY)
        self.amp_bbox_list = [amp.getBBox() for amp in amps]
        amp_bounds_list = [galsim.BoundsI(b.minX, b.maxX, b.minY, b.maxY) for b in self.amp_bbox_list]

        # Generate random deviations from the average gain across amplifiers
        # and adjust them to ensure their sum equals zero. This reflects
        # real-world detectors, with amplifier gains normally distributed due
        # to manufacturing and operational variations.
        deviations = np.random.normal(average_gain, gain_sigma_factor * average_gain, size=num_amps)
        deviations -= np.mean(deviations)

        # Set the gain for amplifiers to be slightly different from each other
        # while averaging to `average_gain`. This is to test the
        # `average_across_amps` option in the `remove_signal_from_variance`
        # function.
        self.amp_gain_list = [average_gain + deviation for deviation in deviations]

        # Define parameters for a mix of source types, including extended
        # sources with assorted profiles as well as point sources simulated
        # with minimal half-light radii to resemble hot pixels
        # post-deconvolution. All flux values are given in electrons and
        # half-light radii in pixels. The goal is for each amplifier to
        # predominantly contain at least one source, enhancing the
        # representativeness of test conditions.
        source_params = [
            {"type": "Sersic", "n": 3, "flux": 1.6e5, "half_light_radius": 3.5, "g1": -0.3, "g2": 0.2},
            {"type": "Sersic", "n": 1, "flux": 9.3e5, "half_light_radius": 2.1, "g1": 0.25, "g2": 0.12},
            {"type": "Sersic", "n": 4, "flux": 1.0e5, "half_light_radius": 1.1, "g1": 0.0, "g2": 0.0},
            {"type": "Sersic", "n": 3, "flux": 1.1e6, "half_light_radius": 4.2, "g1": 0.0, "g2": 0.2},
            {"type": "Sersic", "n": 5, "flux": 1.1e5, "half_light_radius": 3.6, "g1": 0.22, "g2": -0.05},
            {"type": "Sersic", "n": 2, "flux": 4.3e5, "half_light_radius": 2.0, "g1": 0.0, "g2": 0.0},
            {"type": "Sersic", "n": 6, "flux": 1.2e6, "half_light_radius": 11.0, "g1": -0.16, "g2": 0.7},
            {"type": "Exponential", "flux": 1.3e6, "half_light_radius": 1.9, "g1": 0.3, "g2": -0.1},
            {"type": "Exponential", "flux": 1.8e6, "half_light_radius": 5.0, "g1": 0.0, "g2": 0.14},
            {"type": "Exponential", "flux": 6.6e6, "half_light_radius": 4.8, "g1": 0.26, "g2": 0.5},
            {"type": "Exponential", "flux": 7.0e5, "half_light_radius": 3.1, "g1": -0.3, "g2": 0.0},
            {"type": "DeVaucouleurs", "flux": 1.6e5, "half_light_radius": 3.5, "g1": 0.2, "g2": 0.4},
            {"type": "DeVaucouleurs", "flux": 2.0e5, "half_light_radius": 1.6, "g1": -0.06, "g2": -0.2},
            {"type": "DeVaucouleurs", "flux": 8.3e5, "half_light_radius": 5.1, "g1": 0.29, "g2": 0.0},
            {"type": "DeVaucouleurs", "flux": 4.5e5, "half_light_radius": 2.5, "g1": 0.4, "g2": 0.3},
            {"type": "DeVaucouleurs", "flux": 6.2e5, "half_light_radius": 4.9, "g1": -0.08, "g2": -0.01},
            {"type": "Gaussian", "flux": 4.7e6, "half_light_radius": 2.5, "g1": 0.07, "g2": -0.35},
            {"type": "Gaussian", "flux": 5.8e6, "half_light_radius": 3.1, "g1": 0.03, "g2": 0.4},
            {"type": "Gaussian", "flux": 2.3e5, "half_light_radius": 0.5, "g1": 0.0, "g2": 0.0},
            {"type": "Gaussian", "flux": 1.6e6, "half_light_radius": 3.0, "g1": 0.18, "g2": -0.29},
            {"type": "Gaussian", "flux": 3.5e5, "half_light_radius": 4.6, "g1": 0.5, "g2": 0.35},
            {"type": "Gaussian", "flux": 5.9e5, "half_light_radius": 9.5, "g1": 0.1, "g2": 0.55},
            {"type": "Gaussian", "flux": 4.0e5, "half_light_radius": 1.0, "g1": 0.0, "g2": 0.0},
        ]

        # Mapping of profile types to their galsim constructors.
        profile_constructors = {
            "Sersic": galsim.Sersic,
            "Exponential": galsim.Exponential,
            "DeVaucouleurs": galsim.DeVaucouleurs,
            "Gaussian": galsim.Gaussian,
        }

        # Create a galsim image to draw the sources onto. The exposure image
        # that is passed to this method will be modified in place.
        image = galsim.ImageF(exposure.image.array, bounds=image_bounds)

        # Generate random positions within exposure bounds, avoiding edges by a
        # margin.
        margin_x, margin_y = 0.05 * exp_bbox.width, 0.05 * exp_bbox.height
        self.positions = np.random.uniform(
            [exp_bbox.minX + margin_x, exp_bbox.minY + margin_y],
            [exp_bbox.maxX - margin_x, exp_bbox.maxY - margin_y],
            (len(source_params), 2),
        ).tolist()

        # Loop over the sources and draw them onto the image cutout by cutout.
        for i, params in enumerate(source_params):
            # Dynamically get constructor and remove type from params.
            constructor = profile_constructors[params.pop("type")]

            # Get shear parameters and remove them from params.
            g1, g2 = params.pop("g1"), params.pop("g2")

            # The extent of the cutout should be large enough to contain the
            # entire object above the background level. Some empirical factor
            # is used to mitigate artifacts.
            half_extent = 10 * params["half_light_radius"] * (1 + 2 * np.sqrt(g1**2 + g2**2))

            # Pass the remaining params to the constructor and apply shear.
            galsim_object = constructor(**params).shear(galsim.Shear(g1=g1, g2=g2))

            # Retrieve the position of the object.
            x, y = self.positions[i]
            pos = galsim.PositionD(x, y)

            # Get the bounds of the sub-image based on the object position.
            sub_image_bounds = galsim.BoundsI(
                *map(int, [x - half_extent, x + half_extent, y - half_extent, y + half_extent])
            )

            # Identify the overlap region, which could be partially outside the
            # image bounds.
            sub_image_bounds = sub_image_bounds & image.bounds

            # Check that there is some overlap.
            assert sub_image_bounds.isDefined(), "No overlap with image bounds"

            # Get the sub-image cutout.
            sub_image = image[sub_image_bounds]

            # Draw the object onto the image within the the sub-image bounds.
            galsim_object.drawImage(
                image=sub_image,
                offset=pos - sub_image.true_center,
                method="real_space",  # It saves memory, usable w/o convolution
                add_to_image=True,  # Add flux to existing image
                scale=pixel_scale,
            )

        # Add a constant background to the entire image (both in e-/pixel).
        image += self.background

        # Add noise to the image which is in electrons. Note that we won't
        # specify a `sky_level` here to avoid double-counting it, as it's
        # already included as the background.
        image.addNoise(galsim.PoissonNoise(rng))

        # Subtract off the background to get the sky-subtracted image.
        image -= self.background

        # Adjust each amplifier's image segment by its respective gain. After
        # this step, the image will be in ADUs.
        for bounds, gain in zip(amp_bounds_list, self.amp_gain_list):
            image[bounds] /= gain

        # We know that the exposure has already been modified in place, but
        # just to be extra sure, we'll set the exposure image explicitly.
        exposure.image.array = image.array

        # Create a variance plane for the exposure while including signal as a
        # pollutant. Note that the exposure image is pre-adjusted for gain,
        # unlike 'self.background'. Thus, we divide the background by the
        # corresponding gain before adding it to the image. This leads to the
        # variance plane being in units of ADU^2.
        for bbox, gain in zip(self.amp_bbox_list, self.amp_gain_list):
            exposure.variance[bbox].array = (exposure.image[bbox].array + self.background / gain) / gain

        return exposure

    @methodParametersProduct(
        average_gain=[1.4, 1.7],
        predefined_gain_type=["average", "per-amp", None],
        gain_sigma_factor=[0, 0.008],
        sky_level=[2e6, 4e6],
        average_across_amps=[False, True],
        in_place=[False, True],
    )
    def test_variance_signal_removal(
        self, average_gain, predefined_gain_type, gain_sigma_factor, sky_level, average_across_amps, in_place
    ):
        exposure = self.buildExposure(
            average_gain=average_gain, gain_sigma_factor=gain_sigma_factor, sky_level=sky_level
        )

        # Save the original variance plane for comparison, assuming it has
        # Poisson contribution from the source signal.
        signal_polluted_variance = exposure.variance.clone()

        # Check that the variance plane has no negative values.
        self.assertTrue(
            np.all(signal_polluted_variance.array >= 0),
            "Variance plane has negative values (before correction)",
        )

        if predefined_gain_type == "average":
            predefined_gain = average_gain
            predefined_gains = None
        elif predefined_gain_type == "per-amp":
            predefined_gain = None
            predefined_gains = self.amp_gain_list
        elif predefined_gain_type is None:
            # Allow the 'remove_signal_from_variance' function to estimate the
            # gain itself before it attempts to remove the signal from the
            # variance plane.
            predefined_gain = None
            predefined_gains = None

        # Set the relative tolerance for the variance plane checks.
        if predefined_gain_type == "average" or (predefined_gain_type is None and average_across_amps):
            # Relax the tolerance if we are simply averaging across amps to
            # roughly estimate the overall gain.
            rtol = 0.01
            estimate_average_gain = True
        else:
            # Tighten tolerance for the 'predefined_gain_type' of 'per-amp' or
            # for a more accurate per-amp gain estimation strategy.
            rtol = 2e-7
            estimate_average_gain = False

        # Remove the signal from the variance plane.
        signal_free_variance = remove_signal_from_variance(
            exposure,
            gain=predefined_gain,
            gains=predefined_gains,
            average_across_amps=average_across_amps,
            in_place=in_place,
        )

        # Retrieve the variance plane from the exposure in case it was modified
        # in place.
        if in_place:
            assert signal_free_variance is None, "Expected in-place modification but got a return value"
            signal_free_variance = exposure.variance

        # Check that the variance plane has been modified.
        self.assertFloatsNotEqual(signal_polluted_variance.array, signal_free_variance.array)

        # Check that the corrected variance plane has no negative values.
        self.assertTrue(
            np.all(signal_free_variance.array >= 0), "Variance plane has negative values (after correction)"
        )

        for bbox, gain in zip(self.amp_bbox_list, self.amp_gain_list):
            # Calculate the true variance in theoretical terms.
            true_var_amp = self.background / gain**2
            # Pair each variance with the appropriate context manager before
            # looping through them.
            var_context_pairs = [
                # For the signal-free variance, directly execute the checks.
                (signal_free_variance, nullcontext()),
                # For the signal-polluted variance, expect AssertionError
                # unless we are averaging across amps.
                (
                    signal_polluted_variance,
                    nullcontext() if estimate_average_gain else self.assertRaises(AssertionError),
                ),
            ]
            for var, context_manager in var_context_pairs:
                # Extract the segment of the variance plane for the amplifier.
                var_amp = var[bbox]
                with context_manager:
                    if var is signal_polluted_variance and estimate_average_gain:
                        # Skip rigorous checks on the signal-polluted variance,
                        # if we are averaging across amps.
                        pass
                    else:
                        # Get the variance value at the first pixel of the
                        # segment to compare with the rest of the pixels and
                        # the true variance.
                        v00 = var_amp.array[0, 0]
                        # Assert that the variance plane is almost uniform
                        # across the segment because the signal has been
                        # removed from it and the background is constant.
                        self.assertFloatsAlmostEqual(var_amp.array, v00, rtol=rtol)
                        # Assert that the variance plane is almost equal to the
                        # true variance across the segment.
                        self.assertFloatsAlmostEqual(v00, true_var_amp, rtol=rtol)

        if (
            SAVE_PLOT
            and not in_place
            and not average_across_amps
            and gain_sigma_factor in (0, 0.008)
            and sky_level == 4e6
            and average_gain == 1.7
            and predefined_gain_type is None
        ):
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8.5))
            plt.subplots_adjust(wspace=0.17, hspace=0.17)
            colorbar_aspect = 12

            amp_background_variance_ADU_list = [self.background / gain**2 for gain in self.amp_gain_list]
            amp_background_image_ADU_list = [self.background / gain for gain in self.amp_gain_list]
            # Calculate the mean value that corresponds to the background for
            # the variance plane, adjusting for the gain.
            background_mean_variance_ADU = np.mean(
                [self.background / gain**2 for gain in self.amp_gain_list]
            )

            # Extract the variance planes and the image from the exposure.
            arr1 = signal_polluted_variance.array  # Variance with signal
            arr2 = signal_free_variance.array  # Variance without signal
            exp_im = exposure.image.clone()  # Clone of the image plane

            # Incorporate the gain-adjusted background into the image plane to
            # enable combined visualization of sources with the background.
            for gain, bbox in zip(self.amp_gain_list, self.amp_bbox_list):
                exp_im[bbox].array += self.background / gain
            arr3 = exp_im.array

            # Define colors visually distinct from each other for the subplots.
            original_variance_color = "#8A2BE2"  # Periwinkle
            corrected_variance_color = "#618B3C"  # Lush Forest Green
            sky_variance_color = "#c3423f"  # Crimson Red
            amp_colors = [
                "#1f77b4",  # Muted Blue
                "#ff7f0e",  # Vivid Orange
                "#2ca02c",  # Kelly Green
                "#d62728",  # Brick Red
                "#9467bd",  # Soft Purple
                "#8B4513",  # Saddle Brown
                "#e377c2",  # Pale Violet Red
                "#202020",  # Onyx
            ]
            arrowheads_lr = ["$\u25C0$", "$\u25B6$"]  # Left- & right-pointing
            arrowheads_ud = ["$\u25B2$", "$\u25BC$"]  # Up- & down-pointing

            # Set titles for the subplots.
            ax1.set_title("Original variance plane", color=original_variance_color)
            ax2.set_title("Corrected variance plane", color=corrected_variance_color)
            ax3.set_title("Image + background ($\\mathit{uniform}$)")
            ax4.set_title("Histogram of variances")

            # Collect all vertical and horizontal line positions to find the
            # amp boundaries.
            vlines, hlines = set(), set()
            for bbox in self.amp_bbox_list:
                # Adjst by 0.5 for merging of lines at the boundaries.
                vlines.update({bbox.minX - 0.5, bbox.maxX + 0.5})
                hlines.update({bbox.minY - 0.5, bbox.maxY + 0.5})

            # Filter lines at the edges of the overall image bbox.
            image_bbox = exposure.getBBox()
            vlines = {x for x in vlines if image_bbox.minX < x < image_bbox.maxX}
            hlines = {y for y in hlines if image_bbox.minY < y < image_bbox.maxY}

            # Plot image and variance planes.
            for plane, arr, ax in zip(
                ("variance", "variance_corrected", "image"), (arr1, arr2, arr3), (ax1, ax2, ax3)
            ):
                # We skip 'variance_corrected' in the loop below because we use
                # the same normalization and colormap as 'variance' for it.
                if plane in ["variance", "image"]:
                    # Get the normalization.
                    vmin, vmax = arr.min(), arr.max()
                    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

                    # Get the thresholds corresponding to per-amp backgrounds
                    # and their positions in the normalized color scale.
                    thresholds = (
                        amp_background_variance_ADU_list
                        if plane.startswith("variance")
                        else amp_background_image_ADU_list
                    )
                    threshold_positions = [norm(t) for t in thresholds]
                    threshold = np.mean(thresholds)
                    threshold_position = np.mean(threshold_positions)

                    # Create a custom colormap with two distinct colors for the
                    # sky and source contributions.
                    border = (threshold - vmin) / (vmax - vmin)
                    colors1 = plt.cm.Purples_r(np.linspace(0, 1, int(border * 256)))
                    colors2 = plt.cm.Greens(np.linspace(0, 1, int((1 - border) * 256)))
                    colors = np.vstack((colors1, colors2))
                    cmap = mcolors.LinearSegmentedColormap.from_list("cmap", colors)

                # Plot the array with the custom colormap and normalization.
                im = ax.imshow(arr, cmap=cmap, norm=norm)

                # Add colorbars to the plot.
                cbar = fig.colorbar(im, aspect=colorbar_aspect, pad=0)

                # Change the number of ticks on the colorbar for better
                # spacing. Needs to be done before modifying the tick labels.
                cbar.ax.locator_params(nbins=7)

                # Enhance readability by scaling down colorbar tick labels.
                unit = "ADU$^2$" if plane.startswith("variance") else "ADU"
                adjust_tick_scale(cbar.ax, {"y": f"Value [{{scale}} {unit}]"})

                # Mark per-amp thresholds with dotted lines on the colorbar.
                for tp in [min(thresholds), max(thresholds)]:
                    cbar.ax.axhline(tp, color="white", linestyle="-", linewidth=1.5, alpha=0.4)
                    cbar.ax.axhline(tp, color=sky_variance_color, linestyle=":", linewidth=1.5, alpha=0.9)

                # Mark mean threshold with facing arrowheads on the colorbar.
                cbar.ax.annotate(
                    arrowheads_lr[1],  # Right-pointing arrowhead
                    xy=(0, threshold_position),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    xytext=(0, 0),
                    ha="left",
                    va="center",
                    fontsize=6,
                    color=sky_variance_color,
                    clip_on=False,
                    alpha=0.9,
                )
                cbar.ax.annotate(
                    arrowheads_lr[0],  # Left-pointing arrowhead
                    xy=(1, threshold_position),
                    xycoords="axes fraction",
                    textcoords="offset points",
                    xytext=(0, 0),
                    ha="right",
                    va="center",
                    fontsize=6,
                    color=sky_variance_color,
                    clip_on=False,
                    alpha=0.9,
                )

                # Add text inside the colorbar to label the average threshold
                # position.
                sky_level_text = "$\u27E8$" + "Sky" + "$\u27E9$"  # <Sky>
                sky_level_text_artist = cbar.ax.text(
                    0.5,
                    threshold_position,
                    sky_level_text,
                    va="center",
                    ha="center",
                    transform=cbar.ax.transAxes,
                    fontsize=8,
                    color=sky_variance_color,
                    rotation="vertical",
                    alpha=0.9,
                    path_effects=outline_effect(2),
                )

                # Setup renderer and transformation.
                renderer = fig.canvas.get_renderer()
                transform = cbar.ax.transAxes.inverted()

                # Transform the bounding box and calculate adjustment for
                # 'sky_level_text_artist' for when it goes beyond the colorbar.
                sky_level_text_bbox = sky_level_text_artist.get_window_extent(renderer).transformed(transform)
                adjustment = 1.4 * sky_level_text_bbox.height / 2

                if sky_level_text_bbox.ymin < 0:
                    sky_level_text_artist.set_y(adjustment)
                elif sky_level_text_bbox.ymax > 1:
                    sky_level_text_artist.set_y(1 - adjustment)

                # Draw amp boundaries as vertical and/or horizontal lines.
                line_color = "white" if np.mean(norm(arr)) > 0.5 else "#808080"
                for x in vlines:
                    ax.axvline(x=x, color=line_color, linestyle="--", linewidth=1, alpha=0.7)
                for y in hlines:
                    ax.axhline(y=y, color=line_color, linestyle="--", linewidth=1, alpha=0.7)
                # Hide all x and y tick marks.
                ax.tick_params(axis="both", which="both", bottom=False, top=False, left=False, right=False)
                # Hide all x and y tick labels.
                ax.set_xticklabels([])
                ax.set_yticklabels([])

            # Additional ax2 annotations:
            # Labels amplifiers with their respective gains for a visual check.
            for bbox, name, gain, color in zip(
                self.amp_bbox_list, self.amp_names, self.amp_gain_list, amp_colors
            ):
                # Get the center of the bbox to label the gain value.
                bbox_center = (bbox.minX + bbox.maxX) / 2, (bbox.minY + bbox.maxY) / 2
                # Label the gain value at the center of each amplifier segment.
                ax2.text(
                    *bbox_center,
                    f"gain$_{{\\rm \\, {name} \\,}}$: {gain:.3f}",
                    fontsize=9,
                    color=color,
                    alpha=0.95,
                    ha="center",
                    va="center",
                    path_effects=outline_effect(2),
                )

            # Additional ax3 annotations:
            # Label sources with numbers on the image plane.
            for i, pos in enumerate(self.positions, start=1):
                ax3.text(
                    *pos,
                    f"{i}",
                    fontsize=7,
                    color=sky_variance_color,
                    path_effects=outline_effect(1.5),
                    alpha=0.9,
                )

            # Now we use ax4 to plot the histograms of the variance planes for
            # comparison.
            # Plot the histogram of the original variance plane.
            hist_values, bins, _ = ax4.hist(
                arr1.flatten(),
                bins=80,
                histtype="step",
                color=original_variance_color,
                alpha=0.9,
                label="Original variance",
            )
            # Fill the area under the step.
            ax4.fill_between(
                bins[:-1],
                hist_values,
                step="post",
                color=original_variance_color,
                alpha=0.09,
                hatch="/////",
                label=" ",
            )
            # Plot the histogram of the corrected variance plane.
            ax4.hist(
                arr2.flatten(),
                bins=80,
                histtype="bar",
                color=corrected_variance_color,
                alpha=0.9,
                label="Corrected variance",
            )
            adjust_tick_scale(ax4, {"x": "Variance [{scale} ADU$^2$]", "y": "Number of pixels / {scale}"})
            ax4.yaxis.set_label_position("right")
            ax4.yaxis.tick_right()
            ax4.axvline(
                background_mean_variance_ADU,
                color=sky_variance_color,
                linestyle="--",
                linewidth=1,
                alpha=0.9,
                label="Average sky variance\nacross all amps",
            )

            #
            sorted_vars = sorted(amp_background_variance_ADU_list)
            for i, (x, name, gain, color) in enumerate(
                zip(amp_background_variance_ADU_list, self.amp_names, self.amp_gain_list, amp_colors)
            ):
                arrowhead = arrowheads_ud[int(gain < average_gain)]
                arrowhead_text = ax4.annotate(
                    arrowhead,
                    xy=(x, 0),
                    xycoords=("data", "axes fraction"),
                    textcoords="offset points",
                    xytext=(0, 0),
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                    color=color,
                    clip_on=False,
                    alpha=0.85,
                    path_effects=outline_effect(1.5),
                )
                if i == 0:
                    # Draw the canvas once to make sure the renderer is active.
                    fig.canvas.draw()
                    # Get the bounding box of the text annotation in axes
                    # fraction.
                    bbox_axes = arrowhead_text.get_window_extent().transformed(ax4.transAxes.inverted())
                    # Get the height of the text annotation in axes fraction.
                    height = bbox_axes.height
                # Alternate the arrowhead y positions to minimize overlap.
                if np.where(sorted_vars == x)[0] % 2 != 0:
                    arrowhead_text.xy = (x, height)

                # Create a proxy artist for the legend since annotations are
                # not shown in the legend.
                label = "True variance of" if i == 0 else "$\u21AA$"
                ax4.scatter(
                    [],
                    [],
                    color=color,
                    marker=arrowhead,
                    s=15,
                    label=f"{label} {name}",
                    alpha=0.85,
                    path_effects=outline_effect(1.5),
                )

            # Group the legend handles and label them.
            adjust_legend_with_groups(
                ax4, [(0, 1), 2, 3, *range(4, 4 + len(self.amp_names))], colors="match", handlelength=1.9
            )

            # Align the histogram (bottom right panel) with the colorbar of the
            # corrected variance plane (top right panel) for aesthetic reasons.
            pos2 = ax2.get_position()
            pos4 = ax4.get_position()
            fig.canvas.draw()  # Render to ensure accurate colorbar width
            cbar_width = cbar.ax.get_position().width
            ax4.set_position([pos2.x0, pos4.y0, pos2.width + cbar_width, pos4.height])

            # Increase all axes spines' linewidth by 20% for a bolder look.
            for ax in fig.get_axes():
                for spine in ax.spines.values():
                    spine.set_linewidth(spine.get_linewidth() * 1.2)

            # Save the figure.
            filename = f"variance_plane_gain{average_gain}_sigma{gain_sigma_factor}_sky{sky_level}.png"
            fig.savefig(filename, dpi=300)
            print(f"Saved plot of variance plane before and after correction in {filename}")


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
