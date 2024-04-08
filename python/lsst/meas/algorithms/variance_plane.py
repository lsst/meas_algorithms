# This file is part of meas_algorithms.
#
# LSST Data Management System
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See COPYRIGHT file at the top of the source tree.
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
"""Utility functions related to the variance plane of Exposure objects. Tested
in `ip_isr/tests/test_variance_plane.py` to avoid circular dependencies.
"""

import numpy as np

__all__ = ["remove_signal_from_variance"]


def remove_signal_from_variance(exposure, gain=None, gains=None, average_across_amps=False, in_place=False):
    """
    Removes the Poisson contribution from actual sources in the variance plane
    of an Exposure.

    If neither gain nor gains are provided, the function estimates the gain(s).
    If ``average_across_amps`` is True, a single gain value for the entire
    image is estimated. If False, individual gain values for each amplifier are
    estimated. The estimation involves a linear fit of variance versus image
    plane.

    Parameters
    ----------
    exposure : `~lsst.afw.image.Exposure`
        The background-subtracted exposure containing a variance plane to be
        corrected for source contributions.
    gain : `float`, optional
        The gain value for the entire image. This parameter is used if
        ``gains`` is not provided. If both ``gain`` and ``gains`` are None, and
        ``average_across_amps`` is True, ``gain`` is estimated from the image
         and variance planes.
    gains : dict[`str`, `float`], optional
        A dictionary mapping amplifier ID (as a string) to gain value. This
        parameter is used if ``gain`` is not provided. If both ``gain`` and
        ``gains`` are None, and ``average_across_amps`` is False, ``gains`` are
        estimated from the image and variance planes.
    average_across_amps : `bool`, optional
        Determines the gain estimation strategy. If True, the gain for the
        entire image is estimated at once. If False, individual gains for each
        amplifier are estimated. This parameter is ignored if either ``gain``
        or ``gains`` is specified.
    in_place : `bool`, optional
        If True, the variance plane of the input Exposure is modified in place.
        A modified copy of the variance plane is always returned irrespective
        of this.

    Returns
    -------
    variance_plane : `~lsst.afw.image.Image`
        The corrected variance plane, with the signal contribution removed.

    Raises
    ------
    AttributeError
        If amplifiers cannot be retrieved from the exposure.
    ValueError
        If both ``gain`` and ``gains`` are provided, or if the number of
        provided ``gains`` does not match the number of amplifiers.
    """
    variance_plane = exposure.variance if in_place else exposure.variance.clone()
    if gain is None and gains is None:
        if average_across_amps:
            amp_bboxes = [exposure.getBBox()]
        else:
            try:
                amps = exposure.getDetector().getAmplifiers()
                amp_bboxes = [amp.getBBox() for amp in amps]
            except AttributeError:
                raise AttributeError(
                    "Could not retrieve amplifiers from exposure. To compute a simple gain value across the "
                    "entire image, use average_across_amps=True."
                )
        # Fit a straight line to variance vs (sky-subtracted) signal. Then
        # evaluate that line at zero signal to get an estimate of the
        # signal-free variance.
        for amp_bbox in amp_bboxes:
            amp_im_arr = exposure[amp_bbox].image.array
            amp_var_arr = variance_plane[amp_bbox].array
            good = (amp_var_arr != 0) & np.isfinite(amp_var_arr) & np.isfinite(amp_im_arr)
            fit = np.polyfit(amp_im_arr[good], amp_var_arr[good], deg=1)
            # Fit is [1/gain, sky_var].
            amp_gain = 1.0 / fit[0]
            variance_plane[amp_bbox].array[good] -= amp_im_arr[good] / amp_gain
    elif gain is None and gains is not None:
        amps = exposure.getDetector().getAmplifiers()
        namps = len(amps)
        if len(gains) != namps:
            raise ValueError(
                f"Incorrect number of gains provided: {len(gains)} values for {namps} amplifiers."
            )
        for amp in amps:
            amp_bbox = amp.getBBox()
            amp_gain = gains[amp.getName()]
            im_arr = exposure[amp_bbox].image.array
            variance_plane[amp_bbox].array -= im_arr / amp_gain
    elif gain is not None and gains is None:
        im_arr = exposure.image.array
        variance_plane.array -= im_arr / gain
    elif gain is not None and gains is not None:
        raise ValueError(
            "Both 'gain' and 'gains' are provided. Please provide only one of them or none at "
            "all in case of automatic gain estimation from the image and variance planes."
        )
    # Check that the variance plane has no negative values.
    if np.any(variance_plane.array < 0):
        raise ValueError("Corrected variance plane has negative values.")
    return variance_plane
