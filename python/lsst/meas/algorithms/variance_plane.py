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
"""Utility functions related to the variance plane of Exposure objects.
"""

import numpy as np


__all__ = ['remove_signal_from_variance']


def remove_signal_from_variance(exposure, gain=None, average_across_amps=False, in_place=False):
    """Remove the Poisson contribution by actual sources from the variance
    plane of an Exposure.

    If no gain value is provided, it will be estimated as a linear fit of
    variance versus image plane. This estimation can be carried out on the
    whole image at once, or separately for each amplifier.

    Parameters
    -----------
    exposure : `afw.image.Exposure`
        Exposure that contains a variance plane that should be corrected for
        source contributions.
    gain : `float`, optional
        The gain value for the whole image. If not provided (the default),
        will be estimated from the image and variance planes.
    average_across_amps : `bool`, optional
        Whether the gain should be estimated on the whole image at once. If
        False (the default), a different gain value is estimated for each
        amplifier. Ignored if gain is not None.
    in_place : `bool`, optional
        If True, the variance plane is changed in place. Defaults to False.
    """
    variance_plane = exposure.variance if in_place else exposure.variance.clone()
    if average_across_amps:
        amp_bboxes = [exposure.getBBox()]
    else:
        try:
            amps = exposure.getDetector().getAmplifiers()
            amp_bboxes = [amp.getBBox() for amp in amps]
        except AttributeError:
            raise AttributeError("Could not retrieve amplifiers from exposure. To compute a simple gain "
                                 "value across the entire image, use average_across_amps=True.")
    if gain is None:
        # Fit a straight line to variance vs (sky-subtracted) signal.
        # The evaluate that line at zero signal to get an estimate of the
        # signal-free variance.
        for amp_bbox in amp_bboxes:
            amp_im_arr = exposure[amp_bbox].image.array
            amp_var_arr = variance_plane[amp_bbox].array
            good = (amp_var_arr != 0) & np.isfinite(amp_var_arr) & np.isfinite(amp_im_arr)
            fit = np.polyfit(amp_im_arr[good], amp_var_arr[good], deg=1)
            # fit is [1/gain, sky_var]
            gain = 1./fit[0]
            variance_plane[amp_bbox].array[good] -= amp_im_arr[good]/gain
    else:
        im_arr = exposure.image.array
        variance_plane.array -= im_arr/gain
    return variance_plane
