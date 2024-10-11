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


from lsst.pex.config import ChoiceField, Config, Field, FieldValidationError


class TreecorrConfig(Config):
    """A Config class that holds some of the parameters supported by treecorr.

    The fields in this class correspond to the parameters that can be passed to
    any calls to `treecorr` methods, including catalog creation and two-point
    correlation function calculations. The default values set for the fields
    are identical to the default values set in v4.2 of `treecorr`.

    A separate config class is used instead
    of constructing a `~lsst.pex.config.DictField` so that mixed types can be
    supported and the config can be validated at the beginning of the
    execution.

    Notes
    -----
    This is intended to be used in CalcRhoStatistics class. It only supports
    some of the fields that are relevant for rho-statistics calculations.
    """

    nbins = Field[int](
        doc=(
            "How many bins to use. "
            "(Exactly three of nbins, bin_size, min_sep, max_sep "
            "are required. If nbins is not given, it will be "
            "calculated from the values of the other three, "
            "rounding up to the next highest integer. "
            "In this case, bin_size will be readjusted to account "
            "for this rounding up."
        ),
        optional=True,
        check=lambda x: x > 0,
    )

    bin_size = Field[float](
        doc=(
            "The width of the bins in log(separation). "
            "Exactly three of nbins, bin_size, min_sep, max_sep are required. "
            "If bin_size is not given, it will be calculated from the values "
            "of the other three."
        ),
        optional=True,
    )

    min_sep = Field[float](
        doc=(
            "The minimum separation in units of sep_units, if relevant. "
            "Exactly three of nbins, bin_size, min_sep, max_sep are required. "
            "If min_sep is not given, it will be calculated from the values "
            "of the other three."
        ),
        optional=True,
    )

    max_sep = Field[float](
        doc=(
            "The maximum separation in units of sep_units, if relevant. "
            "Exactly three of nbins, bin_size, min_sep, max_sep are required. "
            "If max_sep is not given, it will be calculated from the values "
            "of the other three."
        ),
        optional=True,
    )

    sep_units = ChoiceField[str](
        doc=(
            "The units to use for the separation values, given as a string. "
            "This includes both min_sep and max_sep above, as well as the "
            "units of the output distance values."
        ),
        default="radian",
        optional=True,
        allowed={
            units: units for units in ["arcsec", "arcmin", "degree", "hour", "radian"]
        },
    )

    bin_slop = Field[float](
        doc=(
            "How much slop to allow in the placement of pairs in the bins. "
            "If bin_slop = 1, then the bin into which a particular pair is "
            "placed may be incorrect by at most 1.0 bin widths. "
            r"If None, use a bin_slop that gives a maximum error of 10% on "
            "any bin, which has been found to yield good results for most "
            "applications."
        ),
        default=None,
        optional=True,
    )

    precision = Field[int](
        doc=(
            "The precision to use for the output values. This specifies how many digits to write."
        ),
        default=4,
        optional=True,
        check=lambda x: x > 0,
    )

    metric = ChoiceField[str](
        doc=(
            "Which metric to use for distance measurements. For details, see "
            "https://rmjarvis.github.io/TreeCorr/_build/html/metric.html"
        ),
        default="Euclidean",
        optional=True,
        allowed={
            "Euclidean": "straight-line Euclidean distance between two points",
            "FisherRperp": (
                "the perpendicular component of the distance, "
                "following the definitions in "
                "Fisher et al, 1994 (MNRAS, 267, 927)"
            ),
            "OldRperp": (
                "the perpendicular component of the distance using the "
                "definition of Rperp from TreeCorr v3.x."
            ),
            "Rlens": (
                "Distance from the first object (taken to be a lens) to "
                "the line connecting Earth and the second object "
                "(taken to be a lensed source)."
            ),
            "Arc": "the true great circle distance for spherical coordinates.",
            "Periodic": "Like ``Euclidean``, but with periodic boundaries.",
        },
    )

    bin_type = ChoiceField[str](
        doc="What type of binning should be used?",
        default="Log",
        optional=True,
        allowed={
            "Log": (
                "Logarithmic binning in the distance. The bin steps will "
                "be uniform in log(r) from log(min_sep) .. log(max_sep)."
            ),
            "Linear": (
                "Linear binning in the distance. The bin steps will be "
                "uniform in r from min_sep .. max_sep."
            ),
            "TwoD": (
                "2-dimensional binning from x = (-max_sep .. max_sep) "
                "and y = (-max_sep .. max_sep). The bin steps will be "
                "uniform in both x and y. (i.e. linear in x,y)"
            ),
        },
    )

    var_method = ChoiceField[str](
        doc="Which method to use for estimating the variance",
        default="shot",
        optional=True,
        allowed={
            method: method
            for method in [
                "shot",
                "jackknife",
                "sample",
                "bootstrap",
                "marked_bootstrap",
            ]
        },
    )

    npatch = Field[int](
        doc="How many patches to split the catalog into for the purpose of "
        "jackknife variance or other options that involve running via "
        "patches (boostrap, marked_boostrap etc.)",
        default=1,
        optional=True,
    )

    num_bootstrap = Field[int](
        doc=(
            "How many bootstrap samples to use for the 'bootstrap' and 'marked_bootstrap' var methods."
        ),
        default=500,
        optional=True,
    )

    rng_seed = Field[int](
        doc="Value to seed the treecorr random number generator with. Used to generate patches.",
        default=13579,
    )

    def validate(self):
        # Docs inherited from base class
        super().validate()
        req_params = (self.nbins, self.bin_size, self.min_sep, self.max_sep)
        num_req_params = sum(param is not None for param in req_params)
        if num_req_params != 3:
            msg = (
                "You must specify exactly three of ``nbins``, ``bin_size``, ``min_sep`` and ``max_sep``"
                f" in treecorr_config. {num_req_params} parameters were set instead."
            )
            raise FieldValidationError(self.__class__.bin_size, self, msg)

        if self.min_sep is not None and self.max_sep is not None:
            if self.min_sep > self.max_sep:
                raise FieldValidationError(
                    self.__class__.min_sep, self, "min_sep must be <= max_sep"
                )
