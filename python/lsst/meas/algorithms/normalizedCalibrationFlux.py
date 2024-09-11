# This file is part of lsst.meas.algorithms.
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

__all__ = ["NormalizedCalibrationFluxConfig", "NormalizedCalibrationFluxTask",
           "NormalizedCalibrationFluxError"]

import numpy as np

from lsst.afw.image import ApCorrMap
from lsst.afw.math import ChebyshevBoundedField
import lsst.pex.config
import lsst.pipe.base
from .measureApCorr import MeasureApCorrTask, MeasureApCorrError
from .sourceSelector import sourceSelectorRegistry


class NormalizedCalibrationFluxError(lsst.pipe.base.AlgorithmError):
    """Raised if Aperture Correction fails in a non-recoverable way.

    Parameters
    ----------
    n_initial_sources : `int`
        Number of sources selected by the fallback source selector.
    n_calib_flux_flag : `int`
        Number of selected sources with raw calibration flux flag unset.
    n_ref_flux_flag : `int`
        Number of selected sources with reference flux flag unset.
    """
    def __init__(self, *, n_initial_sources, n_calib_flux_flag, n_ref_flux_flag):
        msg = "There are no valid stars to compute normalized calibration fluxes."
        msg += (f" Of {n_initial_sources} initially selected sources, {n_calib_flux_flag} have good raw"
                f" calibration fluxes and {n_ref_flux_flag} have good reference fluxes.")
        super().__init__(msg)
        self.n_initial_sources = n_initial_sources
        self.n_calib_flux_flag = n_calib_flux_flag
        self.n_ref_flux_flag = n_ref_flux_flag

    @property
    def metadata(self):
        metadata = {"n_init_sources": self.n_initial_sources,
                    "n_calib_flux_flag": self.n_calib_flux_flag,
                    "n_ref_flux_flag": self.n_ref_flux_flag}
        return metadata


class NormalizedCalibrationFluxConfig(lsst.pex.config.Config):
    """Configuration parameters for NormalizedCalibrationFluxTask.
    """
    measure_ap_corr = lsst.pex.config.ConfigurableField(
        target=MeasureApCorrTask,
        doc="Subtask to measure aperture corrections.",
    )
    raw_calibflux_name = lsst.pex.config.Field(
        doc="Name of raw calibration flux to normalize.",
        dtype=str,
        default="base_CompensatedTophatFlux_12",
    )
    normalized_calibflux_name = lsst.pex.config.Field(
        doc="Name of normalized calibration flux.",
        dtype=str,
        default="base_NormalizedCompensatedTophatFlux",
    )
    do_set_calib_slot = lsst.pex.config.Field(
        doc="Set the calib flux slot to the normalized flux?",
        dtype=bool,
        default=True,
    )
    do_measure_ap_corr = lsst.pex.config.Field(
        doc="Measure the aperture correction? (Otherwise, just apply.)",
        dtype=bool,
        default=True,
    )
    fallback_source_selector = sourceSelectorRegistry.makeField(
        doc="Selector that is used as a fallback if the full aperture correction "
            "fails.",
        default="science",
    )

    def setDefaults(self):
        super().setDefaults()

        self.measure_ap_corr.refFluxName = "base_CircularApertureFlux_12_0"

        # This task is meant to be used early when we focus on PSF stars.
        selector = self.measure_ap_corr.sourceSelector["science"]
        selector.doUnresolved = False
        selector.flags.good = ["calib_psf_used"]
        selector.flags.bad = []
        selector.signalToNoise.fluxField = self.raw_calibflux_name + "_instFlux"
        selector.signalToNoise.errField = self.raw_calibflux_name + "_instFluxErr"
        # Do median for this.
        self.measure_ap_corr.fitConfig.orderX = 0
        self.measure_ap_corr.fitConfig.orderY = 0

        fallback_selector = self.fallback_source_selector["science"]
        fallback_selector.doFluxLimit = False
        fallback_selector.doFlags = True
        fallback_selector.doUnresolved = False
        fallback_selector.doSignalToNoise = False
        fallback_selector.doIsolated = False
        fallback_selector.flags.good = ["calib_psf_used"]
        fallback_selector.flags.bad = []


class NormalizedCalibrationFluxTask(lsst.pipe.base.Task):
    """Task to measure the normalized calibration flux.

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Schema for the input table; will be modified in place.
    **kwargs : `dict`
        Additional kwargs to pass to lsst.pipe.base.Task.__init__()

    Raises
    ------
    NormalizedCalibrationFluxError
        Raised if there are not enough sources to calculate normalization.
    """
    ConfigClass = NormalizedCalibrationFluxConfig
    _DefaultName = "normalizedCalibrationFlux"

    def __init__(self, schema, **kwargs):
        lsst.pipe.base.Task.__init__(self, **kwargs)

        if self.config.do_measure_ap_corr:
            self.makeSubtask(
                "measure_ap_corr",
                schema=schema,
                namesToCorrect=[self.config.raw_calibflux_name],
            )

        name = self.config.normalized_calibflux_name
        self.flux_name = name + "_instFlux"
        if self.flux_name not in schema:
            schema.addField(
                self.flux_name,
                type=float,
                doc=f"Normalized calibration flux from {self.config.raw_calibflux_name}.",
            )
        self.err_name = name + "_instFluxErr"
        if self.err_name not in schema:
            schema.addField(
                self.err_name,
                type=float,
                doc=f"Normalized calibration flux error from {self.config.raw_calibflux_name}.",
            )
        self.flag_name = name + "_flag"
        if self.flag_name not in schema:
            schema.addField(
                self.flag_name,
                type="Flag",
                doc=f"Normalized calibration flux failure flag from {self.config.raw_calibflux_name}.",
            )

        if self.config.do_set_calib_slot:
            schema.getAliasMap().set("slot_CalibFlux", name)

        self.makeSubtask("fallback_source_selector")

        self.schema = schema

    def run(self, *, exposure, catalog):
        """Measure the Normalized calibration flux.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure the normalized calibration flux is measured on.
        catalog : `lsst.afw.table.SourceCatalog`
            SourceCatalog containing measurements to be used to compute
            normalized calibration fluxes.  The catalog is modified in-place.

        Returns
        -------
        Struct : `lsst.pipe.base.Struct`
            Contains the following:

            ``ap_corr_map``
                aperture correction map (`lsst.afw.image.ApCorrMap`)
                that contains two entries for the raw flux field:
                - flux field (e.g. config.{raw_calibflux_name}_instFlux): 2d model
                - flux sigma field (e.g. config.{raw_calibflux_name}_instFluxErr): 0 field
        """
        self.log.info("Measuring normalized calibration flux from %s", self.config.raw_calibflux_name)

        raw_name = self.config.raw_calibflux_name
        raw_flux_name = raw_name + "_instFlux"
        raw_fluxerr_name = raw_name + "_instFluxErr"
        norm_name = self.config.normalized_calibflux_name

        if self.config.do_measure_ap_corr:
            ap_corr_field, ap_corr_err_field = self._measure_aperture_correction(exposure, catalog)
        else:
            use_identity = False
            ap_corr_map = exposure.info.getApCorrMap()
            if ap_corr_map is None:
                self.log.warning(
                    "Exposure does not have a valid normalization map; using identity normalization.",
                )
                use_identity = True
            else:
                ap_corr_field = ap_corr_map.get(raw_flux_name)
                ap_corr_err_field = ap_corr_map.get(raw_fluxerr_name)
                if not ap_corr_field or not ap_corr_err_field:
                    self.log.warning(
                        "Exposure aperture correction map is missing %s/%s for normalization; "
                        "using identity normalization.",
                        raw_flux_name,
                        raw_fluxerr_name,
                    )
                    use_identity = True

            if use_identity:
                ap_corr_field = ChebyshevBoundedField(exposure.getBBox(), np.array([[1.0]]))
                ap_corr_err_field = ChebyshevBoundedField(exposure.getBBox(), np.array([[0.0]]))

        corrections = ap_corr_field.evaluate(
            catalog["slot_Centroid_x"],
            catalog["slot_Centroid_y"],
        )

        input_flux_name = raw_flux_name
        input_fluxerr_name = raw_fluxerr_name
        input_flag_name = raw_name + "_flag"
        output_flux_name = norm_name + "_instFlux"
        output_fluxerr_name = norm_name + "_instFluxErr"
        output_flag_name = norm_name + "_flag"

        if catalog.isContiguous():
            catalog[output_flux_name] = catalog[input_flux_name] * corrections
            catalog[output_fluxerr_name] = catalog[input_fluxerr_name] * corrections

            output_flag = catalog[input_flag_name].copy()
            output_flag[corrections <= 0.0] = True
            catalog[output_flag_name] = output_flag
        else:
            # If the catalog is not contiguous we must go row-by-row.
            for i, row in enumerate(catalog):
                row[output_flux_name] = row[input_flux_name] * corrections[i]
                row[output_fluxerr_name] = row[input_fluxerr_name] * corrections[i]

                if row[input_flag_name] or corrections[i] <= 0.0:
                    row[output_flag_name] = True

        ap_corr_map = ApCorrMap()
        ap_corr_map[raw_flux_name] = ap_corr_field
        ap_corr_map[raw_fluxerr_name] = ap_corr_err_field

        return lsst.pipe.base.Struct(
            ap_corr_map=ap_corr_map,
        )

    def _measure_aperture_correction(self, exposure, catalog):
        """Internal method to do the aperture correction measurement.

        This measures the aperture correction with the regular task,
        and if that fails does a fallback median estimate.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure the normalized calibration flux is measured on.
            This is only used for the bounding box.
        catalog : `lsst.afw.table.SourceCatalog`
            SourceCatalog containing measurements to be used to compute
            normalized calibration flux.

        Returns
        -------
        ap_corr_field : `lsst.afw.math.ChebyshevBoundedField`
            Aperture correction field to normalize the calibration flux.
        ap_corr_err_field : `lsst.afw.math.ChebyshevBoundedField`
            Aperture correction to adjust the calibration flux error.
        """
        raw_name = self.config.raw_calibflux_name

        try:
            ap_corr_map = self.measure_ap_corr.run(
                exposure=exposure,
                catalog=catalog,
            ).apCorrMap

            ap_corr_field = ap_corr_map.get(raw_name + "_instFlux")
        except MeasureApCorrError as e:
            self.log.warning("Failed to measure full aperture correction for %s with the following error %s",
                             raw_name, e)

            initSel = self.fallback_source_selector.run(catalog, exposure=exposure).selected
            sel = (initSel & ~catalog[self.config.raw_calibflux_name + "_flag"]
                   & ~catalog[self.config.measure_ap_corr.refFluxName + "_flag"])

            if (n_sel := sel.sum()) == 0:
                # This is a fatal error.
                raise NormalizedCalibrationFluxError(
                    n_initial_sources=initSel.sum(),
                    n_calib_flux_flag=(initSel & ~catalog[self.config.raw_calibflux_name + "_flag"]).sum(),
                    n_ref_flux_flag=(initSel
                                     & ~catalog[self.config.measure_ap_corr.refFluxName + "_flag"]).sum()
                )
            self.log.info("Measuring normalized flux correction with %d stars from fallback selector.",
                          n_sel)

            ratio = np.median(
                catalog[self.config.measure_ap_corr.refFluxName + "_instFlux"][sel]
                / catalog[self.config.raw_calibflux_name + "_instFlux"][sel]
            )

            ap_corr_field = ChebyshevBoundedField(
                exposure.getBBox(),
                np.array([[ratio]]),
            )

            if catalog.isContiguous():
                catalog["apcorr_" + raw_name + "_used"] = sel
            else:
                for i, row in enumerate(catalog):
                    row["apcorr_" + raw_name + "_used"] = sel[i]

        # We are always setting the error field to 0, because we do not
        # have a good model for aperture correction uncertainties.
        ap_corr_err_field = ChebyshevBoundedField(
            exposure.getBBox(),
            np.array([[0.0]]),
        )

        return ap_corr_field, ap_corr_err_field
