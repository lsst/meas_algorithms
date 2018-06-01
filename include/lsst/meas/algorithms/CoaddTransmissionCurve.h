// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2018 LSST/AURA.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#ifndef LSST_MEAS_ALGORITHMS_CoaddTransmissionCurve_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_CoaddTransmissionCurve_h_INCLUDED

#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/TransmissionCurve.h"
#include "lsst/afw/table/Exposure.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 *  Create a TransmissionCurve that represents the effective throughput on a coadd.
 *
 *  @param[in]  coaddWcs       WCS that relates the coadd coordinate system to the sky.
 *  @param[in]  inputSensors   A catalog containing the WCSs, bounding boxes and polygons,
 *                             coaddition weights (in a field called 'weight'), and
 *                             TransmissionCurves of the sensor-level images that went
 *                             into the coadd.
 *
 *  @throws NotFoundError      Thrown if the 'weight' field does not exist in the schema.
 *
 *  @throws InvalidParameterError   Thrown if one or more inputs do not have a TransmissionCurve or
 *                                  a Wcs (ValidPolygons may be null to
 *                                  indicate no spatial restrictions other
 *                                  than the bounding box).
 *
 *  @returns a new TransmissionCurve object.
 */
std::shared_ptr<afw::image::TransmissionCurve const> makeCoaddTransmissionCurve(
        std::shared_ptr<afw::geom::SkyWcs const> coaddWcs, afw::table::ExposureCatalog const& inputSensors);

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_CoaddTransmissionCurve_h_INCLUDED
