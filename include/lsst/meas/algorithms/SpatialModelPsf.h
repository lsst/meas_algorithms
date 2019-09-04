// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_SPATIALMODELPSF_H)
#define LSST_MEAS_ALGORITHMS_SPATIALMODELPSF_H

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

/**
 * @file
 *
 * @brief Class used by SpatialCell for spatial PSF fittig
 *
 * @ingroup algorithms
 */
#include <memory>
#include <utility>
#include <vector>

#include "lsst/afw.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/SpatialCell.h"

namespace lsst {
namespace meas {
namespace algorithms {

template <typename PixelT>
std::pair<std::shared_ptr<afw::math::LinearCombinationKernel>, std::vector<double> >
createKernelFromPsfCandidates(afw::math::SpatialCellSet const& psfCells, geom::Extent2I const& dims,
                              geom::Point2I const& xy0, int const nEigenComponents, int const spatialOrder,
                              int const ksize, int const nStarPerCell = -1, bool const constantWeight = true,
                              int const border = 3);

template <typename PixelT>
int countPsfCandidates(afw::math::SpatialCellSet const& psfCells, int const nStarPerCell = -1);

template <typename PixelT>
std::pair<bool, double> fitSpatialKernelFromPsfCandidates(afw::math::Kernel* kernel,
                                                          afw::math::SpatialCellSet const& psfCells,
                                                          int const nStarPerCell = -1,
                                                          double const tolerance = 1e-5,
                                                          double const lambda = 0.0);
template <typename PixelT>
std::pair<bool, double> fitSpatialKernelFromPsfCandidates(
        afw::math::Kernel* kernel, afw::math::SpatialCellSet const& psfCells, bool const doNonLinearFit,
        int const nStarPerCell = -1, double const tolerance = 1e-5, double const lambda = 0.0);

template <typename ImageT>
double subtractPsf(afw::detection::Psf const& psf, ImageT* data, double x, double y,
                   double psfFlux = std::numeric_limits<double>::quiet_NaN());

template <typename Image>
std::pair<std::vector<double>, afw::math::KernelList> fitKernelParamsToImage(
        afw::math::LinearCombinationKernel const& kernel, Image const& image, geom::Point2D const& pos);

template <typename Image>
std::pair<std::shared_ptr<afw::math::Kernel>, std::pair<double, double> > fitKernelToImage(
        afw::math::LinearCombinationKernel const& kernel, Image const& image, geom::Point2D const& pos);

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif
