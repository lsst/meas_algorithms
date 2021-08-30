// -*- LSST-C++ -*-

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

/*!
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 *
 * @file
 *
 * @ingroup algorithms
 */
#include <cmath>
#include <memory>

#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/PcaPsf.h"
#include "lsst/afw/table/io/Persistable.cc"
#include "lsst/meas/algorithms/KernelPsfFactory.h"

namespace lsst {
namespace afw {
namespace table {
namespace io {

template std::shared_ptr<meas::algorithms::PcaPsf>
PersistableFacade<meas::algorithms::PcaPsf>::dynamicCast(std::shared_ptr<Persistable> const&);

}  // namespace io
}  // namespace table
}  // namespace afw
namespace meas {
namespace algorithms {

PcaPsf::PcaPsf(std::shared_ptr<afw::math::LinearCombinationKernel> kernel, geom::Point2D const& averagePosition)
        : KernelPsf(kernel, averagePosition) {
    if (!kernel) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "PcaPsf kernel must not be null");
    }
}

std::shared_ptr<afw::math::LinearCombinationKernel const> PcaPsf::getKernel() const {
    return std::static_pointer_cast<afw::math::LinearCombinationKernel const>(KernelPsf::getKernel());
}

std::shared_ptr<afw::detection::Psf> PcaPsf::clone() const { return std::make_shared<PcaPsf>(*this); }

std::shared_ptr<afw::detection::Psf> PcaPsf::resized(int width, int height) const {
    std::shared_ptr<afw::math::LinearCombinationKernel>
        kern = std::static_pointer_cast<afw::math::LinearCombinationKernel>(getKernel()->resized(width, height));
    return std::make_shared<PcaPsf>(kern, this->getAveragePosition());
}

namespace {

// registration for table persistence
KernelPsfFactory<PcaPsf, afw::math::LinearCombinationKernel> registration("PcaPsf");

}  // namespace

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

namespace lsst {
namespace afw {
namespace detection {

}
}  // namespace afw
}  // namespace lsst
