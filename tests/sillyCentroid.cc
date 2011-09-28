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
 
// -*- LSST-C++ -*-
/**
 * @file
 */

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/Measure.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that knows how to calculate centroids by guessing the wrong answer
 */
template<typename ExposureT>
class SillyAstrometer : public Algorithm<afwDet::Astrometry, ExposureT>
{
public:
    typedef Algorithm<afwDet::Astrometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<SillyAstrometer> Ptr;
    typedef boost::shared_ptr<SillyAstrometer const> ConstPtr;

    /// Ctor
    SillyAstrometer() : AlgorithmT() {}

    virtual std::string getName() const { return "SILLY"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<SillyAstrometer<ExposureT> >();
    }

    virtual PTR(afwDet::Astrometry) measureNull(void) const {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<afwDet::Astrometry>(NaN, NaN, NaN, NaN);
    }

    virtual void configure(lsst::pex::policy::Policy const& policy) {}

    virtual PTR(afwDet::Astrometry) measureOne(ExposurePatch<ExposureT> const&, afwDet::Source const&) const;
};

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename ExposureT>
PTR(afwDet::Astrometry) SillyAstrometer<ExposureT>::measureOne(ExposurePatch<ExposureT> const& patch,
                                                               afwDet::Source const& source) const
{
    CONST_PTR(afwDet::Peak) peak = patch.getPeak();
    double const NaN = std::numeric_limits<double>::quiet_NaN();
    
    return boost::make_shared<afwDet::Astrometry>(peak->getFx() + 1.0, NaN,
                                                  peak->getFy() + 1.0, NaN);
}

DECLARE_ALGORITHM(SillyAstrometer, afwDet::Astrometry);

}}}}
