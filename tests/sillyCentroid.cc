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
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that knows how to calculate centroids by guessing the wrong answer
 */
class SillyAstrometry : public afwDetection::Astrometry
{
public:
    typedef boost::shared_ptr<SillyAstrometry> Ptr;
    typedef boost::shared_ptr<SillyAstrometry const> ConstPtr;

    /// Ctor
    SillyAstrometry(double x, double xErr, double y, double yErr)
    {
        init();                         // This allocates space for fields added by defineSchema
        set<X>(x);                      // ... if you don't, these set calls will fail an assertion
        set<X_ERR>(xErr);               // the type of the value must match the schema
        set<Y>(y);
        set<Y_ERR>(yErr);
    }

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Astrometry::defineSchema(schema);
    }

    template<typename ExposureT>
    static Astrometry::Ptr doMeasure(CONST_PTR(ExposureT),
                                     CONST_PTR(afwDetection::Peak),
                                     CONST_PTR(afwDetection::Source)
                                    );
};

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename ExposureT>
afwDetection::Astrometry::Ptr SillyAstrometry::doMeasure(CONST_PTR(ExposureT) image,
                                                         CONST_PTR(afwDetection::Peak) peak,
                                                         CONST_PTR(afwDetection::Source)
                                                        )
{
    double const posErr = std::numeric_limits<double>::quiet_NaN();
    if (!peak) {
        double const pos = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<SillyAstrometry>(pos, posErr, pos, posErr);
    }
    
    return boost::make_shared<SillyAstrometry>(peak->getFx() + 1.0, posErr,
                                               peak->getFy() + 1.0, posErr);
}

/*
 * Declare the existence of a "SILLY" algorithm to MeasureAstrometry
 *
 * @cond
 */
#define INSTANTIATE(TYPE) \
    MeasureAstrometry<afwImage::Exposure<TYPE> >::declare("SILLY", \
        &SillyAstrometry::doMeasure<afwImage::Exposure<TYPE> > \
        )

volatile bool isInstance[] = {
    INSTANTIATE(int),
    INSTANTIATE(float),
};

// \endcond

}}}}
