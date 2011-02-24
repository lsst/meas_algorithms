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
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
class NaiveAstrometry : public afwDetection::Astrometry
{
public:
    typedef boost::shared_ptr<NaiveAstrometry> Ptr;
    typedef boost::shared_ptr<NaiveAstrometry const> ConstPtr;

    /// Ctor
    NaiveAstrometry(double x, double xErr, double y, double yErr) : afwDetection::Astrometry(x, xErr, y, yErr) {}

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Astrometry::defineSchema(schema);
    }

    template<typename ExposureT>
    static Astrometry::Ptr doMeasure(typename ExposureT::ConstPtr im, afwDetection::Peak const*);

    static bool doConfigure(lsst::pex::policy::Policy const& policy)
    {
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        
        return true;
    }
private:
    static double _background;
    NaiveAstrometry(void) : afwDetection::Astrometry() { }
    LSST_SERIALIZE_PARENT(afwDetection::Astrometry)
};

LSST_REGISTER_SERIALIZER(NaiveAstrometry)

double NaiveAstrometry::_background = 0.0; // the frame's background level
    
/**
 * @brief Given an image and a pixel position, return a Centroid using a naive 3x3 weighted moment
 */
template<typename ExposureT>
afwDetection::Astrometry::Ptr NaiveAstrometry::doMeasure(typename ExposureT::ConstPtr exposure,
                                                         afwDetection::Peak const* peak)
{
    double const posErr = std::numeric_limits<double>::quiet_NaN();
    if (!peak) {
        double const pos = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<NaiveAstrometry>(pos, posErr, pos, posErr);
    }

    typedef typename ExposureT::MaskedImageT::Image ImageT;
    ImageT const& image = *exposure->getMaskedImage().getImage();

    int x = peak->getIx();
    int y = peak->getIy();

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    if (x < 1 || x >= image.getWidth() - 1 || y < 1 || y >= image.getHeight() - 1) {
         throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                           (boost::format("Object at (%d, %d) is too close to the edge") % x % y).str());
    }

    typename ImageT::xy_locator im = image.xy_at(x, y);

    double const sum =
        (im(-1,  1) + im( 0,  1) + im( 1,  1) +
         im(-1,  0) + im( 0,  0) + im( 1,  0) +
         im(-1, -1) + im( 0, -1) + im( 1, -1)) - 9*_background;

    if (sum == 0.0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has no counts") %
                           peak->getIx() % peak->getIy()).str());
    }

    double const sum_x =
        -im(-1,  1) + im( 1,  1) +
        -im(-1,  0) + im( 1,  0) +
        -im(-1, -1) + im( 1, -1);
    double const sum_y =
        (im(-1,  1) + im( 0,  1) + im( 1,  1)) -
        (im(-1, -1) + im( 0, -1) + im( 1, -1));

    return boost::make_shared<NaiveAstrometry>(
        lsst::afw::image::indexToPosition(x + image.getX0()) + sum_x/sum, posErr,
        lsst::afw::image::indexToPosition(y + image.getY0()) + sum_y/sum, posErr);
}

/*
 * Declare the existence of a "NAIVE" algorithm to MeasureAstrometry
 *
 * \cond
 */
#define INSTANTIATE(TYPE) \
    MeasureAstrometry<afwImage::Exposure<TYPE> >::declare("NAIVE", \
        &NaiveAstrometry::doMeasure<afwImage::Exposure<TYPE> >, \
        &NaiveAstrometry::doConfigure \
        )

volatile bool isInstance[] = {
    INSTANTIATE(int),
    INSTANTIATE(float),
    INSTANTIATE(double)
};

// \endcond

}}}}
