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
#include "all.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/**
 * @brief A class that knows how to calculate centroids as a simple unweighted first moment
 * of the 3x3 region around a pixel
 */
template<typename ExposureT>
class GaussianAstrometer : public Algorithm<afwDet::Astrometry, ExposureT>
{
public:
    typedef Algorithm<afwDet::Astrometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<GaussianAstrometer> Ptr;
    typedef boost::shared_ptr<GaussianAstrometer const> ConstPtr;

    /// Ctor
    GaussianAstrometer() : AlgorithmT() {}

    virtual std::string getName() const { return "GAUSSIAN"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<GaussianAstrometer<ExposureT> >();
    }

    virtual PTR(afwDet::Astrometry) measureNull(void) const {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<afwDet::Astrometry>(NaN, NaN, NaN, NaN);
    }

    virtual void configure(lsst::pex::policy::Policy const& policy) {}

    virtual PTR(afwDet::Astrometry) measureOne(ExposurePatch<ExposureT> const&, afwDet::Source const&) const;
};


/**
 * @brief Given an image and a pixel position, calculate a position using a Gaussian fit
 */
template<typename ExposureT>
PTR(afwDet::Astrometry) GaussianAstrometer<ExposureT>::measureOne(ExposurePatch<ExposureT> const& patch,
                                                                  afwDet::Source const& source) const
{
    CONST_PTR(ExposureT) exposure = patch.getExposure();
    CONST_PTR(afwDet::Peak) peak = patch.getPeak();
    typedef typename ExposureT::MaskedImageT::Image ImageT;
    ImageT const& image = *exposure->getMaskedImage().getImage();

    int x = static_cast<int>(peak->getIx() + 0.5);
    int y = static_cast<int>(peak->getIy() + 0.5);

    x -= image.getX0();                 // work in image Pixel coordinates
    y -= image.getY0();

    FittedModel fit = twodg(image, x, y); // here's the fitter

    if (fit.params[FittedModel::PEAK] <= 0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException,
                          (boost::format("Object at (%d, %d) has a peak of %g") %
                           x % y % fit.params[FittedModel::PEAK]).str());
    }

    double const NaN = std::numeric_limits<double>::quiet_NaN();
    return boost::make_shared<afwDet::Astrometry>(
        lsst::afw::image::indexToPosition(image.getX0()) + fit.params[FittedModel::X0], NaN,
        lsst::afw::image::indexToPosition(image.getY0()) + fit.params[FittedModel::Y0], NaN);
}

// Declare the existence of a "GAUSSIAN" algorithm to MeasureAstrometry
LSST_DECLARE_ALGORITHM(GaussianAstrometer, afwDet::Astrometry);

}}}}
