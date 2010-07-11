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
 
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Centroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

namespace {
/**
 * @brief A class that knows how to calculate centroids
 */
template<typename ImageT>
class SillyMeasureCentroid : public MeasureCentroid<ImageT> {
public:
    typedef MeasureCentroid<ImageT> MeasurePropertyBase;

    SillyMeasureCentroid(typename ImageT::ConstPtr image) : MeasureCentroid<ImageT>(image) {}
private:
    Centroid doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const;
};

/**
 * @brief Given an image and a pixel position, return a Centroid offset by (1, 1) from initial position
 */
template<typename ImageT>
Centroid SillyMeasureCentroid<ImageT>::doApply(ImageT const&, ///< The Image wherein dwells the object
                                               int x,         ///< object's column position
                                               int y,         ///< object's row position
                                               PSF const*,    ///< image's PSF
                                               double         ///< image's background level
                                              ) const
{
    return Centroid(lsst::afw::image::indexToPosition(x) + 1.0,
                    lsst::afw::image::indexToPosition(y) + 1.0);
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasureCentroid
//
// \cond
#define MAKE_CENTROIDER(IMAGE_T) \
    registerMe<SillyMeasureCentroid, lsst::afw::image::Image<IMAGE_T> >("SILLY")
                
volatile bool isInstance[] = {
    MAKE_CENTROIDER(int),
    MAKE_CENTROIDER(float),
};

// \endcond

}}}}
