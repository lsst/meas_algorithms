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
 
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/algorithms/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief Call the concrete photometry algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Photometry MeasurePhotometry<ImageT>::apply(ImageT const& image, ///< The image containing the object
                                            double xcen,         ///< object's column position
                                            double ycen,         ///< object's row position
                                            PSF const* psf,      ///< image's PSF
                                            double background    ///< image's background level
                                           ) const {
    int const x = afwImage::positionToIndex(xcen);
    int const y = afwImage::positionToIndex(ycen);

    if (x - image.getX0() < 1 || x - image.getX0() > image.getWidth() - 2 ||
        y - image.getY0() < 1 || y - image.getY0() > image.getHeight() - 2) {
        throw LSST_EXCEPT(pexExceptions::RangeErrorException,
                          (boost::format("Object at (%d, %d) is too close "
                                         "to the edge of the frame") % x % y).str());
    }
    pexLogging::TTrace<8>("meas.algorithms.photometry", "Photometry object at (%d, %d)", x, y);
    
    return doApply(image, x, y, psf, background);
}

/************************************************************************************************************/
/**
 * Return a MeasurePhotometry of the requested variety
 *
 * @throws std::runtime_error if name can't be found
 */
template<typename ImageT>
MeasurePhotometry<ImageT>* createMeasurePhotometry(std::string const& name ///< desired variety
                                                  ) {
    return MeasurePhotometry<ImageT>::lookup(name).create();
}

/************************************************************************************************************/
//
// Explicit instantiations
// \cond
#define MAKE_PHOTOMETRYS(IMAGE_T) \
    template class MeasurePhotometry<IMAGE_T>; \
    template MeasurePhotometry<IMAGE_T>* \
    createMeasureProperty(std::string const&, IMAGE_T::ConstPtr, MeasurePhotometry<IMAGE_T> const*)

MAKE_PHOTOMETRYS(afwImage::MaskedImage<float>);
#if 0
MAKE_PHOTOMETRYS(afwImage::MaskedImage<double>);
#endif

// \endcond
                
}}}
