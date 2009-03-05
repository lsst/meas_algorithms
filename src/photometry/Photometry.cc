// -*- LSST-C++ -*-
/**
 * @file Photometry.cc
 * @brief Compute aperture and PSF photometry.
 * @ingroup meas/algorithms
 * @author Steve Bickerton (adapted from RHL's Shape class)
 *
 */
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/PhotometryImpl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;

/*
 * Include concrete implementations
 */
#include "NaivePhotometry.h"

namespace lsst { namespace meas { namespace algorithms {
            

/************************************************************************************************************/
/**
 * @brief The mapping between type names (e.g. "SDSS") and an enum (lsst::meas::algorithms::SDSS)
 */
template<typename ImageT>
std::map<std::string, photometryType>* measurePhotometry<ImageT>::_photometryTypes = NULL;

/**
 * @brief Register a (name, enum) pair.
 *
 * This routine should only be called by createMeasurePhotometry
 */
template<typename ImageT>
void measurePhotometry<ImageT>::registerType(std::string const&name, photometryType type) {
    if (_photometryTypes == NULL) {
        _photometryTypes = new(std::map<std::string, photometryType>);
    }

    (*_photometryTypes)[name] = type;
}

/**
 * @brief Return the typename for this measurePhotometry
 *
 * Names are registered using registerType
 */
template<typename ImageT>
photometryType measurePhotometry<ImageT>::lookupType(std::string const& name ///< Name of this type of measurePhotometry
                                           ) {
    assert (_photometryTypes != NULL);
    
    std::map<std::string, photometryType>::const_iterator i = _photometryTypes->find(name);
    if (i == _photometryTypes->end()) {
        throw LSST_EXCEPT(pexExceptions::NotFoundException,
                          (boost::format("Unknown photometry algorithm: %s") % name).str());
    }

    return i->second;
}

/**
 * @brief Call the concrete photometry algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Photometry measurePhotometry<ImageT>::apply(ImageT const& image, ///< The image containing the object
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
                          (boost::format("Object at (%.3f, %.3f) is too close to the edge of the frame") %
                           xcen % ycen).str());
    }
    pexLogging::TTrace<8>("meas.algorithms.photometry", "Measuring photometry of object at (%.3f, %.3f)", xcen, ycen);

    return doApply(image, xcen, ycen, psf, background);
}

/**
 * @brief A factory function to return a measurePhotometry of the specified type, given as a string.
 *
 * The measurePhotometry has a method (apply) that can be used to return a Photometry
 */
template<typename ImageT>
measurePhotometry<ImageT>* createMeasurePhotometry(std::string const& type,
                                                  float const radius
                                                  ) {
    switch (measurePhotometry<ImageT>::lookupType(type)) {
      case NAIVE:
        return measureNaivePhotometry<ImageT>::getInstance(radius);
      default:
        throw LSST_EXCEPT(pexExceptions::NotFoundException, 
                          (boost::format("measurePhotometry of type %d is not implemented") % type).str());
    }
    // NOTREACHED
}

//
// Explicit instantiations
// \cond
#define MAKE_PHOTOMETRYFINDERS(IMAGE_T) \
            template Photometry measurePhotometry<IMAGE_T>::apply(IMAGE_T const&, double, double, PSF const*, double) const; \
            template measurePhotometry<IMAGE_T>* createMeasurePhotometry<IMAGE_T>(std::string const&, float const); \
                template void measurePhotometry<IMAGE_T>::registerType(std::string const&name, photometryType type); \
                template photometryType measurePhotometry<IMAGE_T>::lookupType(std::string const&name);
                
MAKE_PHOTOMETRYFINDERS(lsst::afw::image::MaskedImage<float>)

// \endcond
                
}}}
