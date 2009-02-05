#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/CentroidImpl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

/*
 * Include concrete implementations
 */
#include "lsst/meas/algorithms/NaiveCentroid.h"
#include "lsst/meas/algorithms/SdssCentroid.h"

namespace lsst { namespace meas { namespace algorithms {

/************************************************************************************************************/
/**
 * @brief The mapping between type names (e.g. "SDSS") and an enum (lsst::meas::algorithms::SDSS)
 */
template<typename ImageT>
std::map<std::string, centroidType>* Centroider<ImageT>::_centroidTypes = NULL;

/**
 * @brief The string name of the Centroider (e.g. "SDSS")
 *
 * This name may be used to "format" a Centroider, and recreate it on a remote machine,
 * or to specify the desired algorithm in a Policy file
 */
template<typename ImageT>
std::string Centroider<ImageT>::_typeName = "";

/**
 * @brief Register a (name, enum) pair.
 *
 * This routine should only be called by make_Centroider
 */
template<typename ImageT>
void Centroider<ImageT>::registerType(std::string const&name, centroidType type) {
    if (_centroidTypes == NULL) {
        _centroidTypes = new(std::map<std::string, centroidType>);
    }

    _typeName = name;
    (*_centroidTypes)[name] = type;
}

/**
 * @brief Return the typename for this Centroider
 *
 * Names are registered using registerType
 */
template<typename ImageT>
centroidType Centroider<ImageT>::lookupType(std::string const& name ///< Name of this type of centroider
                                           ) {
    assert (_centroidTypes != NULL);
    
    std::map<std::string, centroidType>::const_iterator i = _centroidTypes->find(name);
    if (i == _centroidTypes->end()) {
        throw LSST_EXCEPT(pexExceptions::NotFoundException,
                          (boost::format("Unknown centroiding algorithm: %s") % name).str());
    }

    return i->second;
}

/**
 * @brief Call the concrete centroiding algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Centroid Centroider<ImageT>::apply(ImageT const& image,
                                   int x,
                                   int y,
                                   PSF const* psf,
                                   double background
                                  ) const {
    if (x < 1 || x > image.getWidth() - 2 || y < 1 || y > image.getHeight() - 2) {
        throw LSST_EXCEPT(pexExceptions::RangeErrorException,
                          (boost::format("Object at (%d, %d) is too close to the edge of the frame") % x % y).str());
    }
    pexLogging::TTrace<8>("meas.algorithms.centroid", "Centroiding object at (%d, %d)", x, y);

    return doApply(image, x, y, psf, background);
}

/**
 * @brief A factory function to return a Centroider of the specified type, given as a string.
 *
 * The Centroider has a method (apply) that can be used to return a Centroid
 */
template<typename ImageT>
Centroider<ImageT>* make_Centroider(std::string const& type) {
    switch (Centroider<ImageT>::lookupType(type)) {
      case NAIVE:
        return NaiveCentroider<ImageT>::getInstance();
      case SDSS:
        return SdssCentroider<ImageT>::getInstance();
      default:
        throw LSST_EXCEPT(pexExceptions::NotFoundException, 
                          (boost::format("Centroider of type %d is not implemented") % type).str());
    }
    // NOTREACHED
}

//
// Explicit instantiations
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
                template Centroid Centroider<IMAGE_T>::apply(IMAGE_T const&, int, int, PSF const*, double) const; \
                template Centroider<IMAGE_T>* make_Centroider<IMAGE_T>(std::string const&); \
                template void Centroider<IMAGE_T>::registerType(std::string const&name, centroidType type); \
                template centroidType Centroider<IMAGE_T>::lookupType(std::string const&name);
                
MAKE_CENTROIDERS(lsst::afw::image::Image<float>)

// \endcond
                
}}}
