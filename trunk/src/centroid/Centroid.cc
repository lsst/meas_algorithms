#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/algorithms/Centroid.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

/**
 * @brief Call the concrete centroiding algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Centroid MeasureCentroid<ImageT>::apply(ImageT const& image,
                                   int x,
                                   int y,
                                   PSF const* psf,
                                   double background
                                  ) const {
    if (x - image.getX0() < 1 || x - image.getX0() > image.getWidth() - 2 ||
        y - image.getY0() < 1 || y - image.getY0() > image.getHeight() - 2) {
            throw LSST_EXCEPT(pexExceptions::RangeErrorException,
                              (boost::format("Object at (%d, %d) is too close "
                                             "to the edge of the frame") % x % y).str());
    }
    pexLogging::TTrace<8>("meas.algorithms.centroid", "Centroiding object at (%d, %d)", x, y);

    return doApply(image, x, y, psf, background);
}

/************************************************************************************************************/
/*
 * Register a factory object by name;  if the factory's NULL, return the named factory
 */
template<typename ImageT>
MeasureCentroidFactoryBase<ImageT>&
        MeasureCentroid<ImageT>::_registry(std::string name,
                                           MeasureCentroidFactoryBase<ImageT>* factory) {
    static std::map<std::string const, MeasureCentroidFactoryBase<ImageT> *> _registry;

    typename std::map<std::string const,
                      MeasureCentroidFactoryBase<ImageT> *>::iterator el = _registry.find(name);

    if (el == _registry.end()) {        // failed to find name
        if (factory) {
            _registry[name] = factory;
        } else {
            throw LSST_EXCEPT(pexExceptions::NotFoundException, 
                              "MeasureCentroid of type \"" + name + "\" is not implemented");
        }
    } else {
        if (!factory) {
            factory = (*el).second;
        } else if(factory == (*el).second) {
            ;                           // OK
        } else {
            throw LSST_EXCEPT(pexExceptions::InvalidParameterException, 
                              "MeasureCentroid of type \"" + name + "\" is already declared");
        }
    }

    return *factory;
}

/************************************************************************************************************/
/**
 * Return a MeasureCentroid of the requested variety
 *
 * @throws std::runtime_error if name can't be found
 */
template<typename ImageT>
MeasureCentroid<ImageT>* createMeasureCentroid(std::string const& name ///< desired variety
                                              ) {
    return MeasureCentroid<ImageT>::lookup(name).create();
}

/************************************************************************************************************/
//
// Explicit instantiations
// \cond
#define MAKE_CENTROIDERS(IMAGE_T) \
                template Centroid MeasureCentroid<IMAGE_T>::apply(IMAGE_T const&, int, int, \
                                                                  PSF const*, double) const; \
                template MeasureCentroid<IMAGE_T>* createMeasureCentroid<IMAGE_T>(std::string const&);
                
MAKE_CENTROIDERS(lsst::afw::image::Image<float>)

// \endcond
                
}}}
