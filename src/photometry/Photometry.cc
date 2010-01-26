// -*- LSST-C++ -*-
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
Photometry MeasurePhotometry<ImageT>::apply(ImageT const& image,
                                            double xcen,
                                            double ycen,
                                            PSF const* psf,
                                            double background
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
/*
 * Register a factory object by name;  if the factory's NULL, return the named factory
 */
template<typename ImageT>
MeasurePhotometryFactoryBase<ImageT>&
MeasurePhotometry<ImageT>::_registry(std::string name,
                                     MeasurePhotometryFactoryBase<ImageT>* factory) {
    static std::map<std::string const, MeasurePhotometryFactoryBase<ImageT> *> _registry;
    
    typename std::map<std::string const,
        MeasurePhotometryFactoryBase<ImageT> *>::iterator el = _registry.find(name);
    
    if (el == _registry.end()) {        // failed to find name
        if (factory) {
            _registry[name] = factory;
        } else {
            throw LSST_EXCEPT(pexExceptions::NotFoundException, 
                              "MeasurePhotometry of type \"" + name + "\" is not implemented");
        }
    } else {
        if (!factory) {
            factory = (*el).second;
        } else if(factory == (*el).second) {
            ;                           // OK
        } else {
            throw LSST_EXCEPT(pexExceptions::InvalidParameterException, 
                              "MeasurePhotometry of type \"" + name + "\" is already declared");
        }
    }
    
    return *factory;
}
    
/************************************************************************************************************/
/**
 * Return a MeasurePhotometry of the requested variety
 *
 * @throws std::runtime_error if name can't be found
 */
template<typename ImageT>
MeasurePhotometry<ImageT>* createMeasurePhotometry(std::string const& name, ///< desired variety
                                                   float const radius       ///< aperture radius
                                              ) {
    return MeasurePhotometry<ImageT>::lookup(name).create(radius);
}

/************************************************************************************************************/
//
// Explicit instantiations
// \cond
#define MAKE_PHOTOMETRYS(IMAGE_T) \
    template class MeasurePhotometry<IMAGE_T>; \
    template MeasurePhotometry<IMAGE_T>* createMeasurePhotometry<IMAGE_T>(std::string const&, float const);

MAKE_PHOTOMETRYS(afwImage::MaskedImage<float>)
#if 0
MAKE_PHOTOMETRYS(afwImage::MaskedImage<double>)
#endif

// \endcond
                
}}}
