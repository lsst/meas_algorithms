#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/algorithms/Shape.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;

namespace lsst { namespace meas { namespace algorithms {

/************************************************************************************************************/
/**
 * @brief the Shape class, to hold object moments and their covariances
 */
/**
 * Return e1 == (m_{xx} - m_{yy})/(m_{xx} + m_{yy})
 */
double Shape::getE1() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    return (_mxx - _myy)/T;
}

/**
 * Return e2 == 2*m_{xy}/(m_{xx} + m_{yy})
 */
double Shape::getE2() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    return 2*_mxy/T;
}

/**
 * Return the object's RMS size, sqrt(0.5*(m_xx + m_yy))
 */
double Shape::getRms() const {
    double const T  = (_mxx + _myy);      // Trace
    return sqrt(0.5*T);
}

/**
 * Return e1's standard deviation
 */
double Shape::getE1Err() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    double const T4 = 4/(T*T*T*T);      // for our convenience
    
    double const e1e1Err = T4*(_covar(1, 1)*_myy*_myy + _covar(2, 2)*_mxx*_mxx - 2.0*_covar(1, 2)*_mxx*_myy);
    return sqrt(e1e1Err);
}

/**
 * Return sign(covar_{e1,e2})*sqrt(|covar_{e1,e2}|)
 */
double Shape::getE1E2Err() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    double const T4 = 4/(T*T*T*T);      // for our convenience
    
    double const var_e1e2 = T4*(-1.0*_myy*_mxy*_covar(1, 1) + _mxx*_mxy*_covar(2, 2) +
                                (_mxx - _myy)*_mxy*_covar(1, 2) + T*(_myy*_covar(1, 3) - _mxx*_covar(2, 3)));

    return (var_e1e2 > 0) ? sqrt(var_e1e2) : -sqrt(-var_e1e2);
}

/**
 * Return e2's standard deviation
 */
double Shape::getE2Err() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    double const T4 = 4/(T*T*T*T);      // for our convenience
    
    double const var_e2e2 = T4*(_mxy*_mxy*(_covar(1, 1) + _covar(2, 2) + 2.0*_covar(1, 2)) -
                                2.0*T*_mxy*(_covar(1, 3) + _covar(2, 3)) + T*T*_covar(3, 3));

    return sqrt(var_e2e2);
}

/**
 * Return the object's RMS size's standard deviation
 */
double Shape::getRmsErr() const {
    double const T  = (_mxx + _myy);      // Trace
    double const var_T = _covar(1, 1) + _covar(2, 2) + 2*_covar(1, 2); // T's variance
    double const ms = 0.5*T;            // rms is sqrt(ms))
    double const var_ms = 0.25*var_T;   // ms's variance

    assert(ms > 0);
    return 0.5*sqrt(var_ms/ms);          // error in sqrt(ms)
}
            
/************************************************************************************************************/
/**
 * @brief Call the concrete shape algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Shape MeasureShape<ImageT>::apply(ImageT const& image, ///< The image containing the object
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
    pexLogging::TTrace<8>("meas.algorithms.shape", "Measuring shape of object at (%.3f, %.3f)", xcen, ycen);

    return doApply(image, xcen, ycen, psf, background);
}

/************************************************************************************************************/
/*
 * Register a factory object by name;  if the factory's NULL, return the named factory
 */
template<typename ImageT>
MeasureShapeFactoryBase<ImageT>& MeasureShape<ImageT>::_registry(std::string name,
                                                                 MeasureShapeFactoryBase<ImageT>* factory) {
    static std::map<std::string const, MeasureShapeFactoryBase<ImageT> *> _registry;

    typename std::map<std::string const,
                      MeasureShapeFactoryBase<ImageT> *>::iterator el = _registry.find(name);

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
 * Return a MeasureShape of the requested variety
 *
 * @throws std::runtime_error if name can't be found
 */
template<typename ImageT>
MeasureShape<ImageT>* createMeasureShape(std::string const& name ///< desired variety
                                        ) {
    return MeasureShape<ImageT>::lookup(name).create();
}

/************************************************************************************************************/
//
// Explicit instantiations
// \cond
#define MAKE_SHAPEFINDERS(IMAGE_T) \
            template Shape MeasureShape<IMAGE_T>::apply(IMAGE_T const&, double, double, \
                                                        PSF const*, double) const; \
            template MeasureShape<IMAGE_T>* createMeasureShape<IMAGE_T>(std::string const&);
                
MAKE_SHAPEFINDERS(lsst::afw::image::MaskedImage<float>)

// \endcond
                
}}}
