#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw.h"

#include "lsst/meas/algorithms/Shape.h"
#include "lsst/meas/algorithms/ShapeImpl.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;

/*
 * Include concrete implementations
 */
#include "SdssShape.h"

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
    
    double const e1e1Err = T4*(_covar.matrix[1][1]*_myy*_myy + _covar.matrix[2][2]*_mxx*_mxx -
                               2.0*_covar.matrix[1][2]*_mxx*_myy);
    return sqrt(e1e1Err);
}

/**
 * Return sign(covar_{e1,e2})*sqrt(|covar_{e1,e2}|)
 */
double Shape::getE1E2Err() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    double const T4 = 4/(T*T*T*T);      // for our convenience
    
    double const var_e1e2 = T4*(-1.0*_myy*_mxy*_covar.matrix[1][1] + _mxx*_mxy*_covar.matrix[2][2] +
                                (_mxx - _myy)*_mxy*_covar.matrix[1][2] +
                                T*(_myy*_covar.matrix[1][3] - _mxx*_covar.matrix[2][3]));

    return (var_e1e2 > 0) ? sqrt(var_e1e2) : -sqrt(-var_e1e2);
}

/**
 * Return e2's standard deviation
 */
double Shape::getE2Err() const {
    double const T  = (_mxx + _myy);      // Trace
    assert(T > 0);
    double const T4 = 4/(T*T*T*T);      // for our convenience
    
    double const var_e2e2 = T4*(_mxy*_mxy*(_covar.matrix[1][1] + _covar.matrix[2][2] + 2.0*_covar.matrix[1][2]) -
                                2.0*T*_mxy*(_covar.matrix[1][3] + _covar.matrix[2][3]) +
                                T*T*_covar.matrix[3][3]);

    return sqrt(var_e2e2);
}

/**
 * Return the object's RMS size's standard deviation
 */
double Shape::getRmsErr() const {
    double const T  = (_mxx + _myy);      // Trace
    double const var_T = _covar.matrix[1][1] + _covar.matrix[2][2] + 2*_covar.matrix[1][2]; // T's variance
    double const ms = 0.5*T;            // rms is sqrt(ms))
    double const var_ms = 0.25*var_T;   // ms's variance

    assert(ms > 0);
    return 0.5*sqrt(var_ms/ms);          // error in sqrt(ms)
}
            
/************************************************************************************************************/
/**
 * @brief The mapping between type names (e.g. "SDSS") and an enum (lsst::meas::algorithms::SDSS)
 */
template<typename ImageT>
std::map<std::string, shapeType>* ShapeFinder<ImageT>::_shapeTypes = NULL;

/**
 * @brief Register a (name, enum) pair.
 *
 * This routine should only be called by createShapeFinder
 */
template<typename ImageT>
void ShapeFinder<ImageT>::registerType(std::string const&name, shapeType type) {
    if (_shapeTypes == NULL) {
        _shapeTypes = new(std::map<std::string, shapeType>);
    }

    (*_shapeTypes)[name] = type;
}

/**
 * @brief Return the typename for this ShapeFinder
 *
 * Names are registered using registerType
 */
template<typename ImageT>
shapeType ShapeFinder<ImageT>::lookupType(std::string const& name ///< Name of this type of shapeFinder
                                           ) {
    assert (_shapeTypes != NULL);
    
    std::map<std::string, shapeType>::const_iterator i = _shapeTypes->find(name);
    if (i == _shapeTypes->end()) {
        throw LSST_EXCEPT(pexExceptions::NotFoundException,
                          (boost::format("Unknown shape algorithm: %s") % name).str());
    }

    return i->second;
}

/**
 * @brief Call the concrete shape algorithm
 *
 * N.b. One purpose of this routine is to provide a place to specify default values for arguments
 */
template<typename ImageT>
Shape ShapeFinder<ImageT>::apply(ImageT const& image,
                                   int x,
                                   int y,
                                   PSF const* psf,
                                   double background
                                  ) const {
    if (x < 1 || x > image.getWidth() - 2 || y < 1 || y > image.getHeight() - 2) {
        throw LSST_EXCEPT(pexExceptions::RangeErrorException,
                          (boost::format("Object at (%d, %d) is too close to the edge of the frame") % x % y).str());
    }
    pexLogging::TTrace<8>("meas.algorithms.shape", "Measuring shape of object at (%d, %d)", x, y);

    return doApply(image, x, y, psf, background);
}

/**
 * @brief A factory function to return a ShapeFinder of the specified type, given as a string.
 *
 * The ShapeFinder has a method (apply) that can be used to return a Shape
 */
template<typename ImageT>
ShapeFinder<ImageT>* createShapeFinder(std::string const& type) {
    switch (ShapeFinder<ImageT>::lookupType(type)) {
      case SDSS:
        return SdssShapeFinder<ImageT>::getInstance();
      default:
        throw LSST_EXCEPT(pexExceptions::NotFoundException, 
                          (boost::format("ShapeFinder of type %d is not implemented") % type).str());
    }
    // NOTREACHED
}

//
// Explicit instantiations
// \cond
#define MAKE_SHAPEFINDERS(IMAGE_T) \
                template Shape ShapeFinder<IMAGE_T>::apply(IMAGE_T const&, int, int, PSF const*, double) const; \
                template ShapeFinder<IMAGE_T>* createShapeFinder<IMAGE_T>(std::string const&); \
                template void ShapeFinder<IMAGE_T>::registerType(std::string const&name, shapeType type); \
                template shapeType ShapeFinder<IMAGE_T>::lookupType(std::string const&name);
                
MAKE_SHAPEFINDERS(lsst::afw::image::Image<float>)

// \endcond
                
}}}
