// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_PHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_PHOTOMETRY_H 1
/**
 * @file
 */
#include <cmath>
#include <utility>
#include <string>
#include <typeinfo>
#include "boost/shared_ptr.hpp"
#include "boost/noncopyable.hpp"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/detail/MeasureFactory.h"

namespace lsst {
namespace meas {
namespace algorithms {
/**
 * @brief Represent a position and its error
 */
class Photometry {
public:
    typedef boost::shared_ptr<Photometry> Ptr;
    typedef boost::shared_ptr<Photometry const> ConstPtr;
    typedef std::pair<double, double> xyAndError;


    Photometry(double apflux = NAN, double psfflux = NAN) :
        _apflux(apflux),
        _psfflux(psfflux) {}
    
    virtual ~Photometry() {}

    void setApFlux(double apflux) { _apflux = apflux; }
    double getApFlux() const { return _apflux; }
    double getApFluxErr() const { return _apfluxerr; }

    void setApMag(double apmag) { _apmag = apmag; }
    double getApMag() const { return _apmag; }
    double getApMagErr() const { return _apmagerr; }

    void setPsfFlux(double psfflux) { _psfflux = psfflux; }
    double getPsfFlux() const { return _psfflux; }
    double getPsfFluxErr() const { return _psffluxerr; }

    void setPsfMag(double psfmag) { _psfmag = psfmag; }
    double getPsfMag() const { return _psfmag; }
    double getPsfMagErr() const { return _psfmagerr; }

private:

    double _apflux, _apfluxerr;
    double _apmag, _apmagerr;
    double _psfflux, _psffluxerr;
    double _psfmag, _psfmagerr;

};

/************************************************************************************************************/
/**
 * @brief A pure virtual base class to calculate a photometry
 *
 * Different implementations will use different algorithms
 */
template<typename T>
class MeasurePhotometry : public MeasureProperty<MeasurePhotometry<T>, T> {
public:
    typedef T ImageT;
    
    typedef boost::shared_ptr<MeasurePhotometry> Ptr;
    typedef boost::shared_ptr<MeasurePhotometry const> ConstPtr;

    MeasurePhotometry(typename ImageT::ConstPtr image=typename ImageT::ConstPtr())
        : MeasureProperty<MeasurePhotometry<T>, T>(image) {}
    virtual ~MeasurePhotometry() {}

    Photometry apply(ImageT const& image, double xcen, double ycen,
                     lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                     double background = 0.0) const;
    Photometry apply(int x, int y,
                     lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                     double background = 0.0) const {
        if (!this->getImage()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "You must provide an image to measure");
        }
        return apply(*this->getImage(), x, y, psf, background);
    }
    
    inline float getRadius() const {
        return _radius;
    }
    inline void setRadius(float const radius) const {
        _radius = radius;
    }

private:
    float mutable _radius;

    virtual Photometry doApply(ImageT const& image,
                               double xcen, double ycen, PSF const* psf, double background) const = 0;
};

/************************************************************************************************************/
/**
 * Provide a convenient wrapper round createMeasureProperty
 */
template<typename ImageT>
MeasurePhotometry<ImageT> *
createMeasurePhotometry(
        std::string const& type,        ///< Algorithm type (e.g. "NAIVE")
        boost::shared_ptr<ImageT const> image=boost::shared_ptr<ImageT const>() ///< The image to process
                       )
{
    MeasurePhotometry<ImageT> const* ptr = NULL;
    return createMeasureProperty(type, image, ptr);
}
    
}}}
#endif
