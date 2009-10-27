#if !defined(LSST_MEAS_ALGORITHMS_PHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_PHOTOMETRY_H 1
/**
 * @file Photometry.h
 * @brief Compute aperture and PSF photometry.
 * @ingroup meas/algorithms
 * @author Steve Bickerton (adapted from RHL's Shape class)
 */
#include <cmath>
#include <utility>
#include <string>
#include "boost/shared_ptr.hpp"
#include "boost/noncopyable.hpp"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/PSF.h"

namespace lsst { namespace meas { namespace algorithms {
/**
 * @class Photometry
 */
class Photometry {
public:
    typedef boost::shared_ptr<Photometry> Ptr;
    typedef boost::shared_ptr<const Photometry> ConstPtr;
    
    Photometry(double apflux=NAN, double psfflux=NAN) :
        _apflux(apflux),
        _psfflux(psfflux) {
    }
    ~Photometry() {}

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

/**
 * Types of supported photometry algorithms
 */
typedef int photometryType;

/**
 * @brief A pure virtual base class to calculate aperture and PSF fluxes
 *
 * Different implementations will use different algorithms
 */
template<typename ImageT>
class measurePhotometry : public boost::noncopyable {
public:
    typedef boost::shared_ptr<measurePhotometry> Ptr;
    typedef boost::shared_ptr<measurePhotometry const> ConstPtr;

    measurePhotometry(float const radius) : _radius(radius) {}
    virtual ~measurePhotometry() {}

    Photometry apply(ImageT const& image, double xcen, double ycen,
                     lsst::meas::algorithms::PSF const* psf=NULL, // fully qualified to make swig happy
                     double background=0.0) const;
    
    static photometryType lookupType(std::string const& name);

    void setRadius(float const radius) const {
        _radius = radius;
    }
        
protected:
    static void registerType(std::string const& name, photometryType type);
    float mutable _radius;
private:
    virtual Photometry doApply(ImageT const& image, double xcen, double ycen,
                               PSF const* psf, double background) const = 0;
    
    static std::map<std::string, photometryType>* _photometryTypes;
};
            
template<typename ImageT>
measurePhotometry<ImageT>* createMeasurePhotometry(std::string const& type, float const radius);
            
}}}
#endif
