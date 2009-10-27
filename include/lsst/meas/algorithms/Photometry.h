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
            
template<typename ImageT> class MeasurePhotometry;
/*
 * Must be here, as it's declared a friend by MeasurePhotometry
 */
template<typename ImageT> MeasurePhotometry<ImageT> *createMeasurePhotometry(std::string const& type,
                                                                             float const radius);
            
/**
 * A polymorphic base class for MeasurePhotometry factories
 */
template<typename ImageT>
class MeasurePhotometryFactoryBase : public lsst::daf::base::Citizen {
public:
    MeasurePhotometryFactoryBase() : lsst::daf::base::Citizen(typeid(this)) {}
    virtual ~MeasurePhotometryFactoryBase() {}
    virtual MeasurePhotometry<ImageT> *create(float const radius) = 0;
};

/**
 * Create a particular sort of MeasurePhotometry
 */
template<typename MeasurePhotometryT>
class MeasurePhotometryFactory : public MeasurePhotometryFactoryBase<typename MeasurePhotometryT::ImageT> {
public:
    /**
     * Return a new MeasurePhotometryT
     */
    MeasurePhotometryT *create(float const radius) {
        return new MeasurePhotometryT(radius);
    }
};

/**
 * @brief A pure virtual base class to calculate a photometry
 *
 * Different implementations will use different algorithms
 */
template<typename T>
class MeasurePhotometry : public boost::noncopyable {
public:
    typedef T ImageT;
    typedef boost::shared_ptr<MeasurePhotometry> Ptr;
    typedef boost::shared_ptr<MeasurePhotometry const> ConstPtr;

    explicit MeasurePhotometry(float const radius) : _radius(radius) {
        static bool _registered = false;

        if (!_registered) {
#if 0                                   // We don't actually declare the (pure virtual) base class
            MeasurePhotometry::declare("base", new MeasurePhotometryFactory<ImageT>());
#endif
            _registered = true;
        }        
    }
    virtual ~MeasurePhotometry() {}

    Photometry apply(T const& image, double xcen, double ycen,
                     lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                     double background = 0.0) const;

protected:
#if !defined(SWIG)
    friend MeasurePhotometry *createMeasurePhotometry<>(std::string const& name, float const radius);
#endif

    /**
     * Declare a MeasurePhotometryFactory for a variety "name"
     *
     * @throws std::runtime_error if name is already declared
     */
    static void declare(std::string name,                        ///< name of variety
                        MeasurePhotometryFactoryBase<ImageT>* factory ///< Factory to make this sort of widget
                       ) {
        (void)_registry(name, factory);
    }

    /**
     * Return the named MeasurePhotometryFactory
     *
     * @throws std::runtime_error if name can't be found
     */
    static MeasurePhotometryFactoryBase<ImageT>& lookup(std::string name ///< desired variety
                                                       ) {
        return _registry(name, NULL);
    }
    
    inline void setRadius(float const radius) const {
        _radius = radius;
    }
    float mutable _radius;

private:
    static MeasurePhotometryFactoryBase<ImageT>& _registry(std::string name,
                                                           MeasurePhotometryFactoryBase<ImageT>* factory);
    virtual Photometry doApply(ImageT const& image,
                               double xcen, double ycen, PSF const* psf, double background) const = 0;
};

}}}
#endif
