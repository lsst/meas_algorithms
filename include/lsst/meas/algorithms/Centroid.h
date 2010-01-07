// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_CENTROID_H)
#define LSST_MEAS_ALGORITHMS_CENTROID_H 1
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
class Centroid {
public:
    typedef boost::shared_ptr<Centroid> Ptr;
    typedef boost::shared_ptr<Centroid const> ConstPtr;
    typedef std::pair<double, double> xyAndError;

    Centroid(double x = NAN, double y = NAN) : _x(x), _xErr(NAN), _y(y), _yErr(NAN), _xyCovar(NAN) {}
    Centroid(xyAndError x, xyAndError y, double xyCovar = NAN) :
        _x(x.first), _xErr(x.second), _y(y.first), _yErr(y.second), _xyCovar(xyCovar) {}
    virtual ~Centroid() {}
    
    double getX() const { return _x; }
    double getXErr() const { return _xErr; }
    xyAndError getX(int dummy) const { return xyAndError(_x, _xErr); }
    void setX(double const x) { _x = x; }
    void setX(xyAndError const& vt) { _x = vt.first; _xErr = vt.second; }
    void setXErr(double const xErr) { _xErr = xErr; }

    double getY() const { return _y; }
    double getYErr() const { return _yErr; }
    xyAndError getY(int dummy) const { return xyAndError(_y, _yErr); }
    void setY(double const y) { _y = y; }
    void setY(xyAndError const& vt) { _y = vt.first; _yErr = vt.second; }
    void setYErr(double const yErr) { _yErr = yErr; }

    double getCovar() const { return _xyCovar; }
    void setCovar(double const xyCovar) { _xyCovar = xyCovar; }
private:
    double _x, _xErr;                   // column position + error (sqrt(variance))
    double _y, _yErr;                   // row position + error
    double _xyCovar;                    // covariance of x and y
};

/************************************************************************************************************/
            
template<typename ImageT> class MeasureCentroid;
/*
 * Must be here, as it's declared a friend by MeasureCentroid
 */
template<typename ImageT> MeasureCentroid<ImageT> *createMeasureCentroid(std::string const& type);
            
/**
 * A polymorphic base class for MeasureCentroid factories
 */
template<typename ImageT>
class MeasureCentroidFactoryBase : public lsst::daf::base::Citizen {
public:
    MeasureCentroidFactoryBase() : lsst::daf::base::Citizen(typeid(this)) {}
    virtual ~MeasureCentroidFactoryBase() {}
    virtual MeasureCentroid<ImageT> *create() = 0;
};

/**
 * Create a particular sort of MeasureCentroid
 */
template<typename MeasureCentroidT>
class MeasureCentroidFactory : public MeasureCentroidFactoryBase<typename MeasureCentroidT::ImageT> {
public:
    /**
     * Return a new MeasureCentroidT
     */
    MeasureCentroidT *create() {
        return new MeasureCentroidT();
    }
};

/**
 * @brief A pure virtual base class to calculate a centroid
 *
 * Different implementations will use different algorithms
 */
template<typename T>
class MeasureCentroid : public boost::noncopyable {
public:
    typedef T ImageT;
    typedef boost::shared_ptr<MeasureCentroid> Ptr;
    typedef boost::shared_ptr<MeasureCentroid const> ConstPtr;

    MeasureCentroid() {
        static bool _registered = false;

        if (!_registered) {
#if 0                                   // We don't actually declare the (pure virtual) base class
            MeasureCentroid::declare("base", new MeasureCentroidFactory<ImageT>());
#endif
            _registered = true;
        }        
    }
    virtual ~MeasureCentroid() {}

    Centroid apply(ImageT const& image, int x, int y,
                   lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                   double background = 0.0) const;
protected:
#if !defined(SWIG)
    friend MeasureCentroid *createMeasureCentroid<>(std::string const& name);
#endif

    /**
     * Declare a MeasureCentroidFactory for a variety "name"
     *
     * @throws std::runtime_error if name is already declared
     */
    static void declare(std::string name,                        ///< name of variety
                        MeasureCentroidFactoryBase<ImageT>* factory ///< Factory to make this sort of widget
                       ) {
        (void)_registry(name, factory);
    }

    /**
     * Return the named MeasureCentroidFactory
     *
     * @throws std::runtime_error if name can't be found
     */
    static MeasureCentroidFactoryBase<ImageT>& lookup(std::string name ///< desired variety
                                                     ) {
        return _registry(name, NULL);
    }
private:
    static MeasureCentroidFactoryBase<ImageT>& _registry(std::string name,
                                                         MeasureCentroidFactoryBase<ImageT>* factory);
    virtual Centroid doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const = 0;
};

}}}
#endif
