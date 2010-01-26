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
#include "lsst/meas/algorithms/detail/MeasureFactory.h"

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
    xyAndError getX(int) const { return xyAndError(_x, _xErr); }
    void setX(double const x) { _x = x; }
    void setX(xyAndError const& vt) { _x = vt.first; _xErr = vt.second; }
    void setXErr(double const xErr) { _xErr = xErr; }

    double getY() const { return _y; }
    double getYErr() const { return _yErr; }
    xyAndError getY(int) const { return xyAndError(_y, _yErr); }
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
/**
 * @brief A pure virtual base class to calculate a centroid
 *
 * Different implementations will use different algorithms
 */
template<typename T>
class MeasureCentroid : public MeasureProperty<MeasureCentroid<T>, T> {
public:
    typedef T ImageT;

    typedef boost::shared_ptr<MeasureCentroid> Ptr;
    typedef boost::shared_ptr<MeasureCentroid const> ConstPtr;

    MeasureCentroid(typename ImageT::ConstPtr image=typename ImageT::ConstPtr())
        : MeasureProperty<MeasureCentroid<T>, T>(image) {}
    virtual ~MeasureCentroid() {}

    Centroid apply(ImageT const& image, int x, int y,
                   lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                   double background = 0.0) const;

    Centroid apply(int x, int y,
                   lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                   double background = 0.0) const {
        if (!this->getImage()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "You must provide an image to measure");
        }
        return apply(*this->getImage(), x, y, psf, background);
    }
private:
    virtual Centroid doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const = 0;
};

/**
 * Provide a convenient wrapper round createMeasureProperty
 */
template<typename ImageT>
MeasureCentroid<ImageT> *
createMeasureCentroid(
        std::string const& type, ///< Algorithm type (e.g. "NAIVE")
        boost::shared_ptr<ImageT const> image=boost::shared_ptr<ImageT const>() ///< The image to process
                     )
{
    lsst::meas::algorithms::MeasureCentroid<ImageT> const* ptr = NULL;
    return createMeasureProperty(type, image, ptr);
}

}}}
#endif
