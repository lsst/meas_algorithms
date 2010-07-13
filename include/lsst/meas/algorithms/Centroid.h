// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_CENTROID_H)
#define LSST_MEAS_ALGORITHMS_CENTROID_H 1
/**
 * @file
 */
#include <cmath>
#include <utility>
#include "boost/shared_ptr.hpp"

namespace lsst {
    namespace afw {
        namespace detection {
            class Psf;
        }
    }
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

}}}
#endif
