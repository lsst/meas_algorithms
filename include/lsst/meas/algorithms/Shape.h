#if !defined(LSST_MEAS_ALGORITHMS_SHAPE_H)
#define LSST_MEAS_ALGORITHMS_SHAPE_H 1
/**
 * @file
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
 * @brief Represent a position and its error
 */
class Shape {
    static int const ncovar = 4;        // dimension of covariance matrix, _covar

public:
#if defined(SWIG)
    #pragma SWIG nowarn=SWIGWARN_PARSE_NESTED_CLASS
#endif

    class Covar {                       // replace me with an LSST matrix class
    public:
        Covar() {
            for (int i = 0; i != ncovar; ++i) {
                for (int j = 0; j != ncovar; ++j) {
                    matrix[i][j] = NAN;
                }
            }
        }

        void invert() {
            ;                           // No this isn't correct, but I don't have a matrix inverter handy
        }

        double matrix[ncovar][ncovar];
    };

    typedef boost::shared_ptr<Shape> Ptr;
    typedef boost::shared_ptr<const Shape> ConstPtr;
    
    typedef std::pair<double, double> xyAndError;

    enum { SHIFT = 0x1, MAXITER = 0x2, UNWEIGHTED = 0x4, UNWEIGHTED_PSF = 0x8, };

    Shape(double m0=NAN, double mxx=NAN, double mxy=NAN, double myy=NAN, Centroid centroid=Centroid()) :
        _centroid(centroid),
        _m0(m0),
        _mxx(mxx), _mxy(mxy), _myy(myy),
        _covar(),
        _mxy4(NAN),
        _flags(0) {
    }
    ~Shape() {}

    Centroid& getCentroid() { return _centroid; } // this is doubtful design...
    Centroid const& getCentroid() const { return _centroid; }

    void setM0(double m0) { _m0 = m0; }
    double getM0() const { return _m0; }
    double getM0Err() const { return _covar.matrix[0][0]; }

    void setMxx(double mxx) { _mxx = mxx; }
    double getMxx() const { return _mxx; }
    double getMxxErr() const { return _covar.matrix[1][1]; }

    void setMxy(double mxy) { _mxy = mxy; }
    double getMxy() const { return _mxy; }
    double getMxyErr() const { return _covar.matrix[2][2]; }

    void setMyy(double myy) { _myy = myy; }
    double getMyy() const { return _myy; }
    double getMyyErr() const { return _covar.matrix[3][3]; }

    void setMxy4(double mxy4) { _mxy4 = mxy4; }
    double getMxy4() const { return _mxy4; }

    void clearFlags() { _flags = 0; }
    void setFlags(int flags) { _flags |= flags; }
    int getFlags() const { return _flags; }

    void setCovar(Covar covar) { _covar = covar; }
    const Covar& getCovar() const { return _covar; }
    
    double getE1() const;
    double getE1Err() const;
    double getE2() const;
    double getE2Err() const;
    double getE1E2Err() const;
    double getRms() const;
    double getRmsErr() const;
private:
    Centroid _centroid;                   // object's centroid, as measured along with moments
    double _m0;                           // 0-th moment
    double _mxx, _mxy, _myy;              // <xx> <xy> <yy>
    Covar _covar;                         // covariance matrix for (_m0, _mxx, _mxy, _myy)
    double _mxy4;                          // 4th moment used for shear calibration
    int _flags;                           // flags describing processing
};

/**
 * Types of supported shape algorithms
 */
typedef int shapeType;

/**
 * @brief A pure virtual base class to calculate a shape
 *
 * Different implementations will use different algorithms
 */
template<typename ImageT>
class measureShape : public boost::noncopyable {
public:
    typedef boost::shared_ptr<measureShape> Ptr;

    measureShape() {}
    virtual ~measureShape() {}

    Shape apply(ImageT const& image, int x, int y,
                   lsst::meas::algorithms::PSF const* psf=NULL, // fully qualified to make swig happy
                   double background=0.0) const;

    static shapeType lookupType(std::string const& name);
protected:
    static void registerType(std::string const& name, shapeType type);
private:
    virtual Shape doApply(ImageT const& image, int x, int y, PSF const* psf, double background) const = 0;

    static std::map<std::string, shapeType>* _shapeTypes;
};

template<typename ImageT>
measureShape<ImageT>* createmeasureShape(std::string const& type);

}}}
#endif
