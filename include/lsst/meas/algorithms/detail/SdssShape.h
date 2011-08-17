#if !defined(LSST_MEAS_ALGORITHMS_DETAIL_H)
#define LSST_MEAS_ALGORITHMS_DETAIL_H 1

namespace lsst { namespace meas { namespace algorithms { namespace detail {

class SdssShapeImpl {
public:
    typedef Eigen::Matrix4d Matrix4;    // type for the 4x4 covariance matrix
    
    SdssShapeImpl(double i0=NAN, double ixx=NAN, double ixy=NAN, double iyy=NAN) :
        _i0(i0),
        _x(NAN), _xErr(NAN), _y(NAN), _yErr(NAN),
        _ixx(ixx), _ixy(ixy), _iyy(iyy),
        _covar(),
        _ixy4(NAN),
        _flags(0) {
        _covar.setConstant(NAN);
    }

    void setI0(double i0) { _i0 = i0; }
    double getI0() const { return _i0; }
    double getI0Err() const { return ::sqrt(_covar(0, 0)); }

    double getX() const { return _x; }
    double getXErr() const { return _xErr; }
    void setX(double const x) { _x = x; }

    double getY() const { return _y; }
    double getYErr() const { return _yErr; }
    void setY(double const y) { _y = y; }
    
    void setIxx(double ixx) { _ixx = ixx; }
    double getIxx() const { return _ixx; }
    double getIxxErr() const { return ::sqrt(_covar(1, 1)); }

    void setIxy(double ixy) { _ixy = ixy; }
    double getIxy() const { return _ixy; }
    double getIxyErr() const { return ::sqrt(_covar(2, 2)); }

    void setIyy(double iyy) { _iyy = iyy; }
    double getIyy() const { return _iyy; }
    double getIyyErr() const { return ::sqrt(_covar(3, 3)); }

    void setIxy4(double ixy4) { _ixy4 = ixy4; }
    double getIxy4() const { return _ixy4; }

    void setFlags(int flags) { _flags = flags; }
    int getFlags() const { return _flags; }

    void setCovar(Matrix4 covar) { _covar = covar; }
    const Matrix4& getCovar() const { return _covar; }
    
#if !defined(SWIG)                      // XXXX
    double getE1() const;
    double getE1Err() const;
    double getE2() const;
    double getE2Err() const;
    double getE1E2Err() const;
    double getRms() const;
    double getRmsErr() const;
#endif

#ifndef SWIG
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#endif

private:
    double _i0;                         // 0-th moment

    double _x, _xErr, _y, _yErr;        // <x>, <y> and errors

    double _ixx, _ixy, _iyy;            // <xx> <xy> <yy>
    Matrix4 _covar;                     // covariance matrix for (_i0, _ixx, _ixy, _iyy)
    double _ixy4;                       // 4th moment used for shear calibration
    int _flags;                         // flags describing processing
};

template<typename ImageT>
bool getAdaptiveMoments(ImageT const& image, double bkgd, double xcen, double ycen, double shiftmax,
                        detail::SdssShapeImpl *shape);
                
}}}}

#endif

