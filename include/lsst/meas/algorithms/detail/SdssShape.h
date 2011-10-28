#if !defined(LSST_MEAS_ALGORITHMS_DETAIL_H)
#define LSST_MEAS_ALGORITHMS_DETAIL_H 1

#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/detection/Shape.h"

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

    explicit SdssShapeImpl(lsst::afw::detection::Shape const& shape) :
        _i0(NAN),
        _x(shape.getX()), _xErr(shape.getXErr()), _y(shape.getY()), _yErr(shape.getYErr()),
        _ixx(shape.getIxx()), _ixy(shape.getIxy()), _iyy(shape.getIyy()),
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

    /// Return multiplier that transforms I0 to a total flux
    double getFluxScale() const {
        /*
         * The shape is an ellipse that's axis-aligned in (u, v) [<uv> = 0] after rotation by theta:
         * <x^2> + <y^2> = <u^2> + <v^2>
         * <x^2> - <y^2> = cos(2 theta)*(<u^2> - <v^2>)
         * 2*<xy>        = sin(2 theta)*(<u^2> - <v^2>)
         */
        double const Mxx = getIxx(); // <x^2>
        double const Mxy = getIxy(); // <xy>
        double const Myy = getIyy(); // <y^2>
        
        double const Muu_p_Mvv = Mxx + Myy;                             // <u^2> + <v^2>
        double const Muu_m_Mvv = ::sqrt(::pow(Mxx - Myy, 2) + 4*::pow(Mxy, 2)); // <u^2> - <v^2>
        double const Muu = 0.5*(Muu_p_Mvv + Muu_m_Mvv);
        double const Mvv = 0.5*(Muu_p_Mvv - Muu_m_Mvv);
        
        return afwGeom::TWOPI * ::sqrt(Muu*Mvv);
    }

    PTR(SdssShapeImpl) transform(afwGeom::AffineTransform const& trans) const {
        PTR(SdssShapeImpl) shape = boost::make_shared<SdssShapeImpl>();

        namespace afwEll = lsst::afw::geom::ellipses;

        double invJacobian = 1.0 / trans.getLinear().computeDeterminant();
        afwGeom::Point2D const& center = trans(afwGeom::Point2D(_x, _y));
        afwEll::Quadrupole const& moments(*afwEll::Quadrupole(_ixx, _iyy, _ixy).
                                          transform(trans.getLinear()).copy());

        shape->setI0(_i0 * invJacobian);
        shape->setX(center.getX());
        shape->setY(center.getY());
        shape->setIxx(moments.getIXX());
        shape->setIxy(moments.getIXY());
        shape->setIyy(moments.getIYY());
        // XXX errors?
        // XXX covar?
        // XXX ixy4?
        shape->setFlags(_flags);

        return shape;
    }

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

template<typename ImageT>
std::pair<double, double>
getFixedMomentsFlux(ImageT const& mimage, double bkgd, double xcen, double ycen,
                    detail::SdssShapeImpl const& shape);

}}}}

#endif

