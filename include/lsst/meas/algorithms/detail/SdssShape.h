// -*- lsst-c++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_DETAIL_H)
#define LSST_MEAS_ALGORITHMS_DETAIL_H 1

#include <bitset>

#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/Angle.h"

namespace lsst { namespace meas { namespace algorithms { namespace detail {

int const SDSS_SHAPE_MAX_ITER = 100;  // Default maximum number of iterations
float const SDSS_SHAPE_TOL1 = 1.0e-5; // Default convergence tolerance for e1,e2
float const SDSS_SHAPE_TOL2 = 1.0e-4; // Default convergence tolerance for FWHM

class SdssShapeImpl {
public:
    typedef Eigen::Matrix<double,4,4,Eigen::DontAlign> Matrix4;    // type for the 4x4 covariance matrix

    enum Flag {
        UNWEIGHTED_BAD=0,
        UNWEIGHTED,
        SHIFT,
        MAXITER,
        N_FLAGS
    };
    
    typedef std::bitset<N_FLAGS> FlagSet;

    SdssShapeImpl(double i0=NAN, double ixx=NAN, double ixy=NAN, double iyy=NAN) :
        _i0(i0),
        _x(NAN), _xErr(NAN), _y(NAN), _yErr(NAN),
        _ixx(ixx), _ixy(ixy), _iyy(iyy),
        _covar(),
        _ixy4(NAN),
        _flags()
    {
        _covar.setConstant(NAN);
    }
    
    explicit SdssShapeImpl(
        afw::geom::Point2D const & centroid,
        afw::geom::ellipses::Quadrupole const & shape
    ) : _i0(NAN),
        _x(centroid.getX()), _xErr(NAN), _y(centroid.getY()), _yErr(NAN),
        _ixx(shape.getIxx()), _ixy(shape.getIxy()), _iyy(shape.getIyy()),
        _covar(),
        _ixy4(NAN),
        _flags()
    {
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

    void setFlag(Flag flag) { _flags.set(flag); }
    void resetFlag(Flag flag) { _flags.reset(flag); }
    bool getFlag(Flag flag) const { return _flags.test(flag); }

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
        
        return lsst::afw::geom::TWOPI * ::sqrt(Muu*Mvv);
    }

    PTR(SdssShapeImpl) transform(lsst::afw::geom::AffineTransform const& trans) const {
        PTR(SdssShapeImpl) shape = boost::make_shared<SdssShapeImpl>();

        double invJacobian = 1.0 / trans.getLinear().computeDeterminant();
        lsst::afw::geom::Point2D const& center = trans(lsst::afw::geom::Point2D(_x, _y));
        lsst::afw::geom::ellipses::Quadrupole const& moments(
            *lsst::afw::geom::ellipses::Quadrupole(_ixx, _iyy, _ixy).transform(trans.getLinear()).copy());

        shape->setI0(_i0 * invJacobian);
        shape->setX(center.getX());
        shape->setY(center.getY());
        shape->setIxx(moments.getIxx());
        shape->setIxy(moments.getIxy());
        shape->setIyy(moments.getIyy());
        // XXX errors?
        // XXX covar?
        // XXX ixy4?

        shape->_flags = _flags;

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
    FlagSet _flags;                     // flags describing processing
};

/// Calculate adaptive moments from an image
///
/// The moments are measured iteratively with a Gaussian window with
/// width equal to the second moments from the previous iteration.
template<typename ImageT>
bool getAdaptiveMoments(
    ImageT const& mimage,               ///< the data to process
    double bkgd,                        ///< background level
    double xcen,                        ///< x-centre of object
    double ycen,                        ///< y-centre of object
    double shiftmax,                    ///< max allowed centroid shift
    detail::SdssShapeImpl *shape,       ///< a place to store desired data
    int maxIter=SDSS_SHAPE_MAX_ITER,    ///< Maximum number of iterations
    float tol1=SDSS_SHAPE_TOL1,         ///< Convergence tolerance for e1,e2
    float tol2=SDSS_SHAPE_TOL2,         ///< Convergence tolerance for FWHM
    bool negative=false                 ///< measured object is a negative detection
                                        // (e.g. from difference imaging)     
    );

template<typename ImageT>
std::pair<double, double>
getFixedMomentsFlux(ImageT const& mimage, double bkgd, double xcen, double ycen,
                    detail::SdssShapeImpl const& shape);

}}}}

#endif

