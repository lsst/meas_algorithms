// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef LSST_MEAS_ALGORITHMS_ImageMoments_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_ImageMoments_h_INCLUDED

#include "lsst/afw/geom/ellipses.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief Class for computing non-adaptive Gaussian moments of images.
 *
 *  Note that all methods of this class return results in a coordinate system shifted
 *  from the one of the image.  This is unavoidable for raw moments, and I opted for
 *  consistency between the various methods rather than using the image coordinate
 *  system only when possible.  Thus the user will generally want to add the fidicial
 *  center passed to the algorithm to the returned centroid.
 */
class ImageMoments {
public:

    /**
     *  @brief Integer constants used to index "raw" moment vectors and matrices.
     *
     *  For weight function @f$w@f$ and data @f$p@f$, the "raw" moments @f$Q@f$ are defined as:
     *  @f{eqnarray*}{
     *    Q_0 &=& \sum_n w(x_n, y_n) p_n \\
     *    Q_{xx} &=& \sum_n w(x_n, y_n) x_n^2 p_n \\
     *    Q_{yy} &=& \sum_n w(x_n, y_n) y_n^2 p_n \\
     *    Q_{xy} &=& \sum_n w(x_n, y_n) x_n y_n p_n \\
     *    Q_x &=& \sum_n w(x_n, y_n) x_n p_n \\
     *    Q_y &=& \sum_n w(x_n, y_n) y_n p_n
     *  @f}
     *
     *  Note that we order the second moments before the first moments, for consistency with
     *  afw::geom::ellipses::Ellipse::ParameterEnum.
     */
    enum RawParameterEnum { Q0=0, QXX, QYY, QXY, QX, QY, N_Q };

    /**
     *  @brief Integer constants used to index ellipse moment vectors and matrices.
     *
     *  @f{eqnarray*}{
     *    M_{xx} &=& Q_{xx} / Q_0 - Q_x^2 \\
     *    M_{xx} &=& Q_{yy} / Q_0 - Q_y^2 \\
     *    M_{xx} &=& Q_{xy} / Q_0 - Q_x Q_y \\
     *    M_x &=& Q_x / Q_0 \\
     *    M_y &=& Q_y / Q_0
     *  @f}
     *
     *  Note that we order the second moments before the first moments, for consistency with
     *  afw::geom::ellipses::Ellipse::ParameterEnum.
     */
    enum EllipseParameterEnum { MXX=0, MYY, MXY, MX, MY, N_M };

    typedef Eigen::Matrix<double,N_Q,1> VectorQ;    ///< Vector type for raw moments
    typedef Eigen::Matrix<double,N_Q,N_Q> MatrixQ;  ///< Matrix type for raw moments
    typedef Eigen::Matrix<double,N_M,1> VectorM;    ///< Vector type for ellipse moments
    typedef Eigen::Matrix<double,N_M,N_M> MatrixM;  ///< Matrix type for ellipse moments
    typedef Eigen::Matrix<double,N_M,N_Q> MatrixMQ; ///< Matrix type for mapping raw to ellipse moments

#ifndef SWIG
    class RawResult;
    class EllipseResult;
#else
    typedef ImageMomentsRawResult RawResult;
    typedef ImageMomentsEllipseResult EllipseResult;
#endif

    /**
     *  Initialize with the given configuration values.
     *
     *  @param[in] nSigmaRegion       Size of pixel region to use when computing moments, in units of the
     *                                weight function RMS size.
     *  @param[in] useApproximateExp  Whether to use fast approximate exponentials when evaluating the
     *                                weight function.
     */
    ImageMoments(double nSigmaRegion=3.0, bool useApproximateExp=false);

    /**
     *  Measure the raw moments of an Image or MaskedImage.
     *
     *  See RawResult and RawParameterEnum for a complete description of the outputs.
     *
     *  @param[in]  image    Image or MaskedImage on which to measure.  If a MaskedImage, the covariance
     *                       matrix of the moments will be computed as well.
     *  @param[in]  weight   Ellipse that defines the Gaussian to be used as a weight function.
     *  @param[in]  origin   Point that determines the origin of the coordinate system for the measured
     *                       moments, and also sets the center of the weight function.
     *
     *  Instantiated for ImageT = (float, double) x (Image, MaskedImage)
     */
    template <typename ImageT>
    RawResult measureRaw(
        ImageT const & image,
        afw::geom::ellipses::Quadrupole const & weight,
        afw::geom::Point2D const & origin
    ) const;

    /**
     *  Measure the ellipse moments of an Image or MaskedImage.
     *
     *  See EllipseResult and EllipseParameterEnum for a complete description of the outputs.
     *
     *  @param[in]  image    Image or MaskedImage on which to measure.  If a MaskedImage, the covariance
     *                       matrix of the moments will be computed as well.
     *  @param[in]  weight   Ellipse that defines the Gaussian to be used as a weight function.
     *  @param[in]  origin   Point that determines the origin of the coordinate system for the measured
     *                       moments, and also sets the center of the weight function.
     *
     *  Instantiated for ImageT = (float, double) x (Image, MaskedImage)
     */
    template <typename ImageT>
    EllipseResult measureEllipse(
        ImageT const & image,
        afw::geom::ellipses::Quadrupole const & weight,
        afw::geom::Point2D const & origin
    ) const;

    /**
     *  @brief Convert linear raw moments into ellipse moments, and optionally compute the derivative
     *         of the conversion.
     *
     *  @note This function is mainly intended for internal use, and is only exposed publically
     *        so it can be unit-tested in Python.
     */
    static VectorM convertRawMoments(VectorQ const & q, MatrixMQ * dm_dq = 0);

    /**
     *  @brief Correct moments measured with a Gaussian weight function by assuming the data was also
     *         an elliptical Gaussian, and optionally compute the derivative of the correction.
     *
     *  @note This function is mainly intended for internal use, and is only exposed publically
     *        so it can be unit-tested in Python.
     *
     *  If we naively measure Gaussian-weighted moments, we'll measure the moments of the product
     *  of the weight function and the data.  What we want is the moments of the data, as if we
     *  had measured them with no weight function (but without sacrificing the S/N benefit that
     *  comes from using a weight function).  To do that, we assume the data is also an elliptical
     *  Gaussian, and "divide" the weight function from the measured moments to compute it.
     *
     *  If @f$W@f$ and @f$M@f$ are the quadruple matrices of the weight function and measurement,
     *  and @f$\eta@f$ is the measured centroid (we work in a coordinate system where the weight
     *  function is centered at the origin), then the corrected quadrupole matrix @f$C@f$ and
     *  centroid are @f$\nu@f$ are:
     *  @f{eqnarray*}{
     *    C &=& \left(M^{-1} - W^{-1}\right)^{-1} \\
     *    \nu &=& C M^{-1} \eta
     *  @f}
     */
    static VectorM correctWeightedMoments(
        afw::geom::ellipses::Quadrupole const & weight,
        VectorM const & m,
        MatrixM * dc_dm = 0
    );

private:

    struct RawMomentAccumulator;

    double _nSigmaRegion;
    bool _useApproximateExp;
};



/**
 *  @brief Struct to hold the results of measureEllipse
 */
#ifndef SWIG // trick swig into thinking this isn't an inner class
struct ImageMoments::EllipseResult
#else
struct ImageMomentsEllipseResult
#endif
{

    /// Return the 2nd moments as a Quadrupole object
    afw::geom::ellipses::Quadrupole getQuadrupole() const {
        return afw::geom::ellipses::Quadrupole(moments[MXX], moments[MYY], moments[MXY]);
    }

    /// Return the 1st moments as a Point object <em>relative to input origin</em>
    afw::geom::Point2D getCentroid() const { return afw::geom::Point2D(moments[MX], moments[MY]); }

    /// Return the moments as an Ellipse object
    afw::geom::ellipses::Ellipse getEllipse() const {
        return afw::geom::ellipses::Ellipse(getQuadrupole(), getCentroid());
    }

#ifndef SWIG

    /// Vector of measured moments (see EllipseParameterEnum for element definitions)
    VectorM moments;

    /// Matrix of uncertainties (see EllipseParameterEnum for element definitions)
    PTR(MatrixM) covariance;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Eigen fixed-size arrays may require 16-bit alignment
#endif
};


#ifndef SWIG // swig can't do anything with this class, so we'll completely it in an %extend block

/**
 *  @brief Struct to hold the results of measureRaw
 */
struct ImageMoments::RawResult {
    VectorQ moments; ///< Vector of raw moments, ordered Q0, Qxx, Qyy, Qxy, Qx, Qy.

    PTR(MatrixQ) covariance; ///< Matrix of uncertainties; ordered Q0, Qxx, Qyy, Qxy, Qx, Qy.

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Eigen fixed-size arrays may require 16-bit alignment
};

#else

struct ImageMomentsRawResult {};

#endif

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_ImageMoments_h_INCLUDED
