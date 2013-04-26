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

#include "lsst/pex/exceptions.h"
#include "lsst/utils/PowFast.h"
#include "lsst/afw/geom/ellipses/PixelRegion.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/meas/algorithms/ImageMoments.h"

namespace lsst { namespace meas { namespace algorithms {

/* ---- Machinery for iterating over elliptical regions ----------------------------------------------------
 *
 * In contrast to FootprintFunctor (which is what's used by most other measurement algorithms), this
 * machinery doesn't require a virtual function call at each pixel, and hence allows the functor to
 * be inlined.
 * Eventually, something like this will replace FootprintFunctor (see #1836)
 */

namespace {

using afw::geom::Span;

Span clipSpan(Span const & span, afw::geom::Box2I const & box) {
    if (span.getY() < box.getMinY() || span.getY() > box.getMaxY()) return Span();
    return Span(
        span.getY(), std::max(span.getMinX(), box.getMinX()), std::min(span.getMaxX(), box.getMaxX())
    );
}

template <typename Function, typename Iterator>
void iterateSpan(Function function, Iterator pixIter, Span const & span) {
    for (
        Span::Iterator pointIter = span.begin(), pointEnd = span.end();
        pointIter != pointEnd;
        ++pointIter, ++pixIter
    ) {
        boost::unwrap_ref(function)(*pointIter, *pixIter);
    }
}

template <typename Function, typename Image, typename Region>
void iterateRegion(Function function, Image const & image, Region const & region) {
    afw::geom::Box2I bbox = image.getBBox(afw::image::PARENT);
    if (bbox.contains(region.getBBox())) {
        // if the box contains the region, there's no need to check each span to make sure it's entirely
        // within the image
        for (
            typename Region::Iterator spanIter = region.begin(), spanEnd = region.end();
            spanIter != spanEnd;
            ++spanIter
        ) {
            iterateSpan(
                function,
                image.x_at(spanIter->getMinX() - image.getX0(), spanIter->getY() - image.getY0()),
                *spanIter
            );
        }
    } else {
        for (
            typename Region::Iterator spanIter = region.begin(), spanEnd = region.end();
            spanIter != spanEnd;
            ++spanIter
        ) {
            Span span = clipSpan(*spanIter, bbox);
            iterateSpan(
                function,
                image.x_at(span.getMinX() - image.getX0(), span.getY() - image.getY0()),
                span
            );
        }
    }
}

} // anonymous

/* ---- Implementation for SimpleShape algorithm ------------------------------------------------------------
 *
 * All of the pixel-based work is done by an accumulating functor that we invoke using the above machinery.
 * The rest is in the static member function measure(), which provides a way to call the algorithm
 * outside the context of the pluggable source measurement algorithm framework.
 */

namespace {

// singleton lookup-table for fast approximate exponentials
utils::PowFast const & powFast = utils::getPowFast<11>();

// metafunction to distinguish between afw::image::Image and afw::image::MaskedImage
template <typename ImageT>
struct IsMaskedImage {
    static bool const value = false;
};
template <typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
struct IsMaskedImage< afw::image::MaskedImage<ImagePixelT,MaskPixelT,VariancePixelT> > {
    static bool const value = true;
};

} // anonymous

struct ImageMoments::RawMomentAccumulator {

    template <typename PixelT>
    void update(
        PixelT const & pixel,
        VectorQ const & q
    ) {
        _result.moments += q * pixel;
    }

    template <typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
    void update(
        afw::image::pixel::Pixel<ImagePixelT,MaskPixelT,VariancePixelT> const & pixel,
        VectorQ const & q
    ) {
        _result.moments += q * pixel.image();
        _result.covariance->selfadjointView<Eigen::Lower>().rankUpdate(q, pixel.variance());
    }

    template <typename PixelT>
    void operator()(afw::geom::Point2I const & pos, PixelT const & pixel) {
        afw::geom::Extent2D d = afw::geom::Point2D(pos) - _center;
        afw::geom::Extent2D gtd = _gt(d);
        float w;
        if (_useApproximateExp) {
            w = powFast.exp(-0.5 * (gtd.getX()*gtd.getX() + gtd.getY()*gtd.getY()));
        } else {
            w = std::exp(-float(0.5 * (gtd.getX()*gtd.getX() + gtd.getY()*gtd.getY())));
        }
        _wSum += w;
        VectorQ q = VectorQ::Constant(w);
        q.coeffRef(QXX) *= d.getX() * d.getX();
        q.coeffRef(QYY) *= d.getY() * d.getY();
        q.coeffRef(QXY) *= d.getX() * d.getY();
        q.coeffRef(QX) *= d.getX();
        q.coeffRef(QY) *= d.getY();
        update(pixel, q);
    }

    template <typename ImageT>
    explicit RawMomentAccumulator(
        afw::geom::ellipses::Quadrupole const & weight,
        afw::geom::Point2D const & center,
        bool useApproximateExp,
        ImageT const *
    ) :
        _useApproximateExp(useApproximateExp),
        _wSum(0.0), _center(center),
        _gt(weight.getGridTransform())
    {
        _result.moments.setZero();
        if (IsMaskedImage<ImageT>::value) {
            _result.covariance.reset(new MatrixQ());
            _result.covariance->setZero();
        }
    }

    RawResult finish() {
        _result.moments /= _wSum;
        if (_result.covariance) {
            (*_result.covariance) /= _wSum * _wSum;
        }
        return _result;
    }

private:
    bool _useApproximateExp;
    double _wSum;
    afw::geom::Point2D _center;
    afw::geom::LinearTransform _gt;
    RawResult _result;
};

ImageMoments::ImageMoments(double nSigmaRegion, bool useApproximateExp) :
    _nSigmaRegion(nSigmaRegion), _useApproximateExp(useApproximateExp)
{}

template <typename ImageT>
ImageMoments::RawResult ImageMoments::measureRaw(
    ImageT const & image,
    afw::geom::ellipses::Quadrupole const & weight,
    afw::geom::Point2D const & center
) const {
    afw::geom::ellipses::Ellipse regionEllipse(weight, center);
    regionEllipse.getCore().scale(_nSigmaRegion);
    afw::geom::ellipses::PixelRegion region(regionEllipse);
    RawMomentAccumulator functor(weight, center, _useApproximateExp, &image);
    iterateRegion(boost::ref(functor), image, region);
    return functor.finish();
}

template <typename ImageT>
ImageMoments::EllipseResult ImageMoments::measureEllipse(
    ImageT const & image,
    afw::geom::ellipses::Quadrupole const & weight,
    afw::geom::Point2D const & center
) const {
    RawResult q = measureRaw(image, weight, center);
    EllipseResult result;
    if (IsMaskedImage<ImageT>::value) {
        MatrixMQ dm_dq;
        VectorM m = convertRawMoments(q.moments, &dm_dq);
        MatrixM mCov = dm_dq * q.covariance->selfadjointView<Eigen::Lower>() * dm_dq.adjoint();
        MatrixM dc_dm;
        result.moments = correctWeightedMoments(weight, m, &dc_dm);
        result.covariance.reset(
            new MatrixM(dc_dm * mCov.selfadjointView<Eigen::Lower>() * dc_dm.adjoint())
        );
    } else {
        VectorM m = convertRawMoments(q.moments);
        result.moments = correctWeightedMoments(weight, m);
    }
    return result;
}

ImageMoments::VectorM ImageMoments::convertRawMoments(VectorQ const & q, MatrixMQ * dm_dq) {
    VectorM m = q.segment<N_M>(1) / q.coeff(Q0);
    if (dm_dq) {
        dm_dq->block<N_M,N_M>(0,1) = MatrixM::Identity() / q.coeff(Q0);
        dm_dq->col(Q0) = -m / q.coeff(Q0);
    }
    m.coeffRef(MXX) -= m.coeff(MX) * m.coeff(MX);
    m.coeffRef(MYY) -= m.coeff(MY) * m.coeff(MY);
    m.coeffRef(MXY) -= m.coeff(MX) * m.coeff(MY);
    if (dm_dq) {
        dm_dq->coeffRef(MXX, QX) = -2.0 * m.coeff(MX) / q.coeff(Q0);
        dm_dq->coeffRef(MYY, QY) = -2.0 * m.coeff(MY) / q.coeff(Q0);
        dm_dq->coeffRef(MXY, QX) = - m.coeff(MY) / q.coeff(Q0);
        dm_dq->coeffRef(MXY, QY) = - m.coeff(MX) / q.coeff(Q0);
        double tmp = 2.0 / (q.coeff(Q0) * q.coeff(Q0) * q.coeff(Q0)); // == d(-Q_0^{-2})/d(Q_0)
        dm_dq->coeffRef(MXX, Q0) += tmp * q.coeff(QX) * q.coeff(QX);
        dm_dq->coeffRef(MYY, Q0) += tmp * q.coeff(QY) * q.coeff(QY);
        dm_dq->coeffRef(MXY, Q0) += tmp * q.coeff(QX) * q.coeff(QY);
    }
    return m;
}

ImageMoments::VectorM ImageMoments::correctWeightedMoments(
    afw::geom::ellipses::Quadrupole const & weight,
    VectorM const & m,
    MatrixM * dc_dm
) {
    Eigen::Matrix2d wMat = weight.getMatrix();
    Eigen::Vector2d mVec(m.coeff(MX), m.coeff(MY));
    Eigen::Matrix2d mMat;
    mMat <<
        m.coeff(MXX), m.coeff(MXY),
        m.coeff(MXY), m.coeff(MYY);
    if (wMat.determinant() <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            (boost::format("Weight moments matrix is singular; |W|=%f") % wMat.determinant()).str()
        );
    }
    if (mMat.determinant() <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            (boost::format("Measured moments matrix is singular; |M|=%f") % mMat.determinant()).str()
        );
    }
    Eigen::Matrix2d mInv = mMat.inverse();
    Eigen::Matrix2d cInv = mInv - wMat.inverse();
    if (cInv.determinant() <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            (boost::format("Corrected moments matrix is singular; |M|^{-1}=%f") % cInv.determinant()).str()
        );
    }
    Eigen::Matrix2d cMat = cInv.inverse();
    Eigen::Matrix2d cMat_mInv = cMat * mInv;
    Eigen::Vector2d cVec = cMat_mInv * mVec;

    if (dc_dm) {
        Eigen::Vector2d mInv_mVec = mInv * mVec;
        Eigen::Matrix2d dcMat_dmxx = cMat_mInv.col(0) * cMat_mInv.col(0).adjoint();
        Eigen::Matrix2d dcMat_dmyy = cMat_mInv.col(1) * cMat_mInv.col(1).adjoint();
        Eigen::Matrix2d dcMat_dmxy = cMat_mInv.col(0) * cMat_mInv.col(1).adjoint()
            + cMat_mInv.col(1) * cMat_mInv.col(0).adjoint();
        Eigen::Vector2d dcVec_dmxx = dcMat_dmxx * mInv_mVec - cMat_mInv.col(0) * mInv_mVec.coeff(0);
        Eigen::Vector2d dcVec_dmyy = dcMat_dmyy * mInv_mVec - cMat_mInv.col(1) * mInv_mVec.coeff(1);
        Eigen::Vector2d dcVec_dmxy = dcMat_dmxy * mInv_mVec - cMat_mInv.col(0) * mInv_mVec.coeff(1)
            - cMat_mInv.col(1) * mInv_mVec.coeff(0);
        Eigen::Matrix2d const & dcVec_dmVec = cMat_mInv;

        // calculations done - now we just need to put these elements back into the 5x5 matrix
        dc_dm->setZero();
        dc_dm->coeffRef(MXX, MXX) = dcMat_dmxx.coeff(0, 0);
        dc_dm->coeffRef(MYY, MXX) = dcMat_dmxx.coeff(1, 1);
        dc_dm->coeffRef(MXY, MXX) = dcMat_dmxx.coeff(0, 1);
        dc_dm->coeffRef(MXX, MYY) = dcMat_dmyy.coeff(0, 0);
        dc_dm->coeffRef(MYY, MYY) = dcMat_dmyy.coeff(1, 1);
        dc_dm->coeffRef(MXY, MYY) = dcMat_dmyy.coeff(0, 1);
        dc_dm->coeffRef(MXX, MXY) = dcMat_dmxy.coeff(0, 0);
        dc_dm->coeffRef(MYY, MXY) = dcMat_dmxy.coeff(1, 1);
        dc_dm->coeffRef(MXY, MXY) = dcMat_dmxy.coeff(0, 1);
        dc_dm->coeffRef(MX, MXX) = dcVec_dmxx.coeff(0);
        dc_dm->coeffRef(MY, MXX) = dcVec_dmxx.coeff(1);
        dc_dm->coeffRef(MX, MYY) = dcVec_dmyy.coeff(0);
        dc_dm->coeffRef(MY, MYY) = dcVec_dmyy.coeff(1);
        dc_dm->coeffRef(MX, MXY) = dcVec_dmxy.coeff(0);
        dc_dm->coeffRef(MY, MXY) = dcVec_dmxy.coeff(1);
        dc_dm->coeffRef(MX, MX) = dcVec_dmVec.coeff(0, 0);
        dc_dm->coeffRef(MX, MY) = dcVec_dmVec.coeff(0, 1);
        dc_dm->coeffRef(MY, MX) = dcVec_dmVec.coeff(1, 0);
        dc_dm->coeffRef(MY, MY) = dcVec_dmVec.coeff(1, 1);
    }

    VectorM c;
    c.coeffRef(MXX) = cMat.coeff(0,0);
    c.coeffRef(MYY) = cMat.coeff(1,1);
    c.coeffRef(MXY) = cMat.coeff(0,1);
    c.coeffRef(MX) = cVec.coeff(0);
    c.coeffRef(MY) = cVec.coeff(1);

    return c;
}

#define INSTANTIATE(IM)                                                 \
    template ImageMoments::RawResult ImageMoments::measureRaw(          \
        afw::image::IM const &,                                         \
        afw::geom::ellipses::Quadrupole const &,                        \
        afw::geom::Point2D const &                                      \
    ) const;                                                                  \
    template ImageMoments::EllipseResult ImageMoments::measureEllipse(  \
        afw::image::IM const &,                                         \
        afw::geom::ellipses::Quadrupole const &,                        \
        afw::geom::Point2D const &                                      \
    ) const

INSTANTIATE(Image<float>);
INSTANTIATE(Image<double>);
INSTANTIATE(MaskedImage<float>);
INSTANTIATE(MaskedImage<double>);

}}} // namespace lsst::meas::algorithms
