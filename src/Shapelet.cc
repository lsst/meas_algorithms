// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#include "Eigen/Core"

#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/deprecated.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"
#include "lsst/meas/algorithms/shapelet/Pixel.h"
#include "lsst/meas/algorithms/shapelet/Bounds.h"
#include "lsst/meas/algorithms/shapelet/MyMatrix.h"
#include "lsst/meas/algorithms/shapelet/Ellipse.h"

namespace lsst {
namespace meas {
namespace algorithms {

    // All of the functionality is imported from shapelet::BVec
    // Just repeat the constructors, destructors, and op=
    class ShapeletImpl : public shapelet::BVec
    {
    public :

        typedef Shapelet::ShapeletVector ShapeletVector;
        typedef Shapelet::ShapeletCovariance ShapeletCovariance;
        typedef shapelet::BVec BVec;

        ShapeletImpl(int order, double sigma) : 
            BVec(order,sigma), _cov()
        {}

        ShapeletImpl(int order, double sigma, const ShapeletVector& vector) : 
            BVec(order,sigma,vector), _cov()
        {}

        ShapeletImpl(int order, double sigma, const ShapeletVector& vector, const ShapeletCovariance& cov) : 
            BVec(order,sigma,vector), _cov(new ShapeletCovariance(cov))
        {}

        ShapeletImpl(const BVec& rhs) : 
            BVec(rhs), _cov()
        {}

        ShapeletImpl(const ShapeletImpl& rhs) : 
            BVec(rhs)
        {
            if (rhs._cov.get()) _cov.reset(new ShapeletCovariance(*rhs._cov));
        }

        ~ShapeletImpl() {}

        using BVec::size;

        void operator=(const ShapeletImpl& rhs) 
        {
            BVec::operator=(rhs);
            if (rhs._cov.get()) {
                if (_cov.get() && 
                    _cov->TMV_colsize() == rhs._cov->TMV_colsize() &&
                    _cov->TMV_rowsize() == rhs._cov->TMV_rowsize()) {
                    *_cov = *rhs._cov;
                } else {
                    _cov.reset(new ShapeletCovariance(*rhs._cov));
                }
            } else {
                _cov.reset();
            }
        }

        boost::shared_ptr<ShapeletCovariance>& getCovariance()
        { return _cov; }

        std::complex<double> getPQ(int p, int q)
        {
            using shapelet::DVector;

            if (p < q) return std::conj(getPQ(q,p));
            else {
                const DVector& v = BVec::getValues();
                if (p == q) {
                    int k = p*(2*p+3);
                    return v(k);
                } else {
                    int k = (p+q)*(p+q+1)/2 + 2*q;
                    return std::complex<double>(v(k),v(k+1));
                }
            }
        }

        double evaluateAt(double x, double y)
        {
            // TODO: This is not efficient, but the functionality
            // is required for Kernel. 
            //
            // It is much more efficient to run makePsi on many
            // positions at once.  Furthermore, makePsi is optimized
            // to be run on many positions at the same time, which
            // involved some coding choices that probably make it even 
            // less efficient for a single position.
            //
            // However, if this is a significant time component,
            // then a more efficient version that is intended for
            // a single position could be written.
            using shapelet::CDVector;
            using shapelet::DMatrix;
            using shapelet::makePsi;

            std::complex<double> z(x,y);
            z /= getSigma();
            DMatrix psi(1,BVec::size());
            CDVector zList(1);
            zList(0) = z;
            makePsi(psi,TMV_vview(zList),BVec::getOrder(),0);
            return (psi * getValues())(0);
        }

    private :
        boost::shared_ptr<ShapeletCovariance> _cov;
    };

    Shapelet::Shapelet(int order, double sigma) :
        pImpl(new ShapeletImpl(order,sigma)) 
    {}

    Shapelet::Shapelet(int order, double sigma, const ShapeletVector& vector) :
        pImpl(new ShapeletImpl(order,sigma,vector)) 
    {}

    Shapelet::Shapelet(int order, double sigma, const ShapeletVector& vector, const ShapeletCovariance& cov) :
        pImpl(new ShapeletImpl(order,sigma,vector,cov)) 
    {}

    Shapelet::Shapelet(const Shapelet& rhs) :
        pImpl(new ShapeletImpl(*rhs.pImpl))
    {}

    Shapelet& Shapelet::operator=(const Shapelet& rhs)
    {
        *pImpl = *rhs.pImpl;
        return *this;
    }

    int Shapelet::getOrder() const 
    { return pImpl->getOrder(); }

    double Shapelet::getSigma() const 
    { return pImpl->getSigma(); }

    int Shapelet::size() const
    { return pImpl->size(); }

    const Shapelet::ShapeletVector& Shapelet::getValues() const 
    { return pImpl->getValues(); }

    bool Shapelet::hasCovariance() const 
    { return pImpl->getCovariance().get(); }

    boost::shared_ptr<const Shapelet::ShapeletCovariance> Shapelet::getCovariance() const 
    { return pImpl->getCovariance(); }

    void Shapelet::setSigma(double sigma)
    { return pImpl->setSigma(sigma); }

    std::complex<double> Shapelet::getPQ(int p, int q)
    { return pImpl->getPQ(p,q); }

    double Shapelet::evaluateAt(const PointD& pos)
    { return evaluateAt(pos.getX(),pos.getY()); }

    double Shapelet::evaluateAt(double x, double y)
    {
        // TODO: MJ - I couldn't reproduce the crash with g++ 4.3.3.
        //       Need to try this on lsst computers I guess.
#if 0                                   // crashes g++ 4.3.5
        return pImpl->evaluateAt(x,y);
#else
        // MJ - This is definitely wrong!  It starts an infinite loop!
        return evaluateAt(lsst::afw::geom::makePointD(x, y));
#endif
    }

    Eigen::Matrix2d getJacobian(
        const lsst::afw::image::Wcs& wcs,
        const lsst::afw::geom::PointD& pos)
    {
        lsst::afw::geom::AffineTransform localTransform = 
            wcs.linearizePixelToSky(pos,lsst::afw::coord::DEGREES);

        // J = ( du/dx  du/dy )
        //     ( dv/dx  dv/dy )
        // where (u,v) are sky coordinates, and (x,y) are chip coordinates.
        Eigen::Matrix2d J;
        // Answer comes out in degrees for u,v.  *3600 to get arcsec.
        J(0,0) = localTransform.getMatrix()(0,0) * 3600.;
        J(0,1) = localTransform.getMatrix()(0,1) * 3600.;
        J(1,0) = localTransform.getMatrix()(1,0) * 3600.;
        J(1,1) = localTransform.getMatrix()(1,1) * 3600.;
        return J;
    }

    // Can't use the version in Pixel.h, since we need to use the 
    // LSST Image and Wcs objects, rather than my Image and Trasnsformation.
    static void getPixList(
        shapelet::PixelList& pix,
        const Shapelet::Source& source, 
        const Shapelet::PointD& cen,
        double aperture,
        const Shapelet::Exposure& exposure,
        lsst::afw::image::MaskPixel okmask)
    {
        using shapelet::Pixel;
        using shapelet::PixelList;
        using lsst::afw::geom::PointD;

        Shapelet::Exposure::MaskedImageT const maskedImage = exposure.getMaskedImage();
        Shapelet::Exposure::MaskedImageT::Image::ConstPtr imagePtr = maskedImage.getImage();
        Shapelet::Exposure::MaskedImageT::Mask::ConstPtr maskPtr = maskedImage.getMask();
        Shapelet::Exposure::MaskedImageT::Variance::ConstPtr variancePtr = maskedImage.getVariance();

        PointD pos = lsst::afw::geom::makePointD(source.getXAstrom(),source.getYAstrom());
        Eigen::Matrix2d J = getJacobian(*(exposure.getWcs()), pos);
        
        double det = std::abs(J.determinant());
        double pixScale = sqrt(det); // arcsec / pixel
        xdbg<<"pixscale = "<<pixScale<<std::endl;

        // xAp,yAp are the maximum deviation from the center in x,y
        // such that u^2+v^2 = aperture^2
        double xAp = aperture / det * 
            sqrt(J(0,0)*J(0,0) + J(0,1)*J(0,1));
        double yAp = aperture / det * 
            sqrt(J(1,0)*J(1,0) + J(1,1)*J(1,1));
        xdbg<<"aperture = "<<aperture<<std::endl;
        xdbg<<"xap = "<<xAp<<", yap = "<<yAp<<std::endl;

        double xCen = cen.getX();
        double yCen = cen.getY();
        xdbg<<"cen = "<<xCen<<"  "<<yCen<<std::endl;

        // Find the square range bounding the aperture
        int i1 = int(floor(xCen-xAp));
        int i2 = int(ceil(xCen+xAp));
        int j1 = int(floor(yCen-yAp));
        int j2 = int(ceil(yCen+yAp));
        xdbg<<"i1,i2,j1,j2 = "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

        // Keep the range used withing the borders of the image.
        // TODO: Need to check what the LSST conventions are for x,y in
        // image(x,y).
        // Are they always relative to whatever the lower left position is?
        // Or relative to the same (0,0) that makes getX0 and getY0
        // possibly non-zero.
        int xMin = exposure.getX0();
        int yMin = exposure.getY0();
        int xMax = xMin + exposure.getWidth(); // no image->getX1() method?
        int yMax = yMin + exposure.getHeight();
        xdbg<<"xMin, yMin = "<<xMin<<"  "<<yMin<<std::endl;
        if (i1 < xMin) { i1 = xMin; }
        if (i2 > xMax) { i2 = xMax; }
        if (j1 < yMin) { j1 = yMin; }
        if (j2 > yMax) { j2 = yMax; }
        xdbg<<"i1,i2,j1,j2 => "<<i1<<','<<i2<<','<<j1<<','<<j2<<std::endl;

        // Do this next loop in two passes.  First figure out which 
        // pixels we want to use.  Then we can resize pix to the full size
        // we will need, and go back through and enter the pixels.
        // This saves us a lot of resizing calls in vector, which are
        // both slow and can fragment the memory.
        xdbg<<"nx = "<<i2-i1+1<<std::endl;
        xdbg<<"ny = "<<j2-j1+1<<std::endl;
        Assert(i2-i1+1 >= 0);
        Assert(j2-j1+1 >= 0);
        std::vector<std::vector<bool> > shouldUsePix(
            i2-i1+1,std::vector<bool>(j2-j1+1,false));
        int nPix = 0;

        const double apsq = aperture*aperture;

        double chipX = i1-xCen;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = j1-yCen;
            double u = J(0,0)*chipX+J(0,1)*chipY;
            double v = J(1,0)*chipX+J(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=J(0,1),v+=J(1,1)) {
                // u,v are in arcsec
                double rsq = u*u + v*v;
                if ( ((*maskPtr)(i,j) & ~okmask) &&
                     (rsq <= apsq) ) {
                    shouldUsePix[i-i1][j-j1] = true;
                    ++nPix;
                }
            }
        }

        xdbg<<"npix = "<<nPix<<std::endl;
        pix.resize(nPix);

        xdbg<<"pixlist size = "<<nPix<<" = "<<nPix*sizeof(Pixel)<<" bytes = "
            <<nPix*sizeof(Pixel)/1024.<<" KB\n";

        // Now the real loop that stores the flux values.
        double sky = source.getSky();

        int k=0;
        chipX = i1-xCen;
        for(int i=i1;i<=i2;++i,chipX+=1.) {
            double chipY = j1-yCen;
            double u = J(0,0)*chipX+J(0,1)*chipY;
            double v = J(1,0)*chipX+J(1,1)*chipY;
            for(int j=j1;j<=j2;++j,u+=J(0,1),v+=J(1,1)) {
                if (shouldUsePix[i-i1][j-j1]) {
                    double flux = (*imagePtr)(i,j)-sky;
                    double variance = (*variancePtr)(i,j);
                    if (variance > 0.0) {
                        double inverseSigma = sqrt(1.0/variance);
                        Assert(k < int(pix.size()));
                        pix[k++] = Pixel(u,v,flux,inverseSigma);
                    }
                }
            }
        }
        Assert(k <= int(pix.size()));
        // Not necessarily k == pix.size() 
        // because we skip pixels with 0.0 variance
        pix.resize(k); // clear off the extras, if any.
        Assert(k == int(pix.size()));
        nPix = pix.size(); // may have changed.
        xdbg<<"npix => "<<nPix<<std::endl;
    }

    bool Shapelet::measureFromImage(
        const Source& source, const PointD& pos,
        bool isCentroidFixed, bool isSigmaFixed, double aperture,
        const Exposure& exposure,
        const lsst::afw::image::MaskPixel okmask)
    {
//        const MaskedImage& image, const Wcs& wcs, 
        using shapelet::Ellipse;
        using shapelet::PixelList;

        std::vector<PixelList> pix(1);
        // Fill PixelList with pixel data around position pos:
        getPixList(pix[0], source, pos, aperture, exposure, okmask);

        double sigma = pImpl->getSigma();
        Ellipse ell;
        ell.fixGam();
        if (isCentroidFixed) ell.fixCen();
        else {
            // Initial crude estimates to get close to the right value
            // in case we have a poor starting point.
            // TODO: These might not be necessary for LSST.  
            // Should compare speed with and without this step.
            ell.peakCentroid(pix[0],aperture/3.);
            ell.crudeMeasure(pix[0],sigma);
        }
        if (isSigmaFixed) ell.fixMu();
        long flag = 0;
        if (!isCentroidFixed || !isSigmaFixed) {
            if (!ell.measure(pix,2,sigma,true,flag)) {
                return false;
            }
            if (flag) return false;
            if (!isSigmaFixed) {
                double mu = ell.getMu();
                sigma *= exp(mu);
                dbg<<"sigma = "<<sigma<<std::endl;
                Assert(sigma > 0);
                pImpl->setSigma(sigma);
            }
        }
        ell.measureShapelet(pix,*pImpl,getOrder(),pImpl->getCovariance().get());
        return true;
    }

    Shapelet::Shapelet(const shapelet::BVec& bvec) : pImpl(new ShapeletImpl(bvec)) {}

    const shapelet::BVec& Shapelet::viewAsBVec() const { return *pImpl; }
    shapelet::BVec& Shapelet::viewAsBVec() { return *pImpl; }

}}}

