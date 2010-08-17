
#include "lsst/meas/algorithms/shapelet/Ellipse.h"
#include "lsst/meas/algorithms/shapelet/EllipseSolver.h"
#include <cmath>
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include <fstream>
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"

#define N_FLUX_ATTEMPTS 0
#define MAXITER 4

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    bool Ellipse::measure(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>& psf,
        int order, double sigma, bool shouldUseInteg, long& flag,
        DMatrix* cov, BVec* bRet, DMatrix* bCov)
    { 
        return doMeasure(pix,&psf,order,sigma,shouldUseInteg,flag,cov,bRet,bCov); 
    }

    bool Ellipse::measure(
        const std::vector<PixelList>& pix,
        int order, double sigma, bool shouldUseInteg, long& flag, 
        DMatrix* cov, BVec* bRet, DMatrix* bCov)
    {
        return doMeasure(pix,0,order,sigma,shouldUseInteg,flag,cov,bRet,bCov);
    }

    void Ellipse::doMeasureShapelet(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>* psf, BVec& b, int order, DMatrix* bCov) const
    {
        xdbg<<"Start MeasureShapelet: order = "<<order<<std::endl;
        xdbg<<"b.order, sigma = "<<b.getOrder()<<", "<<b.getSigma()<<std::endl;
        //xdbg<<"el = "<<*this<<std::endl;

        // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
        // ( v )                         ( -g2   1+g1 ) ( v-vc )
        // 
        // z' = u' + I v' 
        //    = exp(-mu)/sqrt(1-gsq)
        //             [ (1-g1)(u-uc)-g2(v-vc)-Ig2(u-uc)+I(1+g1)(v-vc) ]
        //    = exp(-mu)/sqrt(1-gsq) [ z-zc - g1(z-zc)* -Ig2(z-zc)* ]
        //    = exp(-mu)/sqrt(1-gsq) ( z-zc - g (z-zc)* )

        double sigma = b.getSigma();
        double gsq = std::norm(_gamma);
        Assert(gsq < 1.);

        double m = exp(-_mu)/sqrt(1.-gsq);

        int nTot = 0;
        const int nExp = pix.size();
        for(int i=0;i<nExp;++i) nTot += pix[i].size();
        //dbg<<"ntot = "<<nTot<<" in "<<nExp<<" images\n";

        DVector I(nTot);
        DVector W(nTot);
        CDVector Z(nTot);

        for(int k=0,n=0;k<nExp;++k) {
            double sigma_obs = 
                psf ?
                sqrt(pow(sigma,2)+_fPsf*pow((*psf)[k].getSigma(),2)) :
                sigma;

            const int nPix = pix[k].size();
            for(int i=0;i<nPix;++i,++n) {
                I(n) = pix[k][i].getFlux()*pix[k][i].getInverseSigma();
                W(n) = pix[k][i].getInverseSigma();
                std::complex<double> z1 = pix[k][i].getPos();
                std::complex<double> z2 = m*(z1-_cen) - m*_gamma*conj(z1-_cen);
                Z(n) = z2 / sigma_obs;
            }
        }

        const int bSize = (order+1)*(order+2)/2;

        DMatrix A(nTot,bSize);
        Assert(nTot >= bSize); // Should have been addressed by calling routine.
        makePsi(A,TMV_vview(Z),order,&W);

        if (psf) {
            for(int k=0,n=0,nx;k<nExp;++k,n=nx) {
                const int psfSize = (*psf)[k].size();
                BVec newPsf = (*psf)[k];
                DMatrix S(psfSize,psfSize);
                calculateGTransform(_gamma,newPsf.getOrder(),S);
                newPsf.vec() = S * (*psf)[k].vec();
                DMatrix D(psfSize,psfSize);
                calculateMuTransform(_mu,newPsf.getOrder(),D);
                newPsf.vec() = D * newPsf.vec();
                DMatrix C(bSize,bSize);
                calculatePsfConvolve(newPsf,order,b.getSigma(),C);
                const int nPix = pix[k].size();
                nx = n+nPix;
                TMV_rowRange(A,n,nx) *= C;
            }
        }
#ifdef USE_TMV
        A.divideUsing(tmv::SV);
        A.saveDiv();
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
        A.svd().thresh(sqrtEps);
        dbg<<"For MeasureShapelet: svd = "<<A.svd().getS().diag()<<std::endl;
        if (A.svd().getKMax() < A.rowsize())
            dbg<<"Omitting last "<<A.rowsize()-A.svd().getKMax()<<" singular values\n";
        b.vec().subVector(0,bSize) = I/A;
        if (bCov) {
            A.makeInverseATA(*bCov);
        }
#else
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
        Eigen::SVD<DMatrix> svd = A.svd().sort();
        const DMatrix& svd_u = svd.matrixU();
        const DVector& svd_s = svd.singularValues();
        const DMatrix& svd_v = svd.matrixV();
        double max = svd_s(0);
        int kmax = 0;
        double thresh = sqrtEps * max;
        while (kmax < svd_s.size() && svd_s(kmax) > thresh) ++kmax;
        dbg<<"For MeasureShapelet: svd = "<<svd_s.transpose()<<std::endl;
        if (kmax < svd_s.size())
            dbg<<"Omitting last "<<svd_s.size()-kmax<<" singular values\n";
        // USVtb = I
        // b = VS^-1UtI
        DVector temp = TMV_colRange(svd_u,0,kmax).transpose() * I;
        //dbg<<"after temp = Ut I\n";
        temp = svd_s.TMV_subVector(0,kmax).cwise().inverse().asDiagonal() * temp;
        //dbg<<"after temp = S^-1 * temp\n";
        b.vec().TMV_subVector(0,bSize) = TMV_colRange(svd_v,0,kmax) * temp;
        //dbg<<"after b = v * temp\n";
        if (bCov) {
            //dbg<<"make bCov:\n";
            // (AtA)^-1 = (VSUt USVt)^-1 = (V S^2 Vt)^-1 = V S^-2 Vt
            DMatrix temp2 = 
                svd_s.TMV_subVector(0,kmax).cwise().square().inverse().asDiagonal() *
                TMV_colRange(svd_v,0,kmax).transpose();
            //dbg<<"after temp2 = s^-2 * v\n";
            *bCov = TMV_colRange(svd_v,0,kmax).transpose() * temp2;
            //dbg<<"after bCov = vt * temp2\n";
        }
#endif
        if (order < b.getOrder()) {
            dbg<<"Need to zero the rest of b\n";
            dbg<<"bSize = "<<bSize<<"  b.size = "<<b.size()<<std::endl;
            // Zero out the rest of the shapelet vector:
            b.vec().TMV_subVector(bSize,b.size()).setZero();
        }
        //dbg<<"Done measure Shapelet\n";
    }

    void Ellipse::measureShapelet(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>& psf, BVec& b, int order, DMatrix* bCov) const
    { doMeasureShapelet(pix,&psf,b,order,bCov); }

    void Ellipse::measureShapelet(
        const std::vector<PixelList>& pix, BVec& b, int order, DMatrix* bCov) const
    { doMeasureShapelet(pix,0,b,order,bCov); }

}}}}
