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

#if 0
#include <valarray>
#include <fstream>
#include <CCfits/CCfits>
#endif

#include "lsst/meas/algorithms/shapelet/FittedPsf.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include "lsst/meas/algorithms/shapelet/Function2D.h"
#include "lsst/meas/algorithms/shapelet/Legendre2D.h"
#if 0
#include "Name.h"
#include "WlVersion.h"
#include "WriteParam.h"
#endif

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    static DVector definePXY(
        int order, double x, double xMin, double xMax)
    {
        DVector temp(order+1);
        double newX = (2.*x-xMin-xMax)/(xMax-xMin);
        temp[0] = 1.;
        if(order>0) temp[1] = newX;
        for(int i=2;i<=order;++i) {
            temp[i] = ((2.*i-1.)*newX*temp[i-1] - (i-1.)*temp[i-2])/i;
        }
        return temp;
    }

#ifdef USE_TMV
    static void setPRow(
        int fitOrder, Position pos, const Bounds& bounds, DVectorView pRow)
#else
        static void setPRow(
            int fitOrder, Position pos, const Bounds& bounds, DVector& pRow)
#endif
        {
            Assert(int(pRow.size()) == (fitOrder+1)*(fitOrder+2)/2);
            DVector px = 
                definePXY(fitOrder,pos.getX(),bounds.getXMin(),bounds.getXMax());
            DVector py = 
                definePXY(fitOrder,pos.getY(),bounds.getYMin(),bounds.getYMax());
            int pq = 0;
            for(int n=0;n<=fitOrder;++n) {
                for(int p=n,q=n-p;q<=n;--p,++q) {
                    Assert(pq < int(pRow.size()));
                    pRow(pq) = px[p]*py[q];
                    ++pq;
                }
            }
            Assert(pq == int(pRow.size()));
        }

    void FittedPsf::calculate(
        const std::vector<Position>& pos,
        const std::vector<BVec>& psf,
        const std::vector<double>& nu,
        std::vector<long>& flags)
    {
        const int nStars = pos.size();
        const int psfSize = (_psfOrder+1)*(_psfOrder+2)/2;
        const double nSigmaClip = _params.read("fitpsf_nsigma_outlier",3);

        // This is an empirical fit to the chisq level that corresponds to 
        // a 3-sigma outlier for more than 1 dimension.
        // I calculated the values up to n=100 and for n>30, they form a pretty
        // good approximation to a straight line.
        // This is almost certainly wrong for nSigmaClip != 3, so if we start 
        // choosing other values for nSigmaClip, it might be worth doing this 
        // right.
        // That means calculating the 1-d critical value for the given nSigma.
        // e.g. nSigma = 3 -> alpha = P(chisq > 9) = 0.0027.
        // Then calculate the critical value of chisq for that alpha with 
        // the full degrees of freedom = psfSize
        const double chisqLevel = 0.14*psfSize + 2.13;

        const double outlierThresh = nSigmaClip * nSigmaClip * chisqLevel;
        dbg<<"outlierThresh = "<<outlierThresh<<std::endl;

        int nGoodPsf, nOutliers, dof;
        double chisq;
        do {

            // Calculate the average psf vector
            _avePsf.reset(new DVector(psfSize));
            Assert(psfSize == int(_avePsf->size()));
            _avePsf->setZero();
            nGoodPsf = 0;
            for(int n=0;n<nStars;++n) if ( flags[n]==0 ) {
                Assert(psf[n].getSigma() == _sigma);
                *_avePsf += psf[n].vec();
                ++nGoodPsf;
            }
            if (nGoodPsf == 0) {
                dbg<<"ngoodpsf = 0 in FittedPsf::calculate\n";
                throw ProcessingException("No good stars found for interpolation.");
            }
            *_avePsf /= double(nGoodPsf);

            // Rotate the vectors into their eigen directions.
            // The matrix V is stored to let us get back to the original basis.
            DMatrix mM(nGoodPsf,psfSize);
            DDiagMatrix inverseSigma(nGoodPsf);
            int i=0;
            for(int n=0;n<nStars;++n) if ( flags[n]==0 ) {
                Assert(int(psf[n].size()) == psfSize);
                Assert(i < nGoodPsf);
                mM.row(i) = psf[n].vec() - *_avePsf;
                inverseSigma(i) = nu[n];
                _bounds += pos[n];
                ++i;
            }
            Assert(i == nGoodPsf);
            xdbg<<"bounds = "<<_bounds<<std::endl;
            mM = inverseSigma EIGEN_asDiag() * mM;

            int nPcaTot = std::min(nGoodPsf,psfSize);
            DDiagMatrix mS(nPcaTot);
#ifdef USE_TMV
            DMatrixView mU = mM.colRange(0,nPcaTot);
            _mV.reset(new tmv::Matrix<double,tmv::RowMajor>(nPcaTot,psfSize));
            if (nGoodPsf > psfSize) {
                SV_Decompose(mU.view(),mS.view(),_mV->view(),true);
            } else {
                *_mV = mM;
                SV_Decompose(_mV->transpose(),mS.view(),mU.transpose());
            }
            xdbg<<"In FittedPSF: SVD S = "<<mS.diag()<<std::endl;
#else
            DMatrix mU(mM.TMV_colsize(),nPcaTot);
            _mV_transpose.reset(new DMatrix(psfSize,nPcaTot));
            if (nGoodPsf > psfSize) {
                Eigen::SVD<DMatrix> svd = TMV_colRange(mM,0,nPcaTot).svd();
                mU = svd.matrixU();
                mS = svd.singularValues();
                *_mV_transpose = svd.matrixV();
            } else {
                Eigen::SVD<Eigen::Transpose<DMatrix>::PlainMatrixType > svd = mM.transpose().svd();
                mU = svd.matrixV();
                mS = svd.singularValues();
                *_mV_transpose = svd.matrixU();
            }
            xdbg<<"In FittedPSF: SVD S = "<<EIGEN_Transpose(mS)<<std::endl;
#endif
            if (_params.keyExists("fitpsf_npca")) {
                _nPca = _params["fitpsf_npca"];
                dbg<<"npca = "<<_nPca<<" from parameter file\n";
            } else {
                double thresh = mS(0);
                if (_params.keyExists("fitpsf_pca_thresh")) 
                    thresh *= double(_params["fitpsf_pca_thresh"]);
                else thresh *= std::numeric_limits<double>::epsilon();
                dbg<<"thresh = "<<thresh<<std::endl;
                for(_nPca=1;_nPca<int(mM.TMV_rowsize());++_nPca) {
                    if (mS(_nPca) < thresh) break;
                }
                dbg<<"npca = "<<_nPca<<std::endl;
            }
#ifdef USE_TMV
            mU.colRange(0,_nPca) *= mS.subDiagMatrix(0,_nPca);
#else
            TMV_colRange(mU,0,_nPca) *= mS.TMV_subVector(0,_nPca).asDiagonal();
#endif
            xdbg<<"After U *= S\n";
            // U S = M(orig) * Vt

            while (nGoodPsf <= _fitSize && _fitSize > 1) {
                --_fitOrder;
                _fitSize = (_fitOrder+1)*(_fitOrder+2)/2;
                dbg<<"Too few good stars... reducing order of fit to "<<
                    _fitOrder<<std::endl;
            }
            DMatrix mP(nGoodPsf,_fitSize);
            mP.setZero();
            i=0;
            for(int n=0;n<nStars;++n) if ( flags[n]==0 ) {
                xdbg<<"n = "<<n<<" / "<<nStars<<std::endl;
#ifdef USE_TMV
                setPRow(_fitOrder,pos[n],_bounds,mP.row(i));
#else
                DVector mProwi(mP.TMV_rowsize());
                setPRow(_fitOrder,pos[n],_bounds,mProwi);
                mP.row(i) = mProwi.transpose();
#endif
                ++i;
            }
            Assert(i == nGoodPsf);
            mP = inverseSigma EIGEN_asDiag() * mP;
            xdbg<<"after mP = sigma * mP\n";

#ifdef USE_TMV
            _f.reset(new DMatrix(TMV_colRange(mU,0,_nPca)/mP));
#else
            _f.reset(new DMatrix(_fitSize,_nPca));
            mP.qr().solve(TMV_colRange(mU,0,_nPca),&(*_f));
#endif
            xdbg<<"Done making FittedPSF\n";

            //
            // Remove outliers from the fit using the empirical covariance matrix
            // of the data with respect to the fitted values.
            //
            xdbg<<"Checking for outliers:\n";

            // Calculate the covariance matrix
            DMatrix cov(psfSize,psfSize);
            cov.setZero();
            for(int n=0;n<nStars;++n) if ( flags[n]==0 ) {
                const BVec& data = psf[n];
                DVector fit(psfSize);
                interpolateVector(pos[n],TMV_vview(fit));
                DVector diff = data.vec() - fit;
#ifdef USE_TMV
                cov += diff ^ diff;
#else
                cov += diff * diff.transpose();
#endif
            }

            chisq = (mP * *_f - TMV_colRange(mU,0,_nPca)).TMV_normSq();
            dbg<<"chisq calculation #1 = "<<chisq<<std::endl;
            dof = nGoodPsf - _fitSize;

            if (dof > 0) { 
                cov /= double(dof);
            }
#ifdef USE_TMV
            cov.divideUsing(tmv::SV);
            cov.saveDiv();
            cov.setDiv();
            dbg<<"cov S = "<<cov.svd().getS().diag()<<std::endl;
#else
            Eigen::SVD<DMatrix> cov_svd = cov.svd();
            dbg<<"cov S = "<<EIGEN_Transpose(cov_svd.singularValues())<<std::endl;

#endif

            // Clip out 3 sigma outliers:
            nOutliers = 0;
            chisq = 0;
            for(int n=0;n<nStars;++n) if ( flags[n]==0 ) {
                const BVec& data = psf[n];
                DVector fit(psfSize);
                interpolateVector(pos[n],TMV_vview(fit));
                DVector diff = data.vec() - fit;
#ifdef USE_TMV
                double dev = diff * cov.inverse() * diff;
#else
                DVector temp(psfSize);
                cov_svd.solve(diff,&temp);
                double dev = (diff.transpose() * temp)(0,0);
#endif
                chisq += dev;
                if (dev > outlierThresh) {
                    xdbg<<"n = "<<n<<" is an outlier.\n";
                    xdbg<<"data = "<<data.vec()<<std::endl;
                    xdbg<<"fit = "<<fit<<std::endl;
                    xdbg<<"diff = "<<diff<<std::endl;
#ifdef USE_TMV
                    xdbg<<"diff/cov = "<<diff/cov<<std::endl;
#else
                    xdbg<<"diff/cov = "<<temp<<std::endl;
#endif
                    xdbg<<"dev = "<<dev<<std::endl;
                    ++nOutliers;
                    flags[n] |= PSF_INTERP_OUTLIER;
                }
            }
            dbg<<"ngoodpsf = "<<nGoodPsf<<std::endl;
            dbg<<"nOutliers = "<<nOutliers<<std::endl;
            dbg<<"chisq calculation #2 = "<<chisq<<std::endl;

        } while (nOutliers > 0);
    }

    FittedPsf::FittedPsf(const ConfigFile& params) : 
        _params(params), _psfOrder(_params.read<int>("psf_order")),
        _fitOrder(_params.read<int>("fitpsf_order")),
        _fitSize((_fitOrder+1)*(_fitOrder+2)/2)
    { }

    double FittedPsf::interpolateSingleElement(Position pos, int i) const
    {
        DVector P(_fitSize);
#ifdef USE_TMV
        setPRow(_fitOrder,pos,_bounds,P.view());
        DVector b1 = P * (*_f);
        double bi = b1 * _mV->col(i,0,_nPca);
#else
        setPRow(_fitOrder,pos,_bounds,P);
        DVector b1 = _f->transpose() * P;
        double bi = EIGEN_ToScalar(
            EIGEN_Transpose(b1) * _mV_transpose->TMV_rowpart(i,0,_nPca));
#endif
        bi += (*_avePsf)(i);
        return bi;
    }

    void FittedPsf::interpolateVector(Position pos, DVectorView b) const
    {
        DVector P(_fitSize);
#ifdef USE_TMV
        setPRow(_fitOrder,pos,_bounds,P.view());
        DVector b1 = P * (*_f);
        b = b1 * _mV->rowRange(0,_nPca);
#else
        setPRow(_fitOrder,pos,_bounds,P);
        DVector b1 = _f->transpose() * P;
        b = TMV_colRange(*_mV_transpose,0,_nPca) * b1;
#endif
        b += *_avePsf;
    }

}}}}
