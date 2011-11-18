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

#include <vector>

#include "lsst/meas/algorithms/shapelet/Ellipse.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include "lsst/meas/algorithms/shapelet/NLSolver.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class CrudeSolver : public NLSolver 
    {
    public :
        CrudeSolver(
            const PixelList& pix, double sigma, double I1,
            DVector& xinit);
        ~CrudeSolver() {}

        void calculateF(const DVector& x, DVector& f) const;
        void calculateJ(
            const DVector& x, const DVector& f,
            DMatrix& df) const;

    private : 
        double _sigma;
        DVector _I;
        CDVector _Z;
        DDiagMatrix _W;
        mutable CDVector _Z1;
        mutable DVector _E;
        mutable DVector _Rsq;
        mutable DVector _f1;
        double _I1;
        DVector& _xInit;
    };

    CrudeSolver::CrudeSolver(
        const PixelList& pix,
        double sigma, double I1, DVector& xInit) :
        _sigma(sigma), _I(pix.size()), _Z(pix.size()), _W(pix.size()), 
        _Z1(pix.size()), _E(pix.size()), _Rsq(pix.size()), _f1(pix.size()), 
        _I1(I1), _xInit(xInit)
    {
        const int nPix = pix.size();
        for(int i=0;i<nPix;++i) {
            _Z(i) = pix[i].getPos();
            _I(i) = pix[i].getFlux();
            double mask = exp(-std::norm(pix[i].getPos()/_sigma)/2.);
            _W(i) = pix[i].getInverseSigma()*mask;
        }
    }

    void CrudeSolver::calculateF(const DVector& x, DVector& f) const
    {
        Assert(x.size() == 4);
        Assert(f.size() == _Z.size());

        if ((x-_xInit).TMV_subVector(0,3).TMV_normInf() > 2.) {
            f = 2.e10*_f1;
            f.TMV_addToAll(1.);
            return; 
        }
        if (x(3) < 0.) { 
            f = 2.e10*_f1;
            f.TMV_addToAll(1.);
            return; 
        }

        std::complex<double> zc(x[0],x[1]);
        double mu = x[2];
        double I0 = _I1 * x[3];

        double m0 = exp(-mu)/_sigma;
        // z' = m*(z-zc)
        //    = m*z - m*zc
        // x' = m x - m xc
        // y' = m y - m yc
        _Z1 = m0*_Z;
        _Z1.TMV_addToAll(-m0*zc);

        const int nPix = _Z.size();
        for(int i=0;i<nPix;++i) {
            double rsq = std::norm(_Z1[i]);
            _Rsq[i] = rsq;
            _E[i] = exp(-rsq/2.);
        }
        _f1 = _I - I0*_E;
#ifdef USE_TMV
        _f1 = _W * _f1;
#else
        _f1.array() *= _W.array();
#endif
        f = _f1;
    }

    void CrudeSolver::calculateJ(
        const DVector& x, const DVector& f, DMatrix& df) const
    {
        //xdbg<<"Start J\n";
        Assert(x.size() == 4);
        Assert(f.size() == _Z.size());
        Assert(df.TMV_colsize() == _Z.size());
        Assert(df.TMV_rowsize() == 4);

        double mu = x[2];
        double I0 = _I1 * x[3];
        double m0 = exp(-mu)/_sigma;

        // fi = Wi * (Ii - I0 Ei)
        // dfi/dI0 = -Wi Ei
        // dfi/dmu = Wi I0 Ei (x1 dx1/dxc + y1 dy1/dxc) 
        // likewise for the other 4
        //
        // dx1/dxc = -m         dy1/dxc = 0
        // dx1/dyc = 0          dy1/dyc = -m
        // dx1/dmu = -x1        dy1/dmu = -y1
        //

#ifdef USE_TMV
        df.col(0) = -m0 * _Z1.realPart();
        df.col(1) = -m0 * _Z1.imagPart();
        df.col(2) = -_Rsq;
        df.colRange(0,3) = I0 * DiagMatrixViewOf(_E) * df.colRange(0,3);
        df.col(3) = -_I1 * _E;
        df = _W * df;
#else
        df.col(0).array() = _W.array() * _Z1.real().array();
        df.col(0).array() = _E.array() * df.col(0).array();
        df.col(0) *= -m0 * I0;
        df.col(1).array() = _W.array() * _Z1.imag().array();
        df.col(1).array() = _E.array() * df.col(1).array();
        df.col(1) *= -m0 * I0;
        df.col(2).array() = _W.array() * _Rsq.array();
        df.col(2).array() = _E.array() * df.col(2).array();
        df.col(2) *= -I0;
        df.col(3).array() = _W.array() * _E.array();
        df.col(3) *= -_I1;
#endif
    }

    void Ellipse::crudeMeasure(const PixelList& pix, double sigma)
    {
        // We use as our initial estimate the exact value for a 
        // well-sampled, uniform-variance, undistorted Gaussian intensity pattern
        // with no PSF convolution.
        //
        // That is, we assume the model:
        //
        // I(x,y) = I0 exp( -(|z-zc|^2/ (2 sigma'^2) )
        // where zz = (z-zc) 
        // and sigma' = exp(mu) sigma
        // 
        xdbg<<"Current centroid = "<<_cen<<std::endl;
        xdbg<<"Current mu = "<<_mu<<std::endl;
        xdbg<<"sigma = "<<sigma<<std::endl;
        const int nPix = pix.size();

#if 1
        // With a weight of exp(-|z|^2/(2 sigma^2)),
        // the weighted moments of this function are:
        // Iz/I = zc / (1+exp(2mu))
        // Irr/I = ( |zc|^2 + 2 exp(2mu) sigma^2 ) / (1+exp(2mu))

        std::complex<double> Iz = 0.;
        double I = 0.;
        double sig2 = sigma * exp(real(_mu));
        for(int i=0;i<nPix;++i) {
            double wt = exp(-std::norm((pix[i].getPos()-_cen)/sig2)/2.);
            Iz += wt * pix[i].getFlux() * (pix[i].getPos()-_cen);
            I += wt * pix[i].getFlux();
            if (std::abs(pix[i].getPos()-_cen) < 2.) 
                xdbg<<pix[i].getPos()<<"  "<<pix[i].getFlux()<<std::endl;
        }
        xdbg<<"Iz = "<<Iz<<", I = "<<I<<std::endl;

        // If I <= 0 then this isn't going to work.  Just return and hope
        // the regular measure method might do better.
        if (!(I > 0.)) return;

        std::complex<double> zc = Iz / I;
        // If zc is more than 2, something is probably wrong, so abort now.
        if (!(std::abs(zc) < 2.)) return;

        xdbg<<"Initial offset to centroid = "<<zc<<std::endl;
        if (isFixedCen()) {
            xdbg<<"But centroid is fixed, so don't apply.\n";
            zc = _cen;
        } else {
            zc += _cen;
        }
        xdbg<<"zc = "<<zc<<std::endl;

        double Irr = 0.;
        double W = 0.;
        I = 0.;
        Iz = 0.;
        for(int i=0;i<nPix;++i) {
            double wt = exp(-std::norm((pix[i].getPos()-zc)/sig2)/2.);
            Iz += wt * pix[i].getFlux() * (pix[i].getPos()-zc);
            Irr += wt * pix[i].getFlux() * std::norm(pix[i].getPos()-zc);
            I += wt * pix[i].getFlux();
            W += wt;
        }
        xdbg<<"Iz = "<<Iz<<", Irr = "<<Irr<<", I = "<<I<<", W = "<<W<<std::endl;

        std::complex<double> zc1 = Iz/I;
        double S = Irr/I - norm(zc1);
        // S is now 2 exp(2mu) sigma^2 / (1 + exp(2mu))
        xdbg<<"S = "<<S<<std::endl;
        double exp2mu = S / 2. / (sig2*sig2);
        xdbg<<"exp2mu/(1+exp2mu) = "<<exp2mu<<std::endl;
        // It's actually exp(2mu) / (1+exp(2mu)) at this point
        if (exp2mu < 0.2)
            exp2mu = 0.25;
        else if (exp2mu < 0.8)
            exp2mu = 1./(1./exp2mu-1.); // Now it is really exp(2mu)
        else 
            // The above formula is unstable, and probably inappropriate, since
            // we probably have a failure of our model approximation.
            // So just multiply it by 5 -- the correct factor for exp2mu = 0.8
            exp2mu *= 5.;
        xdbg<<"exp2mu = "<<exp2mu<<std::endl;
        if (exp2mu <= 0.) exp2mu = 1.;

        double m = log(exp2mu)/2.;
        xdbg<<"mu = "<<m<<std::endl;

        if (!isFixedCen()) zc += zc1 * (1.+exp2mu);
        if (isFixedMu()) m = real(_mu);
        else  m += real(_mu);

        xdbg<<"Approx cen = "<<zc<<std::endl;
        xdbg<<"Approx mu = "<<m<<std::endl;

        // I/W = I0 exp(2mu) / (1+exp(2mu)) * exp(-|zc1|^2/2sigma^2*(1+exp(2mu)))
        double I0 = (I/W)*(1.+exp2mu)/exp2mu /
            exp(-norm(zc1)*(1.+exp2mu)/(2.*sig2*sig2));
#else
        std::complex<double> zc = _cen;
        double m = real(_mu);

        double model = 0.;
        double obs = 0.;
        double minx = 1.e100, maxx = -1.e100, miny = 1.e100, maxy = -1.e100;
        for(int i=0;i<nPix;++i) {
            if (real(pix[i].getPos()) < minx) minx = real(pix[i].getPos());
            if (real(pix[i].getPos()) > maxx) maxx = real(pix[i].getPos());
            if (imag(pix[i].getPos()) < miny) miny = imag(pix[i].getPos());
            if (imag(pix[i].getPos()) > maxy) maxy = imag(pix[i].getPos());
            double wt = exp(-std::norm(pix[i].getPos()/sigma)/2.);
            model += mask;
            obs += pix[i].getFlux() * mask;
        }
        double pixarea = (maxx-minx)*(maxy-miny)/nPix;
        xdbg<<"pixarea = "<<pixarea<<std::endl; 
        xdbg<<"obs = "<<obs<<std::endl;
        xdbg<<"model = "<<model<<std::endl;
        double I0 = 2.*obs / model;
#endif

        xdbg<<"Initial I0 estimate = "<<I0<<std::endl;
        if (I0 < 1.e-6) {
            xdbg<<"Warning: small or negative I0: "<<I0<<" -- Use 1.0\n";
            I0 = 1.;
        }

        //if (std::abs(zc) > 1.0) zc /= std::abs(zc);
        //if (std::abs(m) > 1.0) m /= std::abs(m);
        DVector x(4);
        x[0] = std::real(zc); x[1] = std::imag(zc); 
        x[2] = m; 
        x[3] = 1.;
        DVector f(nPix);
        CrudeSolver s(pix,sigma,I0,x);

        s.useHybrid();
        s.setTol(1.e-4,1.e-8);
        s.setMinStep(1.e-15);
        s.setTau(1.0);
        s.setDelta0(0.05);
#ifdef __PGI
        s.noUseCholesky();
#endif
        if (XDEBUG) s.setOutput(*dbgout);
        xdbg<<"Before CrudeSolver: x = "<<EIGEN_Transpose(x)<<std::endl;
        s.solve(x,f);
        xdbg<<"After CrudeSolver: x = "<<EIGEN_Transpose(x)<<std::endl;
        s.setFTol(1.e-4 * (1.e-4 + std::abs(x[3])));
        s.solve(x,f);
        xdbg<<"After 2nd CrudeSolver: x = "<<EIGEN_Transpose(x)<<std::endl;

        std::complex<double> cenNew(x[0],x[1]);
        double muNew = x[2];

        if (std::abs(cenNew-_cen) > 2.) {
            dbg<<"Warning: large centroid shift in CrudeMeasure\n";
            dbg<<"Old centroid = "<<_cen<<", new centroid = "<<cenNew<<std::endl;
            cenNew = _cen + 2.*(cenNew - _cen)/std::abs(cenNew-_cen);
            dbg<<"Scaling back to "<<cenNew<<std::endl;
        }

        if (std::abs(muNew-real(_mu)) > 2.) {
            dbg<<"Warning: large scale change in CrudeMeasure\n";
            dbg<<"Old mu = "<<_mu<<", new mu = "<<muNew<<std::endl;
            muNew = real(_mu) + 2.*(muNew - real(_mu))/std::abs(muNew-real(_mu));
            dbg<<"Scaling back to "<<muNew<<std::endl;
        }

        xdbg<<"Crude cen = "<<cenNew<<std::endl;
        xdbg<<"Crude mu = "<<muNew<<std::endl;

        if (!_isFixedCen) _cen = cenNew;
        if (!_isFixedMu) _mu = muNew;
    }

    void Ellipse::crudeMeasure(
        const std::vector<PixelList>& pix, double sigma)
    {
        int nPix = 0;
        const int nPixList = pix.size();
        for(int i=0;i<nPixList;++i) nPix += pix[i].size();
        PixelList allPix(nPix);
        for(int i=0,k=0;i<nPixList;++i)
            for(size_t j=0;j<pix[i].size();++j,++k)
                allPix[k] = pix[i][j];
        crudeMeasure(allPix,sigma);
    }

    void Ellipse::peakCentroid(const PixelList& pix, double maxR)
    {
        double IPeak = 0.;
        std::complex<double> zPeak = 0.;
        const int nPix = pix.size();
        for(int i=0;i<nPix;++i) if (std::abs(pix[i].getPos()) < maxR) {
            if (pix[i].getFlux() > IPeak) { 
                zPeak = pix[i].getPos();
                IPeak = pix[i].getFlux();
            }
        }
        _cen = zPeak;
    }

}}}}
