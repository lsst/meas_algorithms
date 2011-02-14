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

#include <cmath>
#include <fstream>

#include "lsst/meas/algorithms/shapelet/Ellipse.h"
#include "lsst/meas/algorithms/shapelet/EllipseSolver.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"
#include "lsst/meas/algorithms/shapelet/Params.h"

#define N_FLUX_ATTEMPTS 0
#define MAXITER 4

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    static std::complex<double> addShears(
        const std::complex<double> g1, const std::complex<double> g2)
    {
        double absg1 = std::abs(g1);
        if (absg1 == 0.) return g2;
        double absd1 = tanh(atanh(absg1)*2.);
        std::complex<double> d1 = absd1 / absg1 * g1;
        double absg2 = std::abs(g2);
        if (absg2 == 0.) return g1;
        double absd2 = tanh(atanh(absg2)*2.);
        std::complex<double> d2 = absd2 / absg2 * g2;

        double d2sq = absd2*absd2;
        double x = (1.-std::sqrt(1.-d2sq))/d2sq;
        std::complex<double> d3 = d1+d2*(1.+x*(std::real(d1)*d2-std::real(d2)*d1));
        d3 /= 1. + std::real(d1*std::conj(d2));
        double absd3 = std::abs(d3);
        double absg3 = tanh(atanh(absd3)/2.);
        std::complex<double> g3 = absg3/absd3*d3;
        xdbg<<"GammaAdd: "<<g1<<" + "<<g2<<" = "<<g3<<std::endl;
        return g3;
    }

    bool Ellipse::doMeasure(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>* psf,
        int order, double sigma, bool shouldUseInteg, long& flag, 
        DMatrix* cov, BVec* bRet, DMatrix* bCov)
    {
#if 0
        timeval tp;
        double t1=0.,t2=0.;
#endif

        // The non-linear solver is pretty sensitive to having a good
        // initial estimate.  So start with a simple estimate.
        if (order > 3) {
            doMeasure(pix,psf,order-2,sigma,shouldUseInteg,flag);
        }

        xdbg<<"Start DoMeasure: order = "<<order<<", psf = "<<bool(psf)<<std::endl;
        xdbg<<"fix = "<<_isFixedCen<<"  "<<_isFixedGamma<<"  "<<_isFixedMu<<std::endl;
        xdbg<<"useInteg? "<<shouldUseInteg<<std::endl;
        for(size_t i=0;i<pix.size();++i) xdbg<<"npix["<<i<<"] = "<<pix[i].size()<<std::endl;

        std::auto_ptr<BaseEllipseSolver> solver;

        DVector x(5);
        x(0) = std::real(_cen);  x(1) = std::imag(_cen); 
        x(2) = std::real(_gamma); x(3) = std::imag(_gamma);
        x(4) = _mu;
        xdbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
        DVector xinit = x;
        DVector f(5);

        if (shouldUseInteg && order <= 3) {
#if 0
            if (_shouldDoTimings) {
                gettimeofday(&tp,0);
                t1 = tp.tv_sec + tp.tv_usec/1.e6;
            }
#endif
            // First we make the approximations the the galaxy is well-sampled
            // and has uniform variance.
            // Also, we start by just fitting the centroid.
            if (!_isFixedCen && (!_isFixedGamma || !_isFixedMu)) {
                xdbg<<"Do Integrating centroid solve\n";
                if (psf) {
                    // pixscale doesn't really matter unless we want an accurate B in
                    // the end, so just use 1 here.
                    solver.reset(new EllipseSolver2(
                            pix,*psf,_fPsf,order,sigma,1.,false,true,true));
                } else {
                    solver.reset(new EllipseSolver2(
                            pix,order,sigma,1.,false,true,true));
                }
                xdbg<<"xinit (int fixgam) = "<<EIGEN_Transpose(x)<<std::endl;
                solver->useDogleg();
#ifdef NOTHROW
                solver->noUseCholesky();
#endif
                solver->setTol(1.e-3,1.e-8);
                if (XDEBUG) solver->setOutput(*dbgout);
                solver->setDelta0(0.01);
                solver->setMinStep(1.e-12);
                solver->setMaxIter(60);
                xdbg<<"Integrating solver, centroid only:\n";
                //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
                bool ret = solver->solve(x,f);
                if (!ret) {
                    dbg<<"failed integrating solver, centroid - x = "<<EIGEN_Transpose(x)<<std::endl;
                    dbg<<"f = "<<EIGEN_Transpose(f)<<"  Norm(f) = "<<f.norm()<<std::endl;
                    dbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
                    return false;
                }
                xdbg<<"Done: x (Integ, centroid only) = "<<EIGEN_Transpose(x)<<std::endl;
                xdbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
                xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
            }

            // Next allow the shear and/or mu to be fit as well:
#ifdef N_FLUX_ATTEMPTS
#if N_FLUX_ATTEMPTS > 0
            xdbg<<"Do Integrating regular solve (fixed flux)\n";
            if (psf) {
                // pixscale doesn't really matter unless we want an accurate B in
                // the end, so just use 1 here.
                solver.reset(
                    new EllipseSolver2(
                        pix,*psf,_fPsf,order,sigma,1.,
                        _isFixedCen,_isFixedGamma,_isFixedMu,true));
            } else {
                solver.reset(new EllipseSolver2(
                        pix,order,sigma,1.,
                        _isFixedCen,_isFixedGamma,_isFixedMu,true));
            }
            xdbg<<"xinit (integ) = "<<EIGEN_Transpose(x)<<std::endl;
            solver->useHybrid();
#ifdef NOTHROW
            solver->noUseCholesky();
            solver->noUseDirectH();
#endif
            solver->setTol(3.e-3,1.e-5);
            solver->setTau(1.0);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setMinStep(1.e-12);
            solver->setMaxIter(10);
            xdbg<<"Integrating solver:\n";
            for(int iter = 0; iter < N_FLUX_ATTEMPTS; ++iter) {
                xdbg<<"Attempt #"<<iter<<std::endl;
                //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
                solver->solve(x,f);
                xdbg<<"x => "<<EIGEN_Transpose(x)<<std::endl;
                xdbg<<"f => "<<EIGEN_Transpose(f)<<"  Norm(f) = "<<f.norm()<<std::endl;
                xdbg<<"b => "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
                if (f.norm() < 1.e-2) break;
            } 
#endif
#endif
            // Repeat, but allow flux to change.
            xdbg<<"Do Integrating regular solve with flux\n";
            if (psf) {
                solver.reset(
                    new EllipseSolver2(
                        pix,*psf,_fPsf,order,sigma,1.,
                        _isFixedCen,_isFixedGamma,_isFixedMu));
            } else {
                solver.reset(new EllipseSolver2(
                        pix,order,sigma,1.,
                        _isFixedCen,_isFixedGamma,_isFixedMu));
            }
            xdbg<<"xinit (integ) = "<<EIGEN_Transpose(x)<<std::endl;
            solver->useDogleg();
#ifdef NOTHROW
            solver->noUseCholesky();
#endif
            solver->setTol(1.e-3,1.e-8);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setDelta0(0.01);
            solver->setMinStep(1.e-12);
            solver->setMaxIter(10);
            xdbg<<"Final Integrating solver:\n";
            //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
            solver->solve(x,f);
            xdbg<<"Done: x (Integrating) = "<<EIGEN_Transpose(x)<<std::endl;
            xdbg<<"f (integ) = "<<EIGEN_Transpose(f)<<std::endl;
            xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
#if 0
            if (_shouldDoTimings) {
                gettimeofday(&tp,0);
                t2 = tp.tv_sec + tp.tv_usec/1.e6;
                _times._tInteg += t2-t1;
            }
#endif
        }

#if 0
        // We should be close enough for the exact solver to work.
        // But just to be safe, fit each of centroid, gamma, mu separately first:
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }
#endif
        if (!_isFixedCen) {
            if (psf) {
                solver.reset(new EllipseSolver(
                        pix,*psf,_fPsf,order,sigma,false,true,true));
            } else {
                solver.reset(new EllipseSolver(
                        pix,order,sigma,false,true,true));
            }
            xdbg<<"xinit = "<<EIGEN_Transpose(x)<<std::endl;
            solver->useDogleg();
#ifdef NOTHROW
            solver->noUseCholesky();
#endif
            solver->setTol(1.e-3,1.e-8);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setDelta0(0.01);
            solver->setMinStep(1.e-12);
            solver->setMaxIter(50);
            xdbg<<"ML solver centroid:\n";
            bool rete= solver->solve(x,f);
            if (!ret) {
                dbg<<"ML solver, centroid, failed - x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
                dbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
                if (XDEBUG) solver->testJ(x,f,dbgout,1.e-5);
#if 0
                if (_shouldDoTimings) {
                    gettimeofday(&tp,0);
                    t2 = tp.tv_sec + tp.tv_usec/1.e6;
                    _times._tCentroid += t2-t1;
                }
#endif
                return false;
            }
            xdbg<<"Done: x (ML centroid only) = "<<EIGEN_Transpose(x)<<std::endl;
            xdbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
        }
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            _times._tCentroid += t2-t1;
        }
#endif

#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }
#endif
        if (!_isFixedGamma) {
            // start with an initial guess from the latest b vector:
            if (solver.get()) {
                BVec b1 = solver->getB();
                std::complex<double> g1(b1(3)/b1(0)*sqrt(2.),-b1(4)/b1(0)*sqrt(2.));
                dbg<<"g1 = "<<g1<<std::endl;
                if (std::abs(g1) < 1.) {
                    std::complex<double> g0(x[2],x[3]);
                    std::complex<double> g = addShears(g0,g1);
                    dbg<<"Initial gamma = "<<g<<std::endl;
                    x[2] = std::real(g);
                    x[3] = std::imag(g);
                }
            }

            if (psf) {
                solver.reset(new EllipseSolver(
                        pix,*psf,_fPsf,order,sigma,true,false,true));
            } else {
                solver.reset(new EllipseSolver(
                        pix,order,sigma,true,false,true));
            }
            xdbg<<"xinit = "<<EIGEN_Transpose(x)<<std::endl;
            solver->useDogleg();
#ifdef NOTHROW
            solver->noUseCholesky();
#endif
            solver->setTol(1.e-3,1.e-8);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setDelta0(0.01);
            solver->setMinStep(1.e-12);
            solver->setMaxIter(50);
            xdbg<<"ML solver gamma:\n";
            bool ret = solver->solve(x,f);
            if (!ret) {
                dbg<<"ML solver, gamma, failed - x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
                dbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
                //if (XDEBUG) solver->testJ(x,f,dbgout,1.e-5);
#if 0
                if (_shouldDoTimings) {
                    gettimeofday(&tp,0);
                    t2 = tp.tv_sec + tp.tv_usec/1.e6;
                    _times._tGamma += t2-t1;
                }
#endif
                return false;
            }
            xdbg<<"Done: x (ML gamma only) = "<<EIGEN_Transpose(x)<<std::endl;
            xdbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
        }
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            _times._tGamma += t2-t1;
        }
#endif

#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }
#endif
        if (!_isFixedMu) {
            if (psf) {
                solver.reset(new EllipseSolver(
                        pix,*psf,_fPsf,order,sigma,true,true,false));
            } else {
                solver.reset(new EllipseSolver(
                        pix,order,sigma,true,true,false));
            }
            xdbg<<"xinit = "<<EIGEN_Transpose(x)<<std::endl;
            solver->useDogleg();
#ifdef NOTHROW
            solver->noUseCholesky();
#endif
            solver->setTol(1.e-3,1.e-8);
            if (XDEBUG) solver->setOutput(*dbgout);
            solver->setDelta0(0.01);
            solver->setMinStep(1.e-12);
            solver->setMaxIter(50);
            xdbg<<"ML solver mu:\n";
            bool ret = solver->solve(x,f);
            if (!ret) {
                dbg<<"ML solver, mu, failed - x = "<<EIGEN_Transpose(x)<<std::endl;
                dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
                dbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
                //if (XDEBUG) solver->testJ(x,f,dbgout,1.e-5);
#if 0
                if (_shouldDoTimings) {
                    gettimeofday(&tp,0);
                    t2 = tp.tv_sec + tp.tv_usec/1.e6;
                    _times._tGamma += t2-t1;
                }
#endif
                return false;
            }
            xdbg<<"Done: x (ML mu only) = "<<EIGEN_Transpose(x)<<std::endl;
            xdbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
        }
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            _times._tGamma += t2-t1;
        }
#endif
#endif

        // Now fit everything, but first time, try to maintain the flux level
#ifdef N_FLUX_ATTEMPTS
#if N_FLUX_ATTEMPTS > 0
        xdbg<<"Do regular solve (fixed flux)\n";
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }
        if (psf) {
            solver.reset(new EllipseSolver(
                    pix,*psf,_fPsf,order,sigma,
                    _isFixedCen,_isFixedGamma,_isFixedMu,true));
        } else {
            solver.reset(new EllipseSolver(
                    pix,order,sigma,
                    _isFixedCen,_isFixedGamma,_isFixedMu,true));
        }
        xdbg<<"xinit = "<<EIGEN_Transpose(x)<<std::endl;
        solver->useHybrid();
#ifdef NOTHROW
        solver->noUseCholesky();
        solver->noUseDirectH();
#endif
        solver->setTol(3.e-3,1.e-5);
        solver->setTau(1.0);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->setMinStep(1.e-12);
        solver->setMaxIter(10);
        solver->setDelta0(0.01);
        xdbg<<"ML solver use flux:\n";
        for(int iter = 0; iter < N_FLUX_ATTEMPTS; ++iter) {
            xdbg<<"Attempt #"<<iter<<std::endl;
            //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
            solver->solve(x,f);
            xdbg<<"x => "<<EIGEN_Transpose(x)<<std::endl;
            xdbg<<"f => "<<EIGEN_Transpose(f)<<"  Norm(f) = "<<f.norm()<<std::endl;
            xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
            if (f.norm() < 1.e-2) break;
        }
        xdbg<<"Done: x (ML fixed flux) = "<<EIGEN_Transpose(x)<<std::endl;
        xdbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
        xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            _times._tFixFlux += t2-t1;
        }
#endif
#endif
#endif

        xdbg<<"Do regular solve with flux\n";
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t1 = tp.tv_sec + tp.tv_usec/1.e6;
        }
#endif
        // Finally allow the flux change as needed.
        if (psf) {
            solver.reset(new EllipseSolver(
                    pix,*psf,_fPsf,order,sigma,
                    _isFixedCen,_isFixedGamma,_isFixedMu));
        } else {
            solver.reset(new EllipseSolver(
                    pix,order,sigma,
                    _isFixedCen,_isFixedGamma,_isFixedMu));
        }
        xdbg<<"xinit = "<<EIGEN_Transpose(x)<<std::endl;
        solver->useDogleg();
        if (psf) solver->useNumericJ();
#ifdef NOTHROW
        solver->noUseCholesky();
#endif
        solver->setTol(1.e-8,1.e-25);
        if (XDEBUG) solver->setOutput(*dbgout);
        solver->useVerboseOutput();
        solver->setDelta0(0.01);
        solver->setMinStep(1.e-15);
        solver->setMaxIter(50);
        xdbg<<"Final solver:\n";
        for(int iter = 1;iter<=MAXITER;++iter) {
            //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-6)) exit(1);
            bool ret = solver->solve(x,f);
            if (ret) break;
            else if (iter == MAXITER) {
#if 0
                if (_shouldDoTimings) {
                    gettimeofday(&tp,0);
                    t2 = tp.tv_sec + tp.tv_usec/1.e6;
                    _times._tFinal += t2-t1;
                }
#endif
                return false;
            }
            dbg<<"Final solver pass failed - x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
            //if (XDEBUG) if (!solver->testJ(x,f,dbgout,1.e-5)) exit(1);
#if 1 
            // Try bumping it to get out of a local well and continue:
            BVec b1 = solver->getB();
            if (!_isFixedCen) {
                xdbg<<"Bump z\n";
                dbg<<"Current z = "<<std::complex<double>(x[0],x[1])<<std::endl;
                std::complex<double> z1(2.*b1(1),-2.*b1(2));
                z1 /= b1(0) - b1(5);
                if (std::abs(z1) < 1.) {
                    dbg<<"Add "<<z1<<std::endl;
                    x[0] += std::real(z1); x[1] += std::imag(z1); 
                    dbg<<"New z = "<<std::complex<double>(x[0],x[1])<<std::endl;
                }
            }
#endif
#if 1
            if (!_isFixedGamma) {
                xdbg<<"Bump gamma\n";
                std::complex<double> g0(x[2],x[3]);
                dbg<<"Current gamma = "<<g0<<std::endl;
                std::complex<double> g1(std::sqrt(2.)*b1(3),-std::sqrt(2.)*b1(4));
                g1 /= b1(0) - (b1.getOrder() >= 4 ? b1(14) : 0.);
                if (std::abs(g1) < 1.) {
                    dbg<<"Add "<<g1<<std::endl;
                    g0 = addShears(g1,g0);
                    x[2] = std::real(g0); x[3] = std::imag(g0); 
                    dbg<<"New gamma = "<<std::complex<double>(x[2],x[3])<<std::endl;
                }
            }
#endif
#if 0
            if (!_isFixedMu) {
                xdbg<<"Bump mu\n";
                dbg<<"Current mu = "<<x[4]<<std::endl;
                double m1 = -b1(5);
                m1 /= b1[0] - (b1.getOrder() >= 4 ? 2.*b1(14) : 0.);
                if (std::abs(m1) < 1.) {
                    x[4] += m1; 
                    dbg<<"New mu = "<<x[4]<<std::endl;
                }
            }
            dbg<<"New x = "<<EIGEN_Transpose(x)<<std::endl;
#endif

            if (psf) {
                solver.reset(new EllipseSolver(
                        pix,*psf,_fPsf,order,sigma,
                        _isFixedCen,_isFixedGamma,true));
            } else {
                solver.reset(new EllipseSolver(
                        pix,order,sigma,
                        _isFixedCen,_isFixedGamma,true));
            }

            solver->useHybrid();
#ifdef NOTHROW
            solver->noUseCholesky();
            solver->noUseDirectH();
#endif
            if (iter == 1) {
                solver->setTol(1.e-8,1.e-10);
                if (XDEBUG) solver->setOutput(*dbgout);
                solver->setDelta0(0.01);
                solver->setMinStep(1.e-15);
                solver->setMaxIter(50);
            } else {
                solver->setTol(10.*solver->getFTol(),10.*solver->getGTol());
            }
        }
        xdbg<<"Done: Final solver pass successful: x = "<<EIGEN_Transpose(x)<<std::endl;
        xdbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
        xdbg<<"b = "<<EIGEN_Transpose(solver->getB().vec())<<std::endl;
        if (cov) {
            solver->useSVD();
            solver->getCovariance(*cov);
            xdbg<<"cov = "<<*cov<<std::endl;
        }

        // Check for flags
        xdbg<<"Checking for possible problems\n";
        solver->callF(x,f);
        if (f.norm() > solver->getFTol()) {
            xdbg<<"Local minimum\n";
            flag |= SHEAR_LOCAL_MIN;
        }
        if (solver->getFTol() > 1.e-8) {
            xdbg<<"ftol was raised\n";
            flag |= SHEAR_POOR_FIT;
        }

        if (XDEBUG && psf) {
            for(size_t k=0;k<pix.size();++k) {
                double sig_psf = (*psf)[k].getSigma();
                double D = sigma*sigma / (sig_psf*sig_psf);
                xdbg<<"D = "<<D<<std::endl;
            }
        }
#if 0
        if (_shouldDoTimings) {
            gettimeofday(&tp,0);
            t2 = tp.tv_sec + tp.tv_usec/1.e6;
            _times._tFinal += t2-t1;
        }
#endif

        setCen(std::complex<double>(x(0),x(1)));
        setGamma(std::complex<double>(x(2),x(3))); 
        setMu(x(4));
        if (bRet) {
            *bRet = solver->getB();
            xdbg<<"bret = "<<EIGEN_Transpose(bRet->vec())<<std::endl;
            if (bCov) solver->getBCov(*bCov);
        }
        return true;
    }

}}}}
