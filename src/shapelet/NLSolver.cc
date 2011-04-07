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
// The algorithms contained in this file are taken from the paper
// "Methods for Nonlinear Least-Squares Problems", by Madsen, Nielsen,
// and Tingleff (2004).  
// A copy of this paper should be included with the code in the file
// madsen04.pdf.  Please refer to this paper for more details about
// how the algorithms work.

#include <iostream>
#include <limits>
#include <algorithm>

#include "lsst/meas/algorithms/shapelet/NLSolver.h"

#define dbg if(_nlOut) (*_nlOut)
#define xdbg if(_verbose >= 1 && _nlOut) (*_nlOut)
#define xxdbg if(_verbose >= 2 && _nlOut) (*_nlOut)

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

#ifdef USE_TMV

    NLSolver::NLSolver() : 
        _method(NEWTON),
        _fTol(1.e-8), _gTol(1.e-8), _minStep(1.e-8), _maxIter(200),
        _tau(1.e-3), _delta0(1.), 
        _nlOut(0), _verbose(0),
        _hasDirectH(false), _shouldUseCh(true), _shouldUseSvd(false) 
    {}

    void NLSolver::calculateJ(
        const tmv::Vector<double>& x, const tmv::Vector<double>& f, 
        tmv::Matrix<double>& df) const
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
        // Do a finite difference calculation for J.
        // This function is virtual, so if there is a better way to 
        // calculate J, then you should override this version.

        tmv::Vector<double> x2 = x;
        tmv::Vector<double> f2(f.size());
        tmv::Vector<double> f1(f.size());
        for(size_t j=0;j<x.size();++j) {
            const double dx = sqrtEps * (x.norm() + 1.);
            x2(j) += dx;
            this->calculateF(x2,f2);
            x2(j) -= 2.*dx;
            this->calculateF(x2,f1);
            df.col(j) = (f2-f1)/(2.*dx);
            x2(j) = x(j);
        }
    }

    bool NLSolver::testJ(
        const tmv::Vector<double>& x, tmv::Vector<double>& f,
        std::ostream* os, double relErr) const
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());

        this->calculateF(x,f);
        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        this->calculateJ(x,f,J);
        tmv::Matrix<double> Jn(f.size(),x.size());
        NLSolver::calculateJ(x,f,Jn);
        double err = MaxAbsElement(J-Jn) / Jn.norm();
        if (!relErr) relErr = 10.*sqrtEps;
        if (os) {
            *os << "TestJ:\n";
            if (_verbose >= 1) {
                *os << "x = "<<x<<std::endl;
                *os << "f = "<<f<<std::endl;
                *os << "Direct J = "<<J<<std::endl;
                *os << "Numeric J = "<<Jn<<std::endl;
            }
            *os << "MaxAbsElement(J-J_num) / J.norm() = "<<err<<std::endl;
            *os << "cf. relerr = "<<relErr<<std::endl;
            if (err >= relErr) {
                tmv::Matrix<double> diff = J-Jn;
                *os << "J-J_num = "<<diff;
                double maxel = diff.maxAbsElement();
                *os << "Max element = "<<maxel<<std::endl;
                for(size_t i=0;i<diff.colsize();++i) {
                    for(size_t j=0;j<diff.rowsize();++j) {
                        if (std::abs(diff(i,j)) > 0.9*maxel) {
                            *os<<"J("<<i<<','<<j<<") = "<<J(i,j)<<"  ";
                            *os<<"J_num("<<i<<','<<j<<") = "<<Jn(i,j)<<"  ";
                            *os<<"diff = "<<J(i,j)-Jn(i,j)<<std::endl;
                        }
                    }
                }
            }
        }
        return err < relErr;
    }

    class NoDefinedH : public std::runtime_error
    {
    public :
        NoDefinedH() :
            std::runtime_error("calculateH is undefined in NLSolver") 
        {}
    };

    void NLSolver::calculateH(
        const tmv::Vector<double>& , const tmv::Vector<double>& ,
        const tmv::Matrix<double>& , tmv::SymMatrix<double>& ) const
    { 
#ifdef NOTHROW
        std::cerr<<"H is undefined\n";
        exit(1);
#else
        throw NoDefinedH();
#endif
    }


    // H(i,j) = d^2 Q / dx_i dx_j
    // where Q = 1/2 Sum_k |f_k|^2
    // H = JT J + Sum_k f_k d^2f_k/(dx_i dx_j)
    void NLSolver::calculateNumericH(
        const tmv::Vector<double>& x,
        const tmv::Vector<double>& f, 
        tmv::SymMatrix<double>& h) const
    {
        // Do a finite difference calculation for H.

        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
        const double dx = sqrt(sqrtEps) * (x.norm() + 1.);
        double q0 = 0.5 * f.normSq();

        tmv::Vector<double> x2 = x;
        tmv::Vector<double> f2(f.size());
        for(size_t i=0;i<x.size();++i) {
            x2(i) = x(i) + dx;
            this->calculateF(x2,f2);
            double q2a = 0.5*f2.normSq();
            x2(i) = x(i) - dx;
            this->calculateF(x2,f2);
            double q2b = 0.5*f2.normSq();

            h(i,i) = (q2a + q2b - 2.*q0) / (dx*dx);
            x2(i) = x(i);

            for(size_t j=i+1;j<x.size();++j) {
                x2(i) = x(i) + dx;
                x2(j) = x(j) + dx;
                this->calculateF(x2,f2);
                q2a = 0.5*f2.normSq();

                x2(i) = x(i) + dx;
                x2(j) = x(j) - dx;
                this->calculateF(x2,f2);
                q2b = 0.5*f2.normSq();

                x2(i) = x(i) - dx;
                x2(j) = x(j) + dx;
                this->calculateF(x2,f2);
                double q2c = 0.5*f2.normSq();

                x2(i) = x(i) - dx;
                x2(j) = x(j) - dx;
                this->calculateF(x2,f2);
                double q2d = 0.5*f2.normSq();

                h(i,j) = (q2a - q2b - q2c + q2d) / (4.*dx*dx);
                x2(i) = x(i);
                x2(j) = x(j);
            }
        }
    }


#define CHECKF(normInfF) \
    do { \
        double checkfTemp = (normInfF); \
        if (!(checkfTemp > _fTol)) { \
            dbg<<"Found ||f|| ~= 0\n"; \
            dbg<<"||f||_inf = "<<checkfTemp<<" < "<<_fTol<<std::endl; \
            return true; \
        } \
    } while (false)

#define CHECKG(normInfG) \
    do { \
        double checkgTemp = (normInfG); \
        if (!(checkgTemp > _gTol)) { \
            dbg<<"Found local minimum of ||f||\n"; \
            dbg<<"||g||_inf = "<<checkgTemp<<" < "<<_gTol<<std::endl; \
            return true; \
        } \
    } while (false)

#define SHOWFAILFG \
    do { \
        dbg<<"||f||_inf = "<<f.normInf()<<" !< "<<_fTol<<std::endl; \
        dbg<<"||g||_inf = "<<g.normInf()<<" !< "<<_gTol<<std::endl; \
    } while (false)

#define CHECKSTEP(normH) \
    do { \
        double checkStepTemp1 = (normH); \
        double checkStepTemp2 = _minStep*(x.norm()+_minStep); \
        if (!(checkStepTemp1 > checkStepTemp2)) { \
            dbg<<"Step size became too small\n"; \
            dbg<<"||h|| = "<<checkStepTemp1<<" < "<<checkStepTemp2<<std::endl; \
            SHOWFAILFG; \
            return false; \
        } \
    } while (false)

    bool NLSolver::solveNewton(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // This is a simple descent method which uses either the 
        // Newton direction or the steepest descent direction.
    {
        const double gamma1 = 0.1;
        const double gamma2 = 0.5;
        dbg<<"Start Solve_Newton\n";

        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        tmv::Vector<double> g(x.size());
        tmv::Vector<double> h(x.size());
        tmv::Vector<double> xNew(x.size());
        tmv::Vector<double> fNew(f.size());
        tmv::Vector<double> gNew(x.size());

        xdbg<<"x = "<<x<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f<<std::endl;
        CHECKF(f.normInf());
        double Q = 0.5*f.normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        if (_shouldUseSvd) J.divideUsing(tmv::SV);
        J.saveDiv();
        xdbg<<"J = "<<J<<std::endl;
        g = J.transpose() * f;
        xdbg<<"g = "<<g<<std::endl;
        CHECKG(g.normInf());
        double alpha = Q/g.normSq();
        bool shouldUseNewton = true;

        dbg<<"iter   |f|inf   Q   |g|inf   alpha\n";
        for(int k=0;k<_maxIter;++k) {
            shouldUseNewton = true;
            dbg<<k<<"   "<<f.normInf()<<"   "<<Q<<"   "<<g.normInf()<<"   "<<
                alpha<<std::endl;

            h = -f/J;
            xdbg<<"h = "<<h<<std::endl;
            double normH = h.norm();
            CHECKSTEP(normH);

            // phi(alpha) = Q(x + alpha h)
            // phi'(alpha) = fT J h = gT h
            // where g is measured at xNew, not x
            // m = phi'(0)
            double m = h*g;
            double normG = g.norm();
            double m2 = -normG*normG;

            if ((k%5 == 0 && m >= 0.) || (k%5 != 0 && m/normH >= -0.01*normG)) {
                // i.e. either m >= 0 or |m/normH| < 0.01 * |m2/normH2|
                shouldUseNewton = false;
                xdbg<<"Newton is not a good descent direction - use steepest descent\n";
                h = -g;
                CHECKSTEP(normG);
                m = m2;
            } else {
                xdbg<<"m = "<<m<<", Steepest m = "<<m2<<std::endl;
                xdbg<<"m/h.norm() = "<<m/normH<<", Steepest m/h.norm() = "<<-normG<<std::endl;
            }

            if (shouldUseNewton && alpha > 0.1) alpha = 1.0;
            for(int k2=0;k2<=_maxIter;++k2) {
                if (k2 == _maxIter) { 
                    dbg<<"Maximum iterations exceeded in subloop of Newton method\n";
                    dbg<<"This can happen when there is a singularity (or close to it)\n";
                    dbg<<"along the gradient direction:\n";
                    if (_shouldUseSvd) {
                        dbg<<"J Singular values = \n"<<J.svd().getS().diag()<<std::endl;
                        dbg<<"V = \n"<<J.svd().getV()<<std::endl;
                    }
                    SHOWFAILFG; 
                    return false;
                }
                xNew = x + alpha*h;
                if (alpha < _minStep) {
                    dbg<<"alpha became too small ("<<alpha<<" < "<<_minStep<<")\n";
                    SHOWFAILFG; 
                    return false;
                }
                xdbg<<"xnew = "<<xNew<<std::endl;
                this->calculateF(xNew,fNew);
                xdbg<<"fnew = "<<fNew<<std::endl;
                double QNew = 0.5*fNew.normSq();
                xdbg<<"Qnew = "<<QNew<<std::endl;

                // Check that phi has decreased significantly
                // Require phi(alpha) <= phi(0) + gamma1 phi'(0) alpha
                if (QNew > Q + gamma1 * m * alpha) {
                    alpha /= 2;
                    shouldUseNewton = false;
                    xdbg<<"Qnew not small enough: alpha => "<<alpha<<std::endl;
                    continue;
                }
                this->calculateJ(xNew,fNew,J);
                J.unsetDiv();
                xdbg<<"Jnew = "<<J<<std::endl;
                gNew = J.transpose() * fNew;
                xdbg<<"gnew = "<<gNew<<std::endl;

                // Check that alpha is not too small
                // Require phi'(alpha) >= gamma2 phi'(0)
                double mNew = h*gNew;
                if (mNew < gamma2 * m) {
                    alpha *= 3.;
                    shouldUseNewton = false;
                    xdbg<<"New slope not shallow enough: alpha => "<<alpha<<std::endl;
                    xdbg<<"(m = "<<m<<", mnew = "<<mNew<<")\n";
                    continue;
                }
                xdbg<<"Good choice\n";
                x = xNew; f = fNew; Q = QNew; g = gNew;
                break;
            }
            CHECKF(f.normInf());
            CHECKG(g.normInf());
        }
        dbg<<"Maximum iterations exceeded in Newton method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solveLM(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // This is the Levenberg-Marquardt method
    {
        dbg<<"Start Solve_LM\n";

        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        tmv::Vector<double> h(x.size());
        tmv::Vector<double> xNew(x.size());
        tmv::Vector<double> fNew(f.size());
        tmv::Vector<double> gNew(x.size());

        xdbg<<"x = "<<x<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f<<std::endl;
        CHECKF(f.normInf());
        double Q = 0.5*f.normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        if (_shouldUseSvd) J.divideUsing(tmv::SV);
        J.saveDiv();
        xdbg<<"J = "<<J<<std::endl;
        tmv::Vector<double> g = J.transpose() * f;
        xdbg<<"g = "<<g<<std::endl;
        CHECKG(g.normInf());

        tmv::SymMatrix<double> A = J.transpose() * J;
        xdbg<<"JT J = "<<A<<std::endl;
        if (_shouldUseSvd) A.divideUsing(tmv::SV);
        else if (_shouldUseCh) A.divideUsing(tmv::CH);
        else A.divideUsing(tmv::LU);
        double mu = _tau * A.diag().normInf();
        xdbg<<"initial mu = "<<_tau<<" * "<<A.diag().normInf()<<" = "<<mu<<std::endl;
        A += mu;
        A.saveDiv();
        double nu = 2.;

        dbg<<"iter   |f|inf   Q   |g|inf   mu\n";
        for(int k=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<f.normInf()<<"   "<<Q<<"   "<<g.normInf()<<"   "<<mu<<std::endl;
            xdbg<<"k = "<<k<<std::endl;
            xdbg<<"mu = "<<mu<<std::endl;
            xdbg<<"A = "<<A<<std::endl;
#ifndef NOTHROW
            try {
#endif
                h = -g/A;
#ifndef NOTHROW
            } catch (tmv::NonPosDef) {
                xdbg<<"NonPosDef caught - switching division to LU method.\n";
                // Once the Cholesky decomp fails, just use LU from that point on.
                A.divideUsing(tmv::LU);
                h = -g/A;
            }
#endif
            xdbg<<"h = "<<h<<std::endl;
            CHECKSTEP(h.norm());

            xNew = x + h;
            xdbg<<"xnew = "<<xNew<<std::endl;
            this->calculateF(xNew,fNew);
            xdbg<<"fnew = "<<fNew<<std::endl;
            double QNew = 0.5*fNew.normSq();
            xdbg<<"Qnew = "<<QNew<<std::endl;

            if (QNew < Q) {
                xdbg<<"improved\n";
                x = xNew; f = fNew; 
                CHECKF(f.normInf());

                this->calculateJ(x,f,J);
                J.unsetDiv();
                A = J.transpose() * J;
                gNew = J.transpose() * f;
                xdbg<<"gnew = "<<gNew<<std::endl;
                CHECKG(gNew.normInf());

                // Use g as a temporary for (g - mu*h)
                g -= mu*h;
                double rho = (Q-QNew) / (-0.5*h*g);
                xdbg<<"rho = "<<Q-QNew<<" / "<<(-0.5*h*g)<<" = "<<rho<<std::endl;
                mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
                xdbg<<"mu *= "<<std::max(1./3.,1.-std::pow(2.*rho-1.,3))<<" = "<<mu<<std::endl;
                A += mu;
                A.unsetDiv();
                Q = QNew; g = gNew;
            } else {
                xdbg<<"not improved\n";
                A += mu*(nu-1.); mu *= nu; nu *= 2.;
                A.unsetDiv();
                xdbg<<"mu *= (nu = "<<nu<<") = "<<mu<<std::endl;
            }
        }
        dbg<<"Maximum iterations exceeded in LM method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solveDogleg(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // This is the Dogleg method
    {
        dbg<<"Start Solve_Dogleg\n";
        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        tmv::Vector<double> h(x.size());
        tmv::Vector<double> temp(x.size());
        tmv::Vector<double> xNew(x.size());
        tmv::Vector<double> fNew(f.size());

        xdbg<<"x = "<<x<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f<<std::endl;
        CHECKF(f.normInf());
        double Q = 0.5*f.normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        if (_shouldUseSvd) J.divideUsing(tmv::SV);
        J.saveDiv();
        xdbg<<"J = "<<J<<std::endl;
        xdbg<<"J.svd = "<<J.svd().getS().diag()<<std::endl;

        tmv::Vector<double> g = J.transpose() * f;
        xdbg<<"g = "<<g<<std::endl;
        CHECKG(g.normInf());

        double delta = _delta0;
        int maxnsing = std::min(f.size(),x.size());
        int nsing = maxnsing;

        dbg<<"iter   |f|inf   Q   |g|inf   delta\n";
        for(int k=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<f.normInf()<<"   "<<Q<<"   "<<g.normInf()<<"   "<<delta<<std::endl;
            xxdbg<<"f = "<<f<<std::endl;
            xxdbg<<"g = "<<g<<std::endl;
            xxdbg<<"J = "<<J<<std::endl;
            if (_shouldUseSvd)
                xxdbg<<"J.svd = "<<J.svd().getS().diag()<<std::endl;
            if (_shouldUseSvd && nsing == maxnsing && nsing > 1 &&
                J.svd().isSingular()) {
                xdbg<<"Singular J, so try lowering number of singular values.\n";
                nsing = J.svd().getKMax();
                xdbg<<"J Singular values = \n"<<J.svd().getS().diag()<<std::endl;
                xdbg<<"nsing -> "<<nsing<<std::endl;
            }
            h = -f/J;
            xdbg<<"h_newton = "<<h<<std::endl;

            double normsqg = g.normSq();
            double normH = h.norm();
            double normH1 = normH;
            double rhoDenom;

            if (normH <= delta) {
                xdbg<<"|h| < delta\n";
                rhoDenom = Q;
                xdbg<<"rhodenom = "<<rhoDenom<<std::endl;
            } else {
                xxdbg<<"normsqg = "<<normsqg<<std::endl;
                xxdbg<<"(J*g) = "<<J*g<<std::endl;
                xxdbg<<"normsq = "<<(J*g).normSq()<<std::endl;
                double alpha = normsqg / (J*g).normSq();
                xdbg<<"alpha = "<<alpha<<std::endl;
                double normG = sqrt(normsqg);
                xxdbg<<"normG = "<<normG<<std::endl;
                if (normG >= delta / alpha) {
                    xdbg<<"|g| > delta/alpha\n";
                    h = -(delta / normG) * g;
                    xdbg<<"h_gradient = "<<h<<std::endl;
                    rhoDenom = delta*(2.*alpha*normG-delta)/(2.*alpha);
                    xdbg<<"rhodenom = "<<rhoDenom<<std::endl;
                } else {
                    xdbg<<"dogleg\n";
                    temp = h + alpha*g;
                    double a = temp.normSq();
                    double b = -alpha * g * temp;
                    double c = alpha*alpha*g.normSq()-delta*delta;
                    // beta is the solution of 0 = a beta^2 + 2b beta + c
                    xdbg<<"a,b,c = "<<a<<" "<<b<<" "<<c<<std::endl;
                    double beta = (b <= 0) ?
                        (-b + sqrt(b*b - a*c)) / a :
                        -c / (b + sqrt(b*b - a*c));
                    xdbg<<"beta = "<<beta<<std::endl;
                    h = -alpha*g + beta*temp;
                    xdbg<<"h_dogleg = "<<h<<std::endl;
                    xdbg<<"h.norm() = "<<h.norm()<<"  delta = "<<delta<<std::endl;
                    rhoDenom = 
                        0.5*alpha*std::pow((1.-beta)*normG,2) +
                        beta*(2.-beta)*Q;
                    xdbg<<"rhodenom = "<<rhoDenom<<std::endl;
                }
                normH = h.norm();
            }

            CHECKSTEP(normH);

            xNew = x + h;
            xdbg<<"xnew = "<<xNew<<std::endl;
            this->calculateF(xNew,fNew);
            xdbg<<"fnew = "<<fNew<<std::endl;
            double QNew = 0.5*fNew.normSq();
            xdbg<<"Qnew = "<<QNew<<std::endl;

            bool deltaok = false;
            if (QNew < Q) {
                double rho = (Q-QNew) / rhoDenom;
                xdbg<<"rho = "<<Q-QNew<<" / "<<rhoDenom<<" = "<<rho<<std::endl;
                x = xNew; f = fNew; Q = QNew;
                CHECKF(f.normInf());
                this->calculateJ(x,f,J);
                J.unsetDiv();
                g = J.transpose() * f;
                xdbg<<"g = "<<g<<std::endl;
                CHECKG(g.normInf());
                if (rho > 0.75) {
                    delta = std::max(delta,3.*normH);
                    deltaok = true;
                }
            }
            if (deltaok) {
                nsing = maxnsing;
            } else {
                double normsqh = normH1*normH1;
                if (_shouldUseSvd && 
                    delta < normH1 && 
                    normsqg < 0.01 * normsqh && nsing > 1) {

                    dbg<<"normsqg == "<<normsqg/normsqh<<
                        " * normsqh, so try lowering number of singular values.\n";
                    --nsing;
                    dbg<<"nsing -> "<<nsing<<std::endl;
                    dbg<<"J Singular values = \n"<<J.svd().getS().diag()<<std::endl;
                    J.svd().top(nsing);
                } else {
                    delta /= 2.;
                    double min_delta = _minStep * (x.norm()+_minStep);
                    if (delta < min_delta) {
                        dbg<<"delta became too small ("<<
                            delta<<" < "<<min_delta<<")\n";
                        SHOWFAILFG; 
                        return false;
                    }
                }
            }
        }
        dbg<<"Maximum iterations exceeded in Dogleg method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solveHybrid(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // This is the Hybrid method which starts with the L-M method,
        // but switches to a quasi-newton method if ||f|| isn't approaching 0.
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());

        dbg<<"Start Solve_Hybrid\n";
        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        tmv::Vector<double> h(x.size());
        tmv::Vector<double> xNew(x.size());
        tmv::Vector<double> fNew(f.size());
        tmv::Vector<double> gNew(x.size());
        tmv::Matrix<double> JNew(f.size(),x.size());
        tmv::Vector<double> y(x.size());
        tmv::Vector<double> v(x.size());

        xdbg<<"x = "<<x<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f<<std::endl;
        double normInfF = f.normInf();
        CHECKF(normInfF);
        double Q = 0.5*f.normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        if (_shouldUseSvd) J.divideUsing(tmv::SV);
        J.saveDiv();
        xdbg<<"J = "<<J<<std::endl;

        tmv::SymMatrix<double> A = J.transpose()*J;
        xdbg<<"A = "<<A<<std::endl;
        if (_shouldUseSvd) A.divideUsing(tmv::SV);
        else if (_shouldUseCh) A.divideUsing(tmv::CH);
        else A.divideUsing(tmv::LU);
        A.saveDiv();
        tmv::SymMatrix<double> H(x.size());
        if (_shouldUseSvd) H.divideUsing(tmv::SV);
        else if (_shouldUseCh) H.divideUsing(tmv::CH);
        else H.divideUsing(tmv::LU);
        H.saveDiv();
        bool shouldUseDirectH = _hasDirectH;
        xdbg<<"shouldUseDirectH = "<<shouldUseDirectH<<std::endl;
        if (shouldUseDirectH) {
#ifndef NOTHROW
            try {
#endif
                xdbg<<"try calculateH\n";
                this->calculateH(x,f,J,H);
#ifndef NOTHROW
            } catch(NoDefinedH) {
                dbg<<"No direct H calculation - calculate on the fly\n";
                shouldUseDirectH = false;
                H.setToIdentity();
            }
#endif
        } else {
            xdbg<<"setToIdent\n";
            H.setToIdentity();
        }
        xdbg<<"After calculate H = "<<H<<std::endl;

        tmv::Vector<double> g = J.transpose() * f;
        xdbg<<"g = "<<g<<std::endl;
        double normInfG = g.normInf();
        CHECKG(normInfG);

        double mu = _tau * A.diag().normInf();
        A += mu;
        double nu = 2.;
        double delta = _delta0;
        bool shouldUseQuasiNewton = false;
        int count = 0;

        dbg<<"iter   |f|inf   Q   |g|inf   mu   delta  LM/QN\n";
        for(int k=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<normInfF<<"   "<<Q<<"   "<<normInfG<<"   "<< mu<<"   "<<
                delta<<"   "<<(shouldUseQuasiNewton?"QN":"LM")<<std::endl;
            xdbg<<"k = "<<k<<std::endl;
            xdbg<<"mu = "<<mu<<std::endl;
            xdbg<<"delta = "<<delta<<std::endl;
            xdbg<<"A = "<<A<<std::endl;
            xdbg<<"H = "<<H<<std::endl;
            xdbg<<"method = "<<(shouldUseQuasiNewton ? "quasinewton\n" : "LM\n");
            bool isBetter = false;
            bool shouldSwitchMethod = false;

            if (shouldUseQuasiNewton) {
#ifndef NOTHROW
                try { 
#endif
                    h = -g/H; 
#ifndef NOTHROW
                } catch (tmv::NonPosDef) {
                    xdbg<<"NonPosDef caught - switching division to LU method for H\n";
                    H.divideUsing(tmv::LU); 
                    h = -g/H;
                }
#endif
            } else {
#ifndef NOTHROW
                try { 
#endif
                    h = -g/A; 
#ifndef NOTHROW
                } catch (tmv::NonPosDef) {
                    xdbg<<"NonPosDef caught - switching division to LU method for A\n";
                    A.divideUsing(tmv::LU);
                    h = -g/A; 
                }
#endif
            }

            xdbg<<"h = "<<h<<std::endl;
            double normH = h.norm();
            CHECKSTEP(normH);
            if (shouldUseQuasiNewton && normH > delta) h *= delta/normH;

            xNew = x + h;
            xdbg<<"xnew = "<<xNew<<std::endl;
            this->calculateF(xNew,fNew);
            xdbg<<"fnew = "<<fNew<<std::endl;
            double QNew = 0.5*fNew.normSq();
            xdbg<<"Qnew = "<<QNew<<std::endl;

            double normInfGNew = 0.;
            bool isJNewSet = false;
            bool isGNewSet = false;
            if (!shouldUseDirectH || shouldUseQuasiNewton || QNew < Q) {
                this->calculateJ(xNew,fNew,JNew);
                xdbg<<"Jnew = "<<JNew<<std::endl;
                isJNewSet = true;
            }
            if (shouldUseQuasiNewton || QNew < Q) {
                if (!isJNewSet) dbg<<"Error: JNew should be set!\n";
                gNew = JNew.transpose() * fNew;
                xdbg<<"gnew = "<<gNew<<std::endl;
                normInfGNew = gNew.normInf();
                xdbg<<"NormInf(gnew) = "<<NormInf(gNew)<<std::endl;
                isGNewSet = true;
            }

            if (shouldUseQuasiNewton) {
                xdbg<<"quasinewton\n";
                if (!isGNewSet) dbg<<"Error: gNew should be set!\n";
                isBetter = 
                    (QNew < Q) || 
                    (QNew <= (1.+sqrtEps)*Q && normInfGNew < normInfG);
                xdbg<<"better = "<<isBetter<<std::endl;
                shouldSwitchMethod = (normInfGNew >= normInfG);
                xdbg<<"switchmethod = "<<shouldSwitchMethod<<std::endl;
                if (QNew < Q) {
                    double rho = (Q-QNew) / (-h*g-0.5*(J*h).normSq());
                    if (rho > 0.75) {
                        delta = std::max(delta,3.*normH);
                    } else if (rho < 0.25) {
                        delta /= 2.;
                        double min_delta = _minStep * (x.norm()+_minStep);
                        if (delta < min_delta) {
                            dbg<<"delta became too small ("<<
                                delta<<" < "<<min_delta<<")\n";
                            SHOWFAILFG; 
                            return false;
                        }
                    }
                } else {
                    delta /= 2.;
                    double min_delta = _minStep * (x.norm()+_minStep);
                    if (delta < min_delta) {
                        dbg<<"delta became too small ("<<
                            delta<<" < "<<min_delta<<")\n";
                        SHOWFAILFG; 
                        return false;
                    }
                }
            } else {
                xdbg<<"LM\n";
                if (QNew < Q) {
                    isBetter = true;
                    // we don't need the g vector anymore, so use this space
                    // to calculate g-mu*h
                    //double rho = (Q-QNew) / (0.5*h*(mu*h-g));
                    g -= mu*h;
                    double rho = (Q-QNew) / (-0.5*h*g);
                    mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
                    if (!isGNewSet) dbg<<"Error: gNew should be set!\n";
                    xdbg<<"check1: "<<normInfGNew<<" <? "<<0.02*QNew<<std::endl;
                    xdbg<<"check2: "<<Q-QNew<<" <? "<<0.02*QNew<<std::endl;
                    if (std::min(normInfGNew,Q-QNew) < 0.02 * QNew) {
                        ++count;
                        if (count == 3) shouldSwitchMethod = true;
                    } else {
                        count = 0;
                    }
                    if (count != 3) {
                        if (!isJNewSet) dbg<<"Error: JNew should be set!\n";
                        A = JNew.transpose() * JNew;
                        A += mu;
                    }
                } else {
                    A += mu*(nu-1.); mu *= nu; nu *= 2.;
                    count = 0;
                    // MJ: try this?
                    shouldSwitchMethod = (nu >= 32.);
                }
                A.unsetDiv();
                xdbg<<"better = "<<isBetter<<std::endl;
                xdbg<<"switchmethod = "<<shouldSwitchMethod<<std::endl;
                xdbg<<"count = "<<count<<std::endl;
            }

            if (!shouldUseDirectH) {
                if (!isJNewSet) dbg<<"Error: JNew should be set!\n";
                y = JNew.transpose()*(JNew*h) + (JNew-J).transpose()*fNew;
                double hy = h*y;
                xdbg<<"hy = "<<hy<<std::endl;
                if (hy > 0.) {
                    v = H*h;
                    xdbg<<"v = "<<v<<std::endl;
                    xdbg<<"y = "<<y<<std::endl;
                    xdbg<<"hv = "<<h*v<<std::endl;
                    H -= (1./(h*v)) * (v^v);
                    xdbg<<"H -> "<<H<<std::endl;
                    H += (1./hy) * (y^y);
                    H.unsetDiv();
                    xdbg<<"H -> "<<H<<std::endl;
                }
            }

            if (isBetter) {
                xdbg<<"better"<<std::endl;
                x = xNew; f = fNew; Q = QNew; normInfF = f.normInf(); 
                if (isJNewSet) J = JNew;
                else this->calculateJ(x,f,J);
                if (isGNewSet) { g = gNew; normInfG = normInfGNew; }
                else { g = J.transpose() * f; normInfG = g.normInf(); }
                J.unsetDiv();
                if (shouldUseDirectH && shouldUseQuasiNewton && !shouldSwitchMethod)
                    this->calculateH(x,f,J,H);
                CHECKF(normInfF);
                CHECKG(normInfG);
            }
            if (shouldSwitchMethod) {
                if (shouldUseQuasiNewton) {
                    xdbg<<"switch to LM\n";
                    A = J.transpose() * J;
                    //mu = _tau * A.diag().normInf();
                    A += mu;
                    A.unsetDiv();
                    shouldUseQuasiNewton = false;
                    count = 0;
                } else {
                    xdbg<<"switch to quasinewton\n";
                    delta = std::max(1.5*_minStep*(x.norm()+_minStep),0.2*normH);
                    if (shouldUseDirectH) {
                        this->calculateH(x,f,J,H);
                        H.unsetDiv();
                    }
                    shouldUseQuasiNewton = true;
                }
            }
        }
        dbg<<"Maximum iterations exceeded in Hybrid method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solveSecantLM(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // This is the Secant version of the Levenberg-Marquardt method
    {
        dbg<<"Start Solve_SecantLM\n";
        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        tmv::Vector<double> h(x.size());
        tmv::Vector<double> xNew(x.size());
        tmv::Vector<double> fNew(f.size());
        tmv::Vector<double> gNew(x.size());

        xdbg<<"x = "<<x<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f<<std::endl;
        CHECKF(f.normInf());
        double Q = 0.5*f.normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        if (_shouldUseSvd) J.divideUsing(tmv::SV);
        xdbg<<"J = "<<J<<std::endl;
        tmv::SymMatrix<double> A = J.transpose() * J;
        if (_shouldUseSvd) A.divideUsing(tmv::SV);
        else if (_shouldUseCh) A.divideUsing(tmv::CH);
        else A.divideUsing(tmv::LU);
        tmv::Vector<double> g = J.transpose() * f;
        xdbg<<"g = "<<g<<std::endl;
        CHECKG(g.normInf());

        double mu = _tau * A.diag().normInf();
        A += mu;
        double nu = 2.;

        dbg<<"iter   |f|inf   Q   |g|inf   mu\n";
        for(int k=0,j=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<f.normInf()<<"   "<<Q<<"   "<<g.normInf()<<"   "<<
                mu<<std::endl;
            xdbg<<"k = "<<k<<std::endl;
            xdbg<<"mu = "<<mu<<std::endl;
            xdbg<<"J = "<<J<<std::endl;
#ifndef NOTHROW
            try {
#endif
                h = -g/A;
#ifndef NOTHROW
            } catch (tmv::NonPosDef) {
                xdbg<<"NonPosDef caught - switching division to LU method.\n";
                A.divideUsing(tmv::LU);
                h = -g/A;
            }
#endif
            xdbg<<"h = "<<h<<std::endl;
            double normH = h.norm();
            CHECKSTEP(normH);

            xdbg<<"j = "<<j<<std::endl;
            if (h(j) < 0.8 * normH) {
                xNew = x; 
                double eta = _minStep * (x.norm() + 1.);
                xNew(j) += eta;
                this->calculateF(xNew,fNew);
                J.col(j) = (fNew-f)/eta;
                xdbg<<"J -> "<<J<<std::endl;
            }
            j = (j+1)%J.ncols();

            xNew = x + h;
            xdbg<<"xnew = "<<xNew<<std::endl;
            this->calculateF(xNew,fNew);
            xdbg<<"fnew = "<<fNew<<std::endl;
            double QNew = 0.5*fNew.normSq();
            xdbg<<"Qnew = "<<QNew<<std::endl;
            J += (1./h.normSq()) * ((fNew - f - J*h) ^ h);
            xdbg<<"J -> "<<J<<std::endl;

            if (QNew < Q) {
                x = xNew; f = fNew; 
                CHECKF(f.normInf());

                A = J.transpose() * J;
                gNew = J.transpose() * f;
                CHECKG(g.normInf());

                g -= mu*h;
                double rho = (Q-QNew) / (-0.5*h*g);
                xdbg<<"rho = "<<Q-QNew<<" / "<<(-0.5*h*g)<<" = "<<rho<<std::endl;
                mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
                xdbg<<"mu = "<<mu<<std::endl;
                A += mu;
                Q = QNew; g = gNew;
            } else {
                A += mu*(nu-1.); mu *= nu; nu *= 2.;
            }
        }
        dbg<<"Maximum iterations exceeded in Secant LM method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solveSecantDogleg(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // This is the Secant version of the Dogleg method
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());

        dbg<<"Start Solve_SecantDogleg\n";
        _pJ.reset(new tmv::Matrix<double>(f.size(),x.size()));
        tmv::Matrix<double>& J = *_pJ;
        tmv::Vector<double> h(x.size());
        tmv::Vector<double> temp(x.size());
        tmv::Vector<double> xNew(x.size());
        tmv::Vector<double> fNew(f.size());
        tmv::Vector<double> y(f.size());
        tmv::Vector<double> djodjy(f.size());

        xdbg<<"x = "<<x<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f<<std::endl;
        CHECKF(f.normInf());
        double Q = 0.5*f.normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        if (_shouldUseSvd) J.divideUsing(tmv::SV);
        tmv::Matrix<double> D = J.inverse();

        tmv::Vector<double> g = J.transpose() * f;
        xdbg<<"g = "<<g<<std::endl;
        CHECKG(g.normInf());
        double delta = _delta0;

        dbg<<"iter   |f|inf   Q   |g|inf   delta\n";
        for(int k=0,j=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<f.normInf()<<"   "<<Q<<"   "<<g.normInf()<<"   "<<
                delta<<std::endl;
            h = -D*f;
            xdbg<<"h = "<<h<<std::endl;

            double normsqg = g.normSq();
            double alpha = normsqg / (J*g).normSq();
            double normH = h.norm();
            double rhoDenom;

            if (normH <= delta) {
                xdbg<<"|h| < delta \n";
                rhoDenom = Q;
            } else {
                double normG = sqrt(normsqg);
                if (normG >= delta / alpha) {
                    xdbg<<"|g| > delta/alpha \n";
                    h = -(delta / normG) * g;
                    xdbg<<"h = "<<h<<std::endl;
                    rhoDenom = delta*(2.*alpha*normG-delta)/(2.*alpha);
                } else {
                    xdbg<<"dogleg\n";
                    temp = h + alpha*g;
                    double a = temp.normSq();
                    double b = -alpha * g * temp;
                    double c = alpha*alpha*g.normSq()-delta*delta;
                    // beta is the solution of 0 = a beta^2 + 2b beta + c
                    double beta = (b <= 0) ?
                        (-b + sqrt(b*b - a*c)) / a :
                        -c / (b + sqrt(b*b - a*c));
                    xdbg<<"alpha = "<<alpha<<std::endl;
                    xdbg<<"beta = "<<beta<<std::endl;
                    h = -alpha*g + beta*temp;
                    xdbg<<"h = "<<h<<std::endl;
                    rhoDenom = 
                        0.5*alpha*std::pow((1.-beta)*normG,2) + 
                        beta*(2.-beta)*Q;
                }
                normH = h.norm();
            }

            CHECKSTEP(normH);

            bool resetd = false;
            if (h(j) < 0.8 * normH) {
                xNew = x; 
                double eta = _minStep * (x.norm() + 1.);
                xNew(j) += eta;
                this->calculateF(xNew,fNew);
                y = fNew-f;
                J.col(j) = y/eta;
                double djy = D.row(j)*y;
                if (djy < sqrtEps*eta) {
                    resetd = true;
                } else {
                    djodjy = D.row(j)/djy;
                    D -= ((D*y) ^ djodjy);
                    D.row(j) += eta*djodjy;
                }
            }
            j = (j+1)%J.ncols();

            xNew = x + h;
            this->calculateF(xNew,fNew);
            double QNew = 0.5*fNew.normSq();

            y = fNew - f;
            J += (1./h.normSq()) * ((fNew - f - J*h) ^ h);
            double hDy = h*D*y;
            if (resetd || hDy < sqrtEps*h.norm()) {
                D = J.inverse();
            } else {
                D += 1./(hDy) * ((h-D*y) ^ (h*D));
            }

            if (QNew < Q) {
                double rho = (Q-QNew) / rhoDenom;
                xdbg<<"rho = "<<Q-QNew<<" / "<<rhoDenom<<" = "<<rho<<std::endl;
                x = xNew; f = fNew; Q = QNew;
                CHECKF(f.normInf());
                g = J.transpose() * f;
                xdbg<<"g = "<<g<<std::endl;
                CHECKG(g.normInf());
                if (rho > 0.75) {
                    delta = std::max(delta,3.*normH);
                } else if (rho < 0.25) {
                    delta /= 2.;
                    double min_delta = _minStep * (x.norm()+_minStep);
                    if (delta < min_delta) {
                        dbg<<"delta became too small ("<<
                            delta<<" < "<<min_delta<<")\n";
                        SHOWFAILFG; 
                        return false;
                    }
                }
            } else {
                delta /= 2.;
                double min_delta = _minStep * (x.norm()+_minStep);
                if (delta < min_delta) {
                    dbg<<"delta became too small ("<<delta<<" < "<<min_delta<<")\n";
                    SHOWFAILFG; 
                    return false;
                }
            }
        }
        dbg<<"Maximum iterations exceeded in Secant Dogleg method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solve(
        tmv::Vector<double>& x, tmv::Vector<double>& f) const
        // On input, x is the initial guess
        // On output, if return is true, then
        // x is the solution for which either f.norm() ~= 0
        // or f is a local minimum.
    {
#ifndef NOTHROW
        try {
#endif
            switch (_method) {
              case HYBRID : return solveHybrid(x,f);
              case DOGLEG : return solveDogleg(x,f);
              case LM : return solveLM(x,f);
              case NEWTON : return solveNewton(x,f);
              case SECANT_LM : return solveSecantLM(x,f);
              case SECANT_DOGLEG : return solveSecantDogleg(x,f);
              default : dbg<<"Unknown method\n"; return false;
            }
#ifndef NOTHROW
        } 
#if 0
        catch (int) {}
#else
        catch (tmv::Singular& e) {
            dbg<<"Singular matrix encountered in NLSolver::Solve\n";
            dbg<<e<<std::endl;
        } catch (tmv::Error& e) {
            dbg<<"TMV error encountered in NLSolver::Solve\n";
            dbg<<e<<std::endl;
        } catch (...) {
            dbg<<"Error encountered in NLSolver::Solve\n";
        }
#endif
        return false;
#endif
    }

    void NLSolver::getCovariance(tmv::Matrix<double>& cov) const
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
        if (!_pJ.get()) {
            throw std::runtime_error(
                "J not set before calling getCovariance");
        }
        tmv::Matrix<double>& J = *_pJ;
        // This might have changed between solve and getCovariance:
        // And we need to set the threshold to sqrt(eps) rather than eps
        if (_shouldUseSvd) {
            J.divideUsing(tmv::SV);
            J.svd().thresh(sqrtEps); 
        }
        J.makeInverseATA(cov);
        //dbg<<"getCovariance:\n";
        //dbg<<"J = "<<J<<std::endl;
        //dbg<<"JtJ = "<<J.adjoint()*J<<std::endl;
        //dbg<<"(JtJ)^-1 = "<<(J.adjoint()*J).inverse()<<std::endl;
        //dbg<<"cov = "<<cov<<std::endl;
    }

    void NLSolver::getInverseCovariance(tmv::Matrix<double>& invCov) const
    {
        if (!_pJ.get()) {
            throw std::runtime_error(
                "J not set before calling getInverseCovariance");
        }
        tmv::Matrix<double>& J = *_pJ;
        invCov = J.transpose() * J;
    }

#else // USE_EIGEN

    NLSolver::NLSolver() : 
        _method(HYBRID),
        _fTol(1.e-8), _gTol(1.e-8), _minStep(1.e-8), _maxIter(200),
        _tau(1.e-3), _delta0(1.), 
        _nlOut(0), _verbose(0)
    {}

    void NLSolver::calculateJ(
        const DVector& x, const DVector& f, DMatrix& df) const
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());
        // Do a finite difference calculation for J.
        // This function is virtual, so if there is a better way to 
        // calculate J, then you should override this version.

        DVector x2 = x;
        DVector f2(f.size());
        DVector f1(f.size());
        const int n = x.size();
        for(int j=0;j<n;++j) {
            const double dx = sqrtEps * (x.norm() + 1.);
            x2(j) += dx;
            this->calculateF(x2,f2);
            x2(j) -= 2.*dx;
            this->calculateF(x2,f1);
            df.col(j) = (f2-f1)/(2.*dx);
            x2(j) = x(j);
        }
    }

    bool NLSolver::testJ(
        const DVector& x, DVector& f,
        std::ostream* os, double relErr) const
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());

        this->calculateF(x,f);
        _pJ.reset(new DMatrix(f.size(),x.size()));
        DMatrix& J = *_pJ;
        this->calculateJ(x,f,J);
        DMatrix Jn(f.size(),x.size());
        NLSolver::calculateJ(x,f,Jn);
        double err = (J-Jn).TMV_maxAbsElement() / Jn.norm();
        if (!relErr) relErr = 10.*sqrtEps;
        if (os) {
            *os << "TestJ:\n";
            if (_verbose >= 1) {
                *os << "x = "<<x<<std::endl;
                *os << "f = "<<f<<std::endl;
                *os << "Direct J = "<<J<<std::endl;
                *os << "Numeric J = "<<Jn<<std::endl;
            }
            *os << "MaxAbsElement(J-J_num) / J.norm() = "<<err<<std::endl;
            *os << "cf. relerr = "<<relErr<<std::endl;
            if (err >= relErr) {
                DMatrix diff = J-Jn;
                *os << "J-J_num = "<<diff;
                double maxel = diff.TMV_maxAbsElement();
                *os << "Max element = "<<maxel<<std::endl;
                const int m = diff.TMV_colsize();
                const int n = diff.TMV_rowsize();
                for(int i=0;i<m;++i) {
                    for(int j=0;j<n;++j) {
                        if (std::abs(diff(i,j)) > 0.9*maxel) {
                            *os<<"J("<<i<<','<<j<<") = "<<J(i,j)<<"  ";
                            *os<<"J_num("<<i<<','<<j<<") = "<<Jn(i,j)<<"  ";
                            *os<<"diff = "<<J(i,j)-Jn(i,j)<<std::endl;
                        }
                    }
                }
            }
        }
        return err < relErr;
    }

#define CHECKF(normInfF) \
    do { \
        double checkfTemp = (normInfF); \
        if (!(checkfTemp > _fTol)) { \
            dbg<<"Found ||f|| ~= 0\n"; \
            dbg<<"||f||_inf = "<<checkfTemp<<" < "<<_fTol<<std::endl; \
            return true; \
        } \
    } while (false)

#define CHECKG(normInfG) \
    do { \
        double checkgTemp = (normInfG); \
        if (!(checkgTemp > _gTol)) { \
            dbg<<"Found local minimum of ||f||\n"; \
            dbg<<"||g||_inf = "<<checkgTemp<<" < "<<_gTol<<std::endl; \
            return true; \
        } \
    } while (false)

#define SHOWFAILFG \
    do { \
        dbg<<"||f||_inf = "<<f.TMV_normInf()<<" !< "<<_fTol<<std::endl; \
        dbg<<"||g||_inf = "<<g.TMV_normInf()<<" !< "<<_gTol<<std::endl; \
    } while (false)

#define CHECKSTEP(normH) \
    do { \
        double checkStepTemp1 = (normH); \
        double checkStepTemp2 = _minStep*(x.norm()+_minStep); \
        if (!(checkStepTemp1 > checkStepTemp2)) { \
            dbg<<"Step size became too small\n"; \
            dbg<<"||h|| = "<<checkStepTemp1<<" < "<<checkStepTemp2<<std::endl; \
            SHOWFAILFG; \
            return false; \
        } \
    } while (false)

    bool NLSolver::solveDogleg(DVector& x, DVector& f) const
        // This is the Dogleg method
    {
        dbg<<"Start Solve_Dogleg\n";
        _pJ.reset(new DMatrix(f.size(),x.size()));
        DMatrix& J = *_pJ;
        DVector h(x.size());
        DVector temp(x.size());
        DVector xNew(x.size());
        DVector fNew(f.size());

        xdbg<<"x = "<<x.transpose()<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f.transpose()<<std::endl;
        CHECKF(f.TMV_normInf());
        double Q = 0.5*f.TMV_normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        xdbg<<"J = "<<J<<std::endl;

        DVector g = J.transpose() * f;
        xdbg<<"g = "<<g.transpose()<<std::endl;
        CHECKG(g.TMV_normInf());

        double delta = _delta0;
        int maxnsing = std::min(f.size(),x.size());
        int nsing = maxnsing;

        dbg<<"iter   |f|inf   Q   |g|inf   delta\n";
        for(int k=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<f.TMV_normInf()<<"   "<<Q<<"   "<<g.TMV_normInf()<<"   "<<delta<<std::endl;
            //h = -f/J;
            J.lu().solve(-f,&h);
            xdbg<<"h = "<<h.transpose()<<std::endl;

            double normsqg = g.TMV_normSq();
            double normH = h.norm();
            double rhoDenom;

            if (normH <= delta) {
                xdbg<<"|h| < delta\n";
                rhoDenom = Q;
                xdbg<<"rhodenom = "<<rhoDenom<<std::endl;
            } else {
                double alpha = normsqg / (J*g).TMV_normSq();
                double normG = sqrt(normsqg);
                if (normG >= delta / alpha) {
                    xdbg<<"|g| > delta/alpha\n";
                    h = -(delta / normG) * g;
                    xdbg<<"h = "<<h.transpose()<<std::endl;
                    rhoDenom = delta*(2.*alpha*normG-delta)/(2.*alpha);
                    xdbg<<"rhodenom = "<<rhoDenom<<std::endl;
                } else {
                    xdbg<<"dogleg\n";
                    temp = h + alpha*g;
                    double a = temp.TMV_normSq();
                    double b = -alpha * (g.transpose() * temp)(0,0);
                    double c = alpha*alpha*g.TMV_normSq()-delta*delta;
                    // beta is the solution of 0 = a beta^2 + 2b beta + c
                    xdbg<<"a,b,c = "<<a<<" "<<b<<" "<<c<<std::endl;
                    double beta = (b <= 0) ?
                        (-b + sqrt(b*b - a*c)) / a :
                        -c / (b + sqrt(b*b - a*c));
                    xdbg<<"alpha = "<<alpha<<std::endl;
                    xdbg<<"beta = "<<beta<<std::endl;
                    h = -alpha*g + beta*temp;
                    xdbg<<"h = "<<h.transpose()<<std::endl;
                    xdbg<<"h.norm() = "<<h.norm()<<"  delta = "<<delta<<std::endl;
                    rhoDenom = 
                        0.5*alpha*std::pow((1.-beta)*normG,2) +
                        beta*(2.-beta)*Q;
                    xdbg<<"rhodenom = "<<rhoDenom<<std::endl;
                }
                normH = h.norm();
            }

            CHECKSTEP(normH);

            xNew = x + h;
            xdbg<<"xnew = "<<xNew.transpose()<<std::endl;
            this->calculateF(xNew,fNew);
            xdbg<<"fnew = "<<fNew.transpose()<<std::endl;
            double QNew = 0.5*fNew.TMV_normSq();
            xdbg<<"Qnew = "<<QNew<<std::endl;

            bool deltaok = false;
            if (QNew < Q) {
                double rho = (Q-QNew) / rhoDenom;
                xdbg<<"rho = "<<Q-QNew<<" / "<<rhoDenom<<" = "<<rho<<std::endl;
                x = xNew; f = fNew; Q = QNew;
                CHECKF(f.TMV_normInf());
                this->calculateJ(x,f,J);
                g = J.transpose() * f;
                xdbg<<"g = "<<g.transpose()<<std::endl;
                CHECKG(g.TMV_normInf());
                if (rho > 0.75) {
                    delta = std::max(delta,3.*normH);
                    deltaok = true;
                }
            }
            if (deltaok) {
                nsing = maxnsing;
            } else {
                delta /= 2.;
                double min_delta = _minStep * (x.norm()+_minStep);
                if (delta < min_delta) {
                    dbg<<"delta became too small ("<<
                        delta<<" < "<<min_delta<<")\n";
                    SHOWFAILFG; 
                    return false;
                }
            }
        }
        dbg<<"Maximum iterations exceeded in Dogleg method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solveHybrid(DVector& x, DVector& f) const
        // This is the Hybrid method which starts with the L-M method,
        // but switches to a quasi-newton method if ||f|| isn't approaching 0.
    {
        const double sqrtEps = sqrt(std::numeric_limits<double>::epsilon());

        dbg<<"Start Solve_Hybrid\n";
        _pJ.reset(new DMatrix(f.size(),x.size()));
        DMatrix& J = *_pJ;
        DVector h(x.size());
        DVector xNew(x.size());
        DVector fNew(f.size());
        DVector gNew(x.size());
        DMatrix JNew(f.size(),x.size());
        DVector y(x.size());
        DVector v(x.size());

        xdbg<<"x = "<<x.transpose()<<std::endl;
        this->calculateF(x,f);
        xdbg<<"f = "<<f.transpose()<<std::endl;
        double normInfF = f.TMV_normInf();
        CHECKF(normInfF);
        double Q = 0.5*f.TMV_normSq();
        xdbg<<"Q = "<<Q<<std::endl;
        this->calculateJ(x,f,J);
        xdbg<<"J = "<<J<<std::endl;

        DMatrix A = J.transpose()*J;
        xdbg<<"A = "<<A<<std::endl;
        DMatrix H(x.size(),x.size());
        xdbg<<"setToIdent\n";
        H.TMV_setToIdentity();
        xdbg<<"After calculate H = "<<H<<std::endl;

        DVector g = J.transpose() * f;
        xdbg<<"g = "<<g.transpose()<<std::endl;
        double normInfG = g.TMV_normInf();
        CHECKG(normInfG);

        double mu = _tau * A.TMV_diag().TMV_normInf();
        A.EIGEN_diag() += mu;
        double nu = 2.;
        double delta = _delta0;
        bool shouldUseQuasiNewton = false;
        int count = 0;

        dbg<<"iter   |f|inf   Q   |g|inf   mu   delta  LM/QN\n";
        for(int k=0;k<_maxIter;++k) {
            dbg<<k<<"   "<<normInfF<<"   "<<Q<<"   "<<normInfG<<"   "<< mu<<"   "<<
                delta<<"   "<<(shouldUseQuasiNewton?"QN":"LM")<<std::endl;
            xdbg<<"k = "<<k<<std::endl;
            xdbg<<"mu = "<<mu<<std::endl;
            xdbg<<"delta = "<<delta<<std::endl;
            xdbg<<"A = "<<A<<std::endl;
            xdbg<<"H = "<<H<<std::endl;
            xdbg<<"method = "<<(shouldUseQuasiNewton ? "quasinewton\n" : "LM\n");
            bool isBetter = false;
            bool shouldSwitchMethod = false;

            if (shouldUseQuasiNewton) {
                //h = -g/H; 
                H.ldlt().solve(-g,&h);
            } else {
                //h = -g/A; 
                A.ldlt().solve(-g,&h);
            }

            xdbg<<"h = "<<h.transpose()<<std::endl;
            double normH = h.norm();
            CHECKSTEP(normH);
            if (shouldUseQuasiNewton && normH > delta) h *= delta/normH;

            xNew = x + h;
            xdbg<<"xnew = "<<xNew.transpose()<<std::endl;
            this->calculateF(xNew,fNew);
            xdbg<<"fnew = "<<fNew.transpose()<<std::endl;
            double QNew = 0.5*fNew.TMV_normSq();
            xdbg<<"Qnew = "<<QNew<<std::endl;

            this->calculateJ(xNew,fNew,JNew);
            xdbg<<"Jnew = "<<JNew<<std::endl;
            if (shouldUseQuasiNewton || QNew < Q) {
                gNew = JNew.transpose() * fNew;
                xdbg<<"gnew = "<<gNew.transpose()<<std::endl;
            }
            double normInfGNew = gNew.TMV_normInf();

            if (shouldUseQuasiNewton) {
                xdbg<<"quasinewton\n";
                isBetter = 
                    (QNew < Q) || 
                    (QNew <= (1.+sqrtEps)*Q && normInfGNew < normInfG);
                xdbg<<"better = "<<isBetter<<std::endl;
                shouldSwitchMethod = (normInfGNew >= normInfG);
                xdbg<<"switchmethod = "<<shouldSwitchMethod<<std::endl;
                if (QNew < Q) {
                    double rho = (Q-QNew) / (-(h.transpose()*g)(0,0)-0.5*(J*h).TMV_normSq());
                    if (rho > 0.75) {
                        delta = std::max(delta,3.*normH);
                    } else if (rho < 0.25) {
                        delta /= 2.;
                        double min_delta = _minStep * (x.norm()+_minStep);
                        if (delta < min_delta) {
                            dbg<<"delta became too small ("<<
                                delta<<" < "<<min_delta<<")\n";
                            SHOWFAILFG; 
                            return false;
                        }
                    }
                } else {
                    delta /= 2.;
                    double min_delta = _minStep * (x.norm()+_minStep);
                    if (delta < min_delta) {
                        dbg<<"delta became too small ("<<
                            delta<<" < "<<min_delta<<")\n";
                        SHOWFAILFG; 
                        return false;
                    }
                }
            } else {
                xdbg<<"LM\n";
                if (QNew <= Q) {
                    isBetter = true;
                    // we don't need the g vector anymore, so use this space
                    // to calculate g-mu*h
                    //double rho = (Q-QNew) / (0.5*h*(mu*h-g));
                    g -= mu*h;
                    double rho = (Q-QNew) / (-0.5*(h.transpose()*g)(0,0));
                    mu *= std::max(1./3.,1.-std::pow(2.*rho-1.,3)); nu = 2.;
                    xdbg<<"check1: "<<normInfGNew<<" <? "<<0.02*QNew<<std::endl;
                    xdbg<<"check2: "<<Q-QNew<<" <? "<<0.02*QNew<<std::endl;
                    if (std::min(normInfGNew,Q-QNew) < 0.02 * QNew) {
                        ++count;
                        if (count == 3) shouldSwitchMethod = true;
                    } else {
                        count = 0;
                    }
                    if (count != 3) {
                        A = JNew.transpose() * JNew;
                        A.EIGEN_diag() += mu;
                    }
                } else {
                    A.EIGEN_diag() += mu*(nu-1.); mu *= nu; nu *= 2.;
                    count = 0;
                    // MJ: try this?
                    shouldSwitchMethod = (nu >= 32.);
                }
                xdbg<<"better = "<<isBetter<<std::endl;
                xdbg<<"switchmethod = "<<shouldSwitchMethod<<std::endl;
                xdbg<<"count = "<<count<<std::endl;
            }

            y = JNew.transpose()*(JNew*h) + (JNew-J).transpose()*fNew;
            double hy = (h.transpose()*y)(0,0);
            xdbg<<"hy = "<<hy<<std::endl;
            if (hy > 0.) {
                v = H*h;
                xdbg<<"v = "<<v.transpose()<<std::endl;
                xdbg<<"y = "<<y.transpose()<<std::endl;
                double hv = (h.transpose()*v)(0,0);
                xdbg<<"hv = "<<hv<<std::endl;
                H -= (1./hv) * (v * v.transpose());
                xdbg<<"H -> "<<H<<std::endl;
                H += (1./hy) * (y * y.transpose());
                xdbg<<"H -> "<<H<<std::endl;
            }

            if (isBetter) {
                xdbg<<"better"<<std::endl;
                x = xNew; f = fNew; Q = QNew; J = JNew; g = gNew; 
                normInfF = f.TMV_normInf(); normInfG = normInfGNew;
                CHECKF(normInfF);
                CHECKG(normInfG);
            }
            if (shouldSwitchMethod) {
                if (shouldUseQuasiNewton) {
                    xdbg<<"switch to LM\n";
                    A = J.transpose() * J;
                    //mu = _tau * A.diag().normInf();
                    A.EIGEN_diag() += mu;
                    shouldUseQuasiNewton = false;
                    count = 0;
                } else {
                    xdbg<<"switch to quasinewton\n";
                    delta = std::max(1.5*_minStep*(x.norm()+_minStep),0.2*normH);
                    shouldUseQuasiNewton = true;
                }
            }
        }
        dbg<<"Maximum iterations exceeded in Hybrid method\n";
        SHOWFAILFG; 
        return false;
    }

    bool NLSolver::solve(DVector& x, DVector& f) const
    {
#ifndef NOTHROW
        try {
#endif
            switch (_method) {
              case HYBRID : return solveHybrid(x,f);
              case DOGLEG : return solveDogleg(x,f);
              default : dbg<<"Unknown method\n"; return false;
            }
#ifndef NOTHROW
        } catch (...) {
            dbg<<"Error encountered in NLSolver::Solve\n";
        }
#endif
        return false;
    }

    void NLSolver::getCovariance(DMatrix& cov) const
    {
        const double eps = std::numeric_limits<double>::epsilon();
        if (!_pJ.get()) {
            throw std::runtime_error(
                "J not set before calling getCovariance");
        }
        DMatrix& J = *_pJ;
        // This might have changed between solve and getCovariance:
        // And we need to set the threshold to sqrt(eps) rather than eps
        if (_shouldUseSvd) {
            Eigen::SVD<DMatrix> SV_Solver_J = J.svd();
            SV_Solver_J.sort();
            const DVector& svd_s = SV_Solver_J.singularValues();
            int kmax = svd_s.size();
            while(svd_s(kmax-1) < eps * svd_s(0)) --kmax;
            dbg<<"In NLSolver::getCovariance:\n";
            dbg<<"Using kmax = "<<kmax<<"   size = "<<svd_s.size()<<std::endl;
            //  (JtJ)^-1 = ( (USVt)t (USVt) )^-1
            //          = ( V St Ut U S Vt )^-1
            //          = ( V S^2 Vt )^-1
            //          = V S^-2 Vt
            const DMatrix& svd_v = SV_Solver_J.matrixV();
            DVector sm2 = svd_s.cwise().square().cwise().inverse();
            sm2.TMV_subVector(kmax,svd_s.size()).setZero();
            cov = svd_v * sm2.asDiagonal() * svd_v.transpose();
        } else {
            Eigen::QR<DMatrix> QR_Solver_J = J.qr();
            // (JtJ)^-1 = ( (QR)t (QR) )^-1
            //          = ( Rt Qt Q R ) ^-1
            //          = ( Rt R )^-1
            //          = R^-1 Rt^-1
            cov.setIdentity();
            QR_Solver_J.matrixR().transpose().solveTriangularInPlace(cov);
            QR_Solver_J.matrixR().solveTriangularInPlace(cov);
        }
    }

    void NLSolver::getInverseCovariance(DMatrix& invCov) const
    {
        if (!_pJ.get()) {
            throw std::runtime_error(
                "J not set before calling getInverseCovariance");
        }
        DMatrix& J = *_pJ;
        invCov = J.transpose() * J;
    }

#endif

}}}}
