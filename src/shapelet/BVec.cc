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
#include <iostream>

#include "lsst/meas/algorithms/shapelet/BVec.h"
#include "lsst/meas/algorithms/shapelet/BinomFact.h"
#include "lsst/afw/geom/Angle.h"

namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    BVec& BVec::operator=(const AssignableToBVec& rhs)
    {
        rhs.assignTo(*this);
        return *this; 
    }

    BVec& BVec::operator=(const BVec& rhs)
    {
        rhs.assignTo(*this);
        return *this; 
    }

    void BVec::assignTo(BVec& rhs) const
    {
        rhs.setSigma(_sigma);
        rhs.setValues(_b);
    }

    void BVec::setValues(const DVector& v)
    {
        if (_b.size() == v.size()) {
            _b = v;
        } else if (_b.size() > v.size()) {
            _b.TMV_subVector(0,v.size()) = v;
            _b.TMV_subVector(v.size(),_b.size()).setZero();
        } else {
            xdbg<<"Warning truncating information in BVec::setValues\n";
            _b = v.TMV_subVector(0,_b.size());
        }
    }

    void BVec::conjugateSelf()
    {
        const int N = getOrder();
        for(int n=1,k=1;n<=N;++n) {
            for(int m=n;m>=0;m-=2) {
                if (m==0) ++k;
                else {
                    _b(k+1) *= -1.0;
                    k+=2;
                }
            }
        }
    }

    void calculateZTransform(
        std::complex<double> z, int order1, int order2, DMatrix& T)
    {
        const int size1 = (order1+1)*(order1+2)/2;
        const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(T.TMV_colsize()) >= size1);
        Assert(int(T.TMV_rowsize()) >= size2);

        T.setZero();

        if (z == 0.0) { 
            T.TMV_diagpart(0,0,std::min(size1,size2)).TMV_setAllTo(1.);
            return; 
        }

        // T(st,pq) = f(p,s) f*(q,t)
        // f(p,0) = (-z/2)^p / sqrt(p!) exp(-|z|^2/8)
        // f(p,s+1) = (sqrt(p) f(p-1,s) + 1/2 z* f(p,s) )/sqrt(s+1)

        std::complex<double> zo2 = -z/2.;
        std::vector<std::vector<std::complex<double> > > f(
            order2+1,std::vector<std::complex<double> >(order1+1));

        f[0][0] = exp(-std::norm(z)/8.); // f(0,0)
        for(int p=1;p<=order2;++p) {
            f[p][0] = f[p-1][0]*(-zo2)/sqrtn(p); // f(p,0)
        }
        for(int s=0;s<order1;++s) {
            f[0][s+1] = std::conj(zo2)*f[0][s]/sqrtn(s+1); // f(0,s+1)
            for(int p=1;p<=order2;++p) {
                f[p][s+1] = (sqrtn(p)*f[p-1][s] + std::conj(zo2)*f[p][s])/
                    sqrtn(s+1); // f(p,s+1)
            }
        }

        for(int n=0,pq=0;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                double* Tpq = TMV_ptr(T.col(pq));
                double* Tpq1 = Tpq + TMV_stepj(T);
                const std::vector<std::complex<double> >& fp = f[p];
                const std::vector<std::complex<double> >& fq = f[q];
                for(int nn=0,st=0;nn<=order1;++nn) {
                    for(int s=nn,t=0;s>=t;--s,++t,++st) {
                        std::complex<double> t1 = fp[s] * std::conj(fq[t]);
                        if (p==q) {
                            if (s==t) {
                                Tpq[st] = std::real(t1); // T(st,pq)
                            } else {
                                Tpq[st] = std::real(t1);
                                Tpq[st+1] = std::imag(t1);
                                ++st;
                            }
                        } else if (s==t) {
                            // b_st = t_stpq b_pq + t_stqp b_qp
                            // In this case with s=t, t_stpq == t_stqp*
                            Tpq[st] = 2.*std::real(t1);
                            Tpq1[st] = -2.*std::imag(t1);
                        } else {
                            std::complex<double> t2 = fq[s] * std::conj(fp[t]);
                            Tpq[st] = std::real(t1) + std::real(t2);
                            Tpq[st+1]= std::imag(t1) + std::imag(t2);
                            Tpq1[st] = -std::imag(t1) + std::imag(t2);
                            Tpq1[st+1] = std::real(t1) - std::real(t2);
                            ++st;
                        }
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void augmentZTransformCols(std::complex<double> z, int order, DMatrix& T)
    {
        const int order1 = order;
        const int order2 = order+2;
        const int size1 = (order1+1)*(order1+2)/2;
        //const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(T.TMV_colsize()) >= size1);
        Assert(int(T.TMV_rowsize()) >= (order2+1)*(order2+2)/2);

        if (z == 0.0) return; 

        std::complex<double> zo2 = -z/2.;
        std::vector<std::vector<std::complex<double> > > f(
            order2+1,std::vector<std::complex<double> >(order1+1));

        f[0][0] = exp(-std::norm(z)/8.); // f(0,0)
        for(int p=1;p<=order2;++p) {
            f[p][0] = f[p-1][0]*(-zo2)/sqrtn(p); // f(p,0)
        }
        for(int s=0;s<order1;++s) {
            f[0][s+1] = std::conj(zo2)*f[0][s]/sqrtn(s+1); // f(0,s+1)
            for(int p=1;p<=order2;++p) {
                f[p][s+1] = (sqrtn(p)*f[p-1][s] + std::conj(zo2)*f[p][s])/
                    sqrtn(s+1); // f(p,s+1)
            }
        }

        for(int n=order1+1,pq=size1;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                const std::vector<std::complex<double> >& fp = f[p];
                const std::vector<std::complex<double> >& fq = f[q];
                double* Tpq = TMV_ptr(T.col(pq));
                double* Tpq1 = Tpq + TMV_stepj(T);
                for(int nn=0,st=0;nn<=order1;++nn) {
                    for(int s=nn,t=0;s>=t;--s,++t,++st) {
                        std::complex<double> t1 = fp[s] * std::conj(fq[t]);
                        if (p==q) {
                            if (s==t) {
                                Tpq[st] = std::real(t1); 
                            } else {
                                Tpq[st] = std::real(t1);
                                Tpq[st+1] = std::imag(t1);
                                ++st;
                            }
                        } else if (s==t) {
                            Tpq[st] = 2.*std::real(t1);
                            Tpq1[st] = -2.*std::imag(t1);
                        } else {
                            std::complex<double> t2 = fq[s] * std::conj(fp[t]);
                            Tpq[st] = std::real(t1) + std::real(t2);
                            Tpq[st+1]= std::imag(t1) + std::imag(t2);
                            Tpq1[st] = -std::imag(t1) + std::imag(t2);
                            Tpq1[st+1] = std::real(t1) - std::real(t2);
                            ++st;
                        }
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void applyZ(std::complex<double> z, BVec& b)
    {
        if (z != 0.0) {
            DMatrix T(int(b.size()),int(b.size()));
            calculateZTransform(z/b.getSigma(),b.getOrder(),T);
            b.vec() = T * b.vec();
        }
    }

    void calculateMuTransform(double mu, int order1, int order2, DMatrix& D)
    {
        const int size1 = (order1+1)*(order1+2)/2;
        const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(D.TMV_colsize()) >= size1);
        Assert(int(D.TMV_rowsize()) >= size2);
        // First I should point out an error in BJ02.  It lists
        // D(pq,00) = e^mu sech(mu) tanh(mu)^p delta_pq
        // This is wrong.
        // This is the formula for D(00,pq).
        //
        // Second, this formula has a different normalization than what we 
        // want.  This gives the transformation I(x,y) -> I(exp(-mu)x,exp(-mu)y).
        // However, we want the transformation for the same I, but measured in
        // a different coordinate system.  The difference is in the fact that
        // the flux scales as sigma^2, so this gives another factor of exp(-2mu).
        // So the formula we need is:
        // D(00,pq) = e^-mu sech(mu) tanh(mu)^p delta_pq
        //
        // Third, the recursion suggested for calculating D(st,pq) 
        // is a bit cumbersome since it has a q+1 on the rhs, so the 
        // calculation of D(30,30) for example, would require knowledge 
        // of D(00,33) which is a higher order than we really need.
        // 
        // An easier recursion is the following:
        // D(00,pq) = e^-mu sech(mu) tanh(mu)^p delta_pq
        // D(s0,pq) = 1/sqrt(s) sech(mu) sqrt(p) D(s-10,p-1q)
        // D(st,pq) = 1/sqrt(t) [ sech(mu) sqrt(q) D(st-1,pq-1)
        //                        - tanh(mu) sqrt(s) D(s-1t-1,pq) ]
        // Note: Only s-t = p-q terms are nonzero.
        //       This means that this can be significantly sped up
        //       by forming smaller matrices for each m, and using 
        //       permutations to rotate b so that these elements are
        //       continuous and then rotate back after multiplying.
        //       However, I'll wait to implement this until such speed 
        //       is found to be necessary.

        D.setZero();

        if (mu == 0.0) { 
            D.TMV_diagpart(0,0,std::min(size1,size2)).TMV_setAllTo(1.);
            return; 
        }

        double tmu = tanh(mu);
        double smu = 1./cosh(mu);

        int minorder = std::min(order1,order2);
        double Dqq = exp(-mu)*smu;
        // Dqq = exp(-mu) sech(mu) (-tanh(mu))^q
        // This variable keeps the latest Dqq value:
        for(int n=0,pq=0;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                double* Dpq = TMV_ptr(D.col(pq));
                double* Dpq_n = TMV_ptr(D.col(pq-n));
                double* Dpq_n_2 = q>0?TMV_ptr(D.col(pq-n-2)):0;
                double* Dpq1 = p>q?TMV_ptr(D.col(pq+1)):0;
                if (p==q) {
                    if (p > 0) Dqq *= tmu;
                    Dpq[0] = Dqq;  // D(0,pq)
                }
                for(int m=1,st=1;m<=minorder;++m) {
                    for(int s=m,t=0;s>=t;--s,++t,++st) {
                        if (p-q==s-t) {
                            double temp;
                            if (t == 0) {
                                temp = smu*sqrtn(p)*Dpq_n[st-m]/sqrtn(s);
                            } else {
                                temp = -tmu*sqrtn(s)*Dpq[st-2*m-1];
                                if (q > 0) {
                                    temp += smu*sqrtn(q)*Dpq_n_2[st-m-2];
                                }
                                temp /= sqrtn(t);
                            }
                            if (s == t) {
                                Dpq[st] = temp; // D(st,pq)
                            } else {
                                Dpq[st] = temp; // D(st,pq)
                                Dpq1[st+1] = temp; // D(st+1,pq+1)
                            }
                        }
                        if (s!=t) ++st;
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void augmentMuTransformCols(double mu, int order, DMatrix& D)
    {
        const int order1 = order;
        const int order2 = order+2;
        const int size1 = (order1+1)*(order1+2)/2;
        //const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(D.TMV_colsize()) >= size1);
        Assert(int(D.TMV_rowsize()) >= (order2+1)*(order2+2)/2);

        if (mu == 0.0) return;

        double tmu = tanh(mu);
        double smu = 1./cosh(mu);

        double Dqq = exp(-mu)*smu;
        for(int q=order1/2;q>0;--q) Dqq *= tmu;
        for(int n=order1+1,pq=size1;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                double* Dpq = TMV_ptr(D.col(pq));
                double* Dpq_n = TMV_ptr(D.col(pq-n));
                double* Dpq_n_2 = q>0?TMV_ptr(D.col(pq-n-2)):0;
                double* Dpq1 = p>q?TMV_ptr(D.col(pq+1)):0;
                if (p==q) {
                    if (p > 0) Dqq *= tmu;
                    Dpq[0] = Dqq;  // D(0,pq)
                }
                for(int m=1,st=1;m<=order1;++m) {
                    for(int s=m,t=0;s>=t;--s,++t,++st) {
                        if (p-q==s-t) {
                            double temp;
                            if (t == 0) {
                                temp = smu*sqrtn(p)*Dpq_n[st-m]/sqrtn(s);
                            } else {
                                temp = -tmu*sqrtn(s)*Dpq[st-2*m-1];
                                if (q > 0) {
                                    temp += smu*sqrtn(q)*Dpq_n_2[st-m-2];
                                }
                                temp /= sqrtn(t);
                            }
                            Dpq[st] = temp; 
                            if (s!=t) Dpq1[st+1] = temp; 
                        }
                        if (s!=t) ++st;
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void augmentMuTransformRows(double mu, int order, DMatrix& D)
    {
        const int order1 = order+2;
        const int order2 = order;
        //const int size1 = (order1+1)*(order1+2)/2;
        const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(D.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(D.TMV_rowsize()) == size2);

        if (mu == 0.0) return;

        double tmu = tanh(mu);
        double smu = 1./cosh(mu);

        for(int n=0,pq=0;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                double* Dpq = TMV_ptr(D.col(pq));
                double* Dpq_n = TMV_ptr(D.col(pq-n));
                double* Dpq_n_2 = q>0?TMV_ptr(D.col(pq-n-2)):0;
                double* Dpq1 = p>q?TMV_ptr(D.col(pq+1)):0;
                for(int m=order2+1,st=size2;m<=order1;++m) {
                    for(int s=m,t=0;s>=t;--s,++t,++st) {
                        if (p-q==s-t) {
                            double temp;
                            if (t == 0) {
                                temp = smu*sqrtn(p)*Dpq_n[st-m]/sqrtn(s);
                            } else {
                                temp = -tmu*sqrtn(s)*Dpq[st-2*m-1];
                                if (q > 0) {
                                    temp += smu*sqrtn(q)*Dpq_n_2[st-m-2];
                                }
                                temp /= sqrtn(t);
                            }
                            Dpq[st] = temp;
                            if (s!=t) Dpq1[st+1] = temp;
                        }
                        if (s!=t) ++st;
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void applyMu(double mu, BVec& b)
    {
        if (mu != 0.0) {
            DMatrix D(int(b.size()),int(b.size()));
            calculateMuTransform(mu,b.getOrder(),D);
            b.vec() = D * b.vec();
        }
    }

    void calculateThetaTransform(
        double theta, int order1, int order2, DBandMatrix& R)
    {
        const int size1 = (order1+1)*(order1+2)/2;
        const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(R.TMV_colsize()) >= size1);
        Assert(int(R.TMV_rowsize()) >= size2);

        R.setZero();

        if (theta == 0.0) { 
            R.TMV_diagpart(0,0,std::min(size1,size2)).TMV_setAllTo(1.);
            return; 
        }

        int minorder = std::min(order1,order2);
        std::vector<std::complex<double> > expimt(minorder+1);
        expimt[0] = 1.;
        if (minorder > 0) expimt[1] = std::polar(1.,theta);
        for(int m=2;m<=minorder;++m) expimt[m] = expimt[m-1] * expimt[1];

        for(int n=0,pq=0;n<=minorder;++n) {
            for(int p=n,q=0,m=n;p>=q;--p,++q,++pq,m-=2) {
                if (m==0) {
                    R(pq,pq) = 1.;
                } else {
                    R(pq,pq) = real(expimt[m]);
                    R(pq,pq+1) = -imag(expimt[m]);
                    R(pq+1,pq) = imag(expimt[m]);
                    R(pq+1,pq+1) = real(expimt[m]);
                    ++pq;
                }
            }
        }
    }

    void applyTheta(double theta, BVec& b)
    {
        if (theta != 0.0) {
#ifdef USE_TMV
            DBandMatrix R(int(b.size()),int(b.size()),1,1);
#else
            DBandMatrix R(int(b.size()),int(b.size()));
#endif
            calculateThetaTransform(theta,b.getOrder(),R);
            b.vec() = R * b.vec();
        }
    }

    void calculateGTransform(
        std::complex<double> g, int order1, int order2, DMatrix& S)
    {
        const int size1 = (order1+1)*(order1+2)/2;
        const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(S.TMV_colsize()) >= size1);
        Assert(int(S.TMV_rowsize()) >= size2);
        // S(st,pq) = f(p,s) f(q,t) (eta/|eta|)^(s-t-p+q)
        // f(p,0) = sqrt(p!)/(p/2)! sqrt(sech(|eta|/2)) (-tanh(|eta|/2)/2)^p/2
        // f(p,s+1) = (sqrt(p) sech(|eta|/2) f(p-1,s) +
        //                      sqrt(s) tanh(|eta|/2) f(p,s-1))/sqrt(s+1)
        //
        // tanh(|eta|/2) = |g|
        // sech(|eta|/2) = sqrt(1-|g|^2)
        //
        // Note: Like with the matrix in applyMu, this one is also fairly
        // sparse.  Could get a speedup by expoiting that, but currently don't.
        // I'll wait until the speedup is found to be necessary.

        S.setZero();

        if (g == 0.0) { 
            S.TMV_diagpart(0,0,std::min(size1,size2)).TMV_setAllTo(1.);
            return; 
        }

        double absg = std::abs(g);
        double normg = std::norm(g);
        std::vector<std::complex<double> > phase(order1+order2+1);
        phase[0] = 1.;
        std::complex<double> ph = -std::conj(g)/absg;
        // I'm not sure why conj was needed here.  Maybe there is an
        // error in the phase indices below that this corrects.
        for(int i=1;i<=order1+order2;++i) phase[i] = phase[i-1]*ph;

        double te = absg;
        double se = sqrt(1.-normg);

        std::vector<std::vector<double> > f(
            order2+1,std::vector<double>(order1+1,0.));

        f[0][0] = sqrt(se); // f(0,0)
        // only terms with p+s even are non-zero.
        for(int p=2;p<=order2;p+=2) {
            f[p][0] = f[p-2][0]*(-te/2.)*sqrtn(p*(p-1))/double(p/2); // f(p,0)
        }

        for(int s=0;s<order1;++s) {
            if (s%2==1) {
                f[0][s+1] = sqrtn(s)*te*f[0][s-1]/sqrtn(s+1); // f(0,s+1)
            }
            for(int p=s%2+1;p<=order2;p+=2) {
                double temp = sqrtn(p)*se*f[p-1][s]; 
                if (s>0) temp += sqrtn(s)*te*f[p][s-1];
                temp /= sqrtn(s+1); 
                f[p][s+1] = temp; // f(p,s+1)
            }
        }

        for(int n=0,pq=0;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                const std::vector<double>& fp = f[p];
                const std::vector<double>& fq = f[q];
                double* Spq = TMV_ptr(S.col(pq));
                double* Spq1 = p>q?TMV_ptr(S.col(pq+1)):0;
                for(int nn=n%2,st=(nn==0?0:1);nn<=order1;nn+=2,st+=nn) {
                    for(int s=nn,t=0;s>=t;--s,++t,++st) {

                        double s0 = fp[s] * fq[t];

                        int iphase = s-t-p+q;
                        std::complex<double> s1 = s0 * 
                            (iphase >= 0 ?
                             phase[iphase/2] :
                             std::conj(phase[-iphase/2]));

                        if (p==q) {
                            if (s==t) {
                                Spq[st] = std::real(s1);
                            } else {
                                Spq[st] = std::real(s1);
                                Spq[st+1] = std::imag(s1);
                                ++st;
                            }
                        } else if (s==t) {
                            // b_st = t_stpq b_pq + t_stqp b_qp
                            // In this case with s=t, t_stpq == t_stqp*
                            Spq[st] = 2.*std::real(s1);
                            Spq1[st] = -2.*std::imag(s1);
                        } else {
                            s0 = fq[s] * fp[t];
                            iphase = s-t-q+p;
                            std::complex<double> s2 = s0 *
                                (iphase >= 0 ?
                                 phase[iphase/2] :
                                 std::conj(phase[-iphase/2]));
                            Spq[st] = std::real(s1) + std::real(s2);
                            Spq[st+1] = std::imag(s1) + std::imag(s2);
                            Spq1[st] = -std::imag(s1) + std::imag(s2);
                            Spq1[st+1] = std::real(s1) - std::real(s2);
                            ++st;
                        }
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void augmentGTransformCols(std::complex<double> g, int order, DMatrix& S)
    {
        const int order1 = order;
        const int order2 = order+2;
        const int size1 = (order1+1)*(order1+2)/2;
        //const int size2 = (order2+1)*(order2+2)/2;
        Assert(int(S.TMV_colsize()) == size1);
        Assert(int(S.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        if (g == 0.0) return; 

        double absg = std::abs(g);
        double normg = std::norm(g);
        std::vector<std::complex<double> > phase(order1+order2+1);
        phase[0] = 1.;
        std::complex<double> ph = -std::conj(g)/absg;
        for(int i=1;i<=order1+order2;++i) phase[i] = phase[i-1]*ph;

        double te = absg;
        double se = sqrt(1.-normg);

        std::vector<std::vector<double> > f(
            order2+1,std::vector<double>(order1+1,0.));

        f[0][0] = sqrt(se); // f(0,0)
        // only terms with p+s even are non-zero.
        for(int p=2;p<=order2;p+=2) {
            f[p][0] = f[p-2][0]*(-te/2.)*sqrtn(p*(p-1))/double(p/2); // f(p,0)
        }

        for(int s=0;s<order1;++s) {
            if (s%2==1) {
                f[0][s+1] = sqrtn(s)*te*f[0][s-1]/sqrtn(s+1); // f(0,s+1)
            }
            for(int p=s%2+1;p<=order2;p+=2) {
                double temp = sqrtn(p)*se*f[p-1][s]; 
                if (s>0) temp += sqrtn(s)*te*f[p][s-1];
                temp /= sqrtn(s+1); 
                f[p][s+1] = temp; // f(p,s+1)
            }
        }

        for(int n=order1+1,pq=size1;n<=order2;++n) {
            for(int p=n,q=0;p>=q;--p,++q,++pq) {
                const std::vector<double>& fp = f[p];
                const std::vector<double>& fq = f[q];
                double* Spq = TMV_ptr(S.col(pq));
                double* Spq1 = p>q?TMV_ptr(S.col(pq+1)):0;
                for(int nn=n%2,st=(nn==0?0:1);nn<=order1;nn+=2,st+=nn) {
                    for(int s=nn,t=0;s>=t;--s,++t,++st) {

                        double s0 = fp[s] * fq[t];
                        int iphase = s-t-p+q;
                        std::complex<double> s1 = s0 *
                            (iphase >= 0 ?
                             phase[iphase/2] :
                             std::conj(phase[-iphase/2]));

                        if (p==q) {
                            if (s==t) {
                                Spq[st] = std::real(s1);
                            } else {
                                Spq[st] = std::real(s1);
                                Spq[st+1] = std::imag(s1);
                                ++st;
                            }
                        } else if (s==t) {
                            Spq[st] = 2.*std::real(s1);
                            Spq1[st] = -2.*std::imag(s1);
                        } else {
                            s0 = fq[s] * fp[t];
                            iphase = s-t-q+p;
                            std::complex<double> s2 = s0 *
                                (iphase >= 0 ?
                                 phase[iphase/2] :
                                 std::conj(phase[-iphase/2]));
                            Spq[st]= std::real(s1) + std::real(s2);
                            Spq1[st] = -std::imag(s1) + std::imag(s2);
                            Spq[st+1] = std::imag(s1) + std::imag(s2);
                            Spq1[st+1] = std::real(s1) - std::real(s2);
                            ++st;
                        }
                    }
                }
                if (p!=q) ++pq;
            }
        }
    }

    void applyG(std::complex<double> g, BVec& b)
    {
        if (g != 0.0) {
            DMatrix S(int(b.size()),int(b.size()));
            calculateGTransform(g,b.getOrder(),S);
            b.vec() = S * b.vec();
        }
    }

    void calculatePsfConvolve(
        const BVec& bpsf, int order1, int order2, double sigma, DMatrix& C)
    {
        //xdbg<<"Start calculatePsfConvolve\n";
        //xdbg<<"bpsf = "<<bpsf<<std::endl;
        //xdbg<<"order = "<<order<<std::endl;
        //xdbg<<"sigma = "<<sigma<<std::endl;
        // Just make the matrix which multiplies b to effect the convolution.
        // b_obs = C * b_init
        // However, we do not actually use it this way.  Rather, we use it to 
        // switch the model from:
        // I = Sum psi_pq b_obs_pq
        // to 
        // I = Sum psi_pq C b_init_pq
        // We use this to solve for the ML b_init.
        const int order3 = bpsf.getOrder();
        Assert(int(C.TMV_colsize()) >= (order1+1)*(order2+2)/2);
        Assert(int(C.TMV_rowsize()) >= (order2+1)*(order2+2)/2);
        Assert(int(bpsf.size()) == (order3+1)*(order3+2)/2);

        // C(st,pq) = 1/sqrt(pi) Sum_uv sqrt(p!u!q!v!/s!t!)/w! 
        //                               G(s,p,u) G(t,q,v) bpsf_uv
        // The sum is only over terms for which  p+u-s == q+v-t,
        // and w = p+u-s = q+v-t >= 0
        // 
        // G(0,p,u) = binom(p+u,u) (-A)^u B^p
        // G(s+1,p,u) = A G(s,p-1,u) + B G(s,p,u-1)
        // where A = sigma_init / sigma_obs, B = sigma_psf / sigma_obs
        // D = A^2 = 1-B^2
        //
        // It is more efficient to combine the sqrt(p!u!/s!w!) into the G(s,p,u).
        // Call the product H(s,p,u).  We need to translate the above formulae:
        // H(0,p,u) = sqrt(p!u!/(p+u)!) (p+u)!/(p!u!) (-A)^u B^p
        //          = sqrt((p+u)!/(p!u!)) (-A)^u B^p
        // H(s+1,p,u) = sqrt(p!u!/(s+1)!(p+u-s-1)!) * [
        //                    A H(s,p-1,u) / sqrt((p-1)!u!/s!(p+u-s-1)!) +
        //                    B H(s,p,u-1) / sqrt(p!(u-1)!/s!(p+u-s-1)!) ]
        //            = A sqrt(p)/sqrt(s+1) H(s,p-1,u) + 
        //              B sqrt(u)/sqrt(s+1) H(s,p,u-1)

        C.setZero();
        
        // sigma^2 = exp(mu)
        // D = sigma_i^2 / (sigma_i^2 + sigma_psf^2) 
        //   = 1 / (1 + (sigma_psf/sigma_i)^2)
        double D = 1./(1.+pow(bpsf.getSigma()/sigma,2));
        double A = sqrt(D);
        double B = sqrt(1-D);

        std::vector<std::vector<std::vector<double> > > H(
            order1+1,std::vector<std::vector<double> >(
                order2+1,std::vector<double>(order3+1)));

        H[0][0][0] = 1.; // H[0](0,0)
        for(int u=0;u<=order3;++u) {
            if (u>0) H[0][0][u] = -A * H[0][0][u-1];
            for(int p=1;p<=order2;++p) 
                H[0][p][u] = B*sqrtn(p+u)/sqrtn(p)*H[0][p-1][u]; // H[0](p,u)
        }
        for(int s=0;s<order1;++s) {
            H[s+1][0][0] = 0.;
            for(int p=1;p<=order2;++p) 
                H[s+1][p][0] = A*sqrtn(p)*H[s][p-1][0]/sqrtn(s+1);
            for(int u=1;u<=order3;++u) {
                H[s+1][0][u] = B*sqrtn(u)*H[s][0][u-1]/sqrtn(s+1);
                for(int p=1;p<=order2;++p) 
                    H[s+1][p][u] = (A*sqrtn(p)*H[s][p-1][u] + 
                                    B*sqrtn(u)*H[s][p][u-1])/ sqrtn(s+1);
            }
        }

        int pq = 0;
        const double* bpsfv = TMV_cptr(bpsf.vec());
        for(int n=0;n<=order2;++n) {
            for(int p=n,q=0;p>=q;(p==q?++pq:pq+=2),--p,++q) {
                //xdbg<<"Start column "<<pq<<" = "<<p<<','<<q<<std::endl;
                int pmq = p-q;
                double* Cpq = TMV_ptr(C.col(pq));
                double* Cpq1 = p>q?TMV_ptr(C.col(pq+1)):0;
                int st = 0;
                for(int nn=0;nn<=order1;++nn) {
                    for(int s=nn,t=0;s>=t;(s==t?++st:st+=2),--s,++t) {
                        //xdbg<<"st = "<<st<<" = "<<s<<','<<t<<std::endl;
                        int smt = s-t;
                        double Cpqst = 0.;
                        double Cpq1st = 0.;
                        double Cpqst1 = 0.;
                        double Cpq1st1 = 0.;
                        const std::vector<double>& Hsp = H[s][p];
                        const std::vector<double>& Hsq = H[s][q];
                        const std::vector<double>& Htp = H[t][p];
                        const std::vector<double>& Htq = H[t][q];
                        int uv0 = 0;
                        int parity = (n+nn)%2;
                        for(int upv=0;upv<=order3;++upv,uv0+=upv) {
                            if (upv % 2 != parity) continue;
                            // There are three values of u-v that are worth considering:
                            // u-v = (s-t) - (p-q) >= 0
                            // u-v = (s-t) + (p-q) > 0
                            // u-v = (p-q) - (s-t) < 0

                            int umv = smt-pmq;
                            if (umv >= 0 && umv <= upv) {
                                // First do terms with p>=q, u>=v  (always keep s>=t)
                                // s-t = p-q + u-v
                                int u = (upv+umv)/2;
                                int v = (upv-umv)/2;

                                int w = p+u-s;
                                Assert(q+v-t == w);
                                //Assert((w >= 0) == (umv >= 0));
                                if (w >= 0) {
                                    int uv = uv0 + 2*v;
                                    if (umv == 0) {
                                        double temp = Hsp[u]*Htq[v]*bpsfv[uv];
                                        if (s==t) {
                                            Assert(p==q);
                                            Cpqst += temp;
                                        } else {
                                            Cpqst += temp;
                                            Cpq1st1 += temp;
                                        }
                                    } else {
                                        Assert(s>t);
                                        double tempr = Hsp[u]*Htq[v];
                                        double tempi = tempr * bpsfv[uv+1];
                                        tempr *= bpsfv[uv];
                                        if (p==q) {
                                            Cpqst += tempr;
                                            Cpqst1 += tempi;
                                        } else {
                                            Assert(p>q);
                                            Cpqst += tempr;
                                            Cpqst1 += tempi;
                                            Cpq1st -= tempi;
                                            Cpq1st1 += tempr;
                                        }
                                    }
                                }
                            }

                            umv = smt+pmq;
                            if (pmq != 0 && umv <= upv) {
                                // Next p<q, u>v.  Implement by swapping p,q
                                // These terms account for the fact that 
                                // b_init_qp = b_init_pq*
                                // s-t = q-p + u-v

                                Assert(umv > 0);
                                int u = (upv+umv)/2;
                                int v = (upv-umv)/2;

                                int w = q+u-s;
                                Assert(p+v-t == w);
                                //Assert((w >= 0) == (umv > 0));
                                if (w >= 0) {
                                    int uv = uv0 + 2*v;
                                    //Assert(w > 0);
                                    //Assert((w > 0) == (umv > 0));
                                    Assert(u>v);
                                    double tempr = Hsq[u]*Htp[v];
                                    double tempi = tempr * bpsfv[uv+1];
                                    tempr *= bpsfv[uv];
                                    if (smt==0) {
                                        Cpqst += tempr;
                                        Cpq1st += tempi;
                                    } else {
                                        Cpqst += tempr;
                                        Cpqst1 += tempi;
                                        Cpq1st += tempi;
                                        Cpq1st1 -= tempr;
                                    }
                                }
                            }

                            umv = pmq-smt;
                            if (umv > 0 && umv <= upv) {
                                // Next p>q, u<v.
                                // These terms account for b_psf_vu = b_psf_uv*
                                int u = (upv+umv)/2;
                                int v = (upv-umv)/2;

                                int w = p+v-s;
                                Assert(q+u-t == w);
                                if (w >= 0) {
                                    int uv = uv0 + 2*v;
                                    // s-t = p-q + v-u
                                    Assert(p>q);
                                    double tempr = Hsp[v]*Htq[u];
                                    double tempi = -tempr * bpsfv[uv+1];
                                    tempr *= bpsfv[uv];
                                    if (smt==0) {
                                        Cpqst += tempr;
                                        Cpq1st -= tempi;
                                    } else {
                                        Cpqst += tempr;
                                        Cpqst1 += tempi;
                                        Cpq1st -= tempi;
                                        Cpq1st1 += tempr;
                                    }
                                }
                            }
                        }
                        Assert(uv0 == int(bpsf.size()));
                        //xdbg<<"Cpq[st] = "<<Cpqst<<std::endl;
                        Cpq[st] = Cpqst;
                        if (smt != 0) Cpq[st+1] = Cpqst1;
                        if (pmq != 0) {
                            Cpq1[st] = Cpq1st;
                            if (smt != 0) Cpq1[st+1] = Cpq1st1;
                        }
                    }
                }
                Assert(st == int(C.TMV_colsize()));
            }
        }
        Assert(pq == int(C.TMV_rowsize()));
        C /= afwGeom::SQRTPI;
    }

    void applyPsf(const BVec& bpsf, BVec& b)
    {
        DMatrix C(int(b.size()),int(b.size()));
        calculatePsfConvolve(bpsf,b.getOrder(),b.getSigma(),C);
        b.vec() = C * b.vec();
    }

}}}}
