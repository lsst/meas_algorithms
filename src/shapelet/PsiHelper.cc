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

#include "lsst/meas/algorithms/shapelet/PsiHelper.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include "lsst/afw/geom/Angle.h"

namespace afwGeom = lsst::afw::geom;

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    // Here is the order of p,q along the indices of psi:
    //
    // k 0 1,2 3,4 5 6,7 8,9 10,11 12,13 14 15,16 17,18 19,20 21,22 23,24 25,26 27
    // p 0  1   2  1  3   2    4     3    2   5     4     3     6     5     4    3
    // q 0  0   0  1  0   1    0     1    2   0     1     2     0     1     2    3
    // n 0  1   2  2  3   3    4     4    4   5     5     5     6     6     6    6
    // m 0  1   2  0  3   1    4     2    0   5     3     1     6     4     2    0
    //

#ifdef USE_TMV
    // The Eigen version required pretty significant rewriting before it was
    // remotely efficient, so I separate this out with an ifdef, rather than
    // use the TMV and EIGEN macros.

    void makePsi(DMatrix& psi, CDVectorView z, int order, const DVector* coeff)
    {
        // For p>=q:
        //
        // psi_pq = (pi p! q!)^-1/2 z^m exp(-r^2/2) K_pq(r^2)
        //
        // where K_pq(r^2) satisfies the recursion relation:
        // K_pq = (r^2 - (N-1)) K_p-1,q-1 - (p-1)(q-1) K_p-2,q-2
        // with K_p0 = 1
        //
        // The recursion for psi_pq can then be derived as:
        //
        // psi_pq = (pq)^-1/2 (r^2-(N-1)) psi_p-1,q-1
        //          - sqrt( (p-1)(q-1)/(pq) ) psi_p-2,q-2
        //
        // psi_p0 = (z/sqrt(p)) psi_p-1,0
        // 
        // psi_00 = 1/sqrt(pi) exp(-r^2/2)

        Assert(int(psi.rowsize()) >= (order+1)*(order+2)/2);
        Assert(psi.colsize() == z.size());
        if (coeff) Assert(psi.colsize() == coeff->size());
        Assert(psi.iscm());
        Assert(!psi.isconj());

        // Setup rsq, z vectors and set psi_00
        DVector rsq(z.size());
        double* rsqit = rsq.ptr();
        double* psi00it = psi.ptr();
        const std::complex<double>* zit = z.cptr();
        const int zsize = z.size();
        for(int i=0;i<zsize;++i) {
            rsqit[i] = std::norm(zit[i]);
            psi00it[i] = afwGeom::INVSQRTPI * exp(-(rsqit[i])/2.);
        }
        if (coeff) psi.col(0) *= DiagMatrixViewOf(*coeff);

        DVector zr = z.realPart();
        DVector zi = z.imagPart();
        if (order >= 1) {
            // Set psi_10
            // All m > 0 elements are intrinsically complex.
            // However, we are fitting to a real intensity pattern
            // with complex coefficients of the complex shapelets.
            // Since psi_pq = psi_qp* (* = complex conjugate),
            // we know that b_pq must be b_qp*
            // b_pq psi_pq + b_pq* psi_pq* = 2 b_pqr psi_pqr - 2 b_pqi psi_pqi
            // So the values we want for the real fitter are
            // really 2 real(psi_pq) and -2 imag(psi_pq)
            // Putting the 2's here carries through to the rest of the 
            // elements via the recursion.
            psi.col(1) = 2. * DiagMatrixViewOf(zr) * psi.col(0);
            psi.col(2) = -2. * DiagMatrixViewOf(zi) * psi.col(0);
        }
        for(int N=2,k=3;N<=order;++N) {
            // Set psi_N0
            // The signs of these are not what you naively think due to 
            // the +2, -2 discussed above.  You just have to follow through
            // what the complex psi_N0 is, and what value is stored in the
            // psi_N-1,0 location, and what needs to get stored here.
            psi.col(k) = sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N);
            psi.col(k) += sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N+1);
            psi.col(k+1) = -sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N);
            psi.col(k+1) += sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N+1);
            k+=2;

            // Set psi_pq with q>0
            // The rsq part of this calculation can be done in batch, which 
            // speeds things up a bit.
            psi.colRange(k,k+N-1) =
                DiagMatrixViewOf(rsq) * psi.colRange(k-2*N-1,k-N-2);
            psi.colRange(k,k+N-1) -= (N-1.) * psi.colRange(k-2*N-1,k-N-2);
            // The other calculation steps are different for each component:
            for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
                double pq = p*q;
                if (m==0) {
                    psi.col(k) /= sqrt(pq);
                    if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
                    ++k;
                } else {
                    psi.colRange(k,k+2) /= sqrt(pq);
                    if (q > 1)
                        psi.colRange(k,k+2) -= 
                            sqrt(1.-(N-1.)/pq)*psi.colRange(k+2-4*N,k+4-4*N);
                    k+=2;
                }
            }
        }
    }

    void makePsi(DVector& psi, std::complex<double> z, int order)
    {
        double rsq = std::norm(z);
        psi(0) = afwGeom::INVSQRTPI * exp(-rsq/2.);

        double zr = std::real(z);
        double zi = std::imag(z);
        if (order >= 1) {
            psi(1) = 2. * zr * psi(0);
            psi(2) = -2. * zi * psi(0);
        }
        for(int N=2,k=3;N<=order;++N) {
            psi(k) = sqrt(1./N) * zr * psi(k-N);
            psi(k) += sqrt(1./N) * zi * psi(k-N+1);
            psi(k+1) = -sqrt(1./N) * zi * psi(k-N);
            psi(k+1) += sqrt(1./N) * zr * psi(k-N+1);
            k+=2;

            psi.subVector(k,k+N-1) = rsq * psi.subVector(k-2*N-1,k-N-2);
            psi.subVector(k,k+N-1) -= (N-1.) * psi.subVector(k-2*N-1,k-N-2);
            for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
                double pq = p*q;
                if (m==0) {
                    psi(k) /= sqrt(pq);
                    if (q > 1) psi(k) -= sqrt(1.-(N-1.)/pq)*psi(k+2-4*N);
                    ++k;
                } else {
                    psi.subVector(k,k+2) /= sqrt(pq);
                    if (q > 1)
                        psi.subVector(k,k+2) -= 
                            sqrt(1.-(N-1.)/pq)*psi.subVector(k+2-4*N,k+4-4*N);
                    k+=2;
                }
            }
        }
    }

    void augmentPsi(DMatrix& psi, CDVectorView z, int order)
    {
        Assert(int(psi.rowsize()) >= (order+3)*(order+4)/2);
        Assert(psi.colsize() == z.size());
        Assert(order >= 1);
        //Assert(psi.iscm());
        //Assert(!psi.isconj());

        DVector rsq(z.size());
        double* rsqit = rsq.ptr();
        const std::complex<double>* zit = z.cptr();
        const int zsize = z.size();
        for(int i=0;i<zsize;++i) {
            rsqit[i] = std::norm(zit[i]);
        }

        DVector zr = z.realPart();
        DVector zi = z.imagPart();
        for(int N=order+1,k=N*(N+1)/2;N<=order+2;++N) {
            psi.col(k) = sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N);
            psi.col(k) += sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N+1);
            psi.col(k+1) = -sqrt(1./N) * DiagMatrixViewOf(zi) * psi.col(k-N);
            psi.col(k+1) += sqrt(1./N) * DiagMatrixViewOf(zr) * psi.col(k-N+1);
            k+=2;

            psi.colRange(k,k+N-1) =
                DiagMatrixViewOf(rsq) * psi.colRange(k-2*N-1,k-N-2);
            psi.colRange(k,k+N-1) -= (N-1.) * psi.colRange(k-2*N-1,k-N-2);
            for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
                double pq = p*q;
                if (m==0) {
                    psi.col(k) /= sqrt(pq);
                    if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
                    ++k;
                } else {
                    psi.colRange(k,k+2) /= sqrt(pq);
                    if (q > 1)
                        psi.colRange(k,k+2) -= 
                            sqrt(1.-(N-1.)/pq)*psi.colRange(k+2-4*N,k+4-4*N);
                    k+=2;
                }
            }
        }
    }

#else

    void makePsi(DMatrix& psi, CDVectorView z, int order, const DVector* coeff)
    {
        // For p>=q:
        //
        // psi_pq = (pi p! q!)^-1/2 z^m exp(-r^2/2) K_pq(r^2)
        //
        // where K_pq(r^2) satisfies the recursion relation:
        // K_pq = (r^2 - (N-1)) K_p-1,q-1 - (p-1)(q-1) K_p-2,q-2
        // with K_p0 = 1
        //
        // The recursion for psi_pq can then be derived as:
        //
        // psi_pq = (pq)^-1/2 (r^2-(N-1)) psi_p-1,q-1
        //          - sqrt( (p-1)(q-1)/(pq) ) psi_p-2,q-2
        //
        // psi_p0 = (z/sqrt(p)) psi_p-1,0
        // 
        // psi_00 = 1/sqrt(pi) exp(-r^2/2)

        Assert(int(psi.TMV_rowsize()) >= (order+1)*(order+2)/2);
        Assert(psi.TMV_colsize() == z.size());
        if (coeff) Assert(psi.TMV_colsize() == coeff->size());
        //Assert(psi.iscm());
        //Assert(!psi.isconj());

        // Setup rsq, z vectors and set psi_00
        DVector rsq(z.size());
        double* rsqit = TMV_ptr(rsq);
        double* psi00it = TMV_ptr(psi);
        const std::complex<double>* zit = TMV_cptr(z);
        const int zsize = z.size();
        for(int i=0;i<zsize;++i) {
            rsqit[i] = std::norm(zit[i]);
            psi00it[i] = afwGeom::INVSQRTPI * exp(-(rsqit[i])/2.);
        }
        if (coeff) psi.col(0).array() = coeff->array() * psi.col(0).array();

        DVector zr = z.TMV_realPart();
        DVector zi = z.TMV_imagPart();
        if (order >= 1) {
            // Set psi_10
            // All m > 0 elements are intrinsically complex.
            // However, we are fitting to a real intensity pattern
            // with complex coefficients of the complex shapelets.
            // Since psi_pq = psi_qp* (* = complex conjugate),
            // we know that b_pq must be b_qp*
            // b_pq psi_pq + b_pq* psi_pq* = 2 b_pqr psi_pqr - 2 b_pqi psi_pqi
            // So the values we want for the real fitter are
            // really 2 real(psi_pq) and -2 imag(psi_pq)
            // Putting the 2's here carries through to the rest of the 
            // elements via the recursion.
            psi.col(1).array() = zr.array() * psi.col(0).array();
            psi.col(2).array() = (-zi).array() * psi.col(0).array();
            psi.col(1) *= 2.;
            psi.col(2) *= 2.;
        }
        for(int N=2,k=3;N<=order;++N) {
            // Set psi_N0
            // The signs of these are not what you naively think due to 
            // the +2, -2 discussed above.  You just have to follow through
            // what the complex psi_N0 is, and what value is stored in the
            // psi_N-1,0 location, and what needs to get stored here.
            psi.col(k).array() = zr.array() * psi.col(k-N).array();
            psi.col(k).array() += zi.array() * psi.col(k-N+1).array();
            psi.col(k+1).array() = zr.array() * psi.col(k-N+1).array();
            psi.col(k+1).array() -= zi.array() * psi.col(k-N).array();
            double sqrt_1_N = sqrt(1./N);
            psi.col(k) *= sqrt_1_N;
            psi.col(k+1) *= sqrt_1_N;
            k+=2;

            // Set psi_pq with q>0
            // The rsq part of this calculation can be done in batch, which 
            // speeds things up a bit.
            TMV_colRange(psi,k,k+N-1) = rsq.asDiagonal() * TMV_colRange(psi,k-2*N-1,k-N-2);
            TMV_colRange(psi,k,k+N-1) -= (N-1.) * TMV_colRange(psi,k-2*N-1,k-N-2);
            // The other calculation steps are different for each component:
            for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
                double pq = p*q;
                if (m==0) {
                    psi.col(k) /= sqrt(pq);
                    if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
                    ++k;
                } else {
                    TMV_colRange(psi,k,k+2) /= sqrt(pq);
                    if (q > 1)
                        TMV_colRange(psi,k,k+2) -= 
                            sqrt(1.-(N-1.)/pq)*TMV_colRange(psi,k+2-4*N,k+4-4*N);
                    k+=2;
                }
            }
        }
    }

    void augmentPsi(DMatrix& psi, CDVectorView z, int order)
    {
        Assert(int(psi.TMV_rowsize()) >= (order+3)*(order+4)/2);
        Assert(psi.TMV_colsize() == z.size());
        Assert(order >= 1);
        //Assert(psi.iscm());
        //Assert(!psi.isconj());

        DVector rsq(z.size());
        double* rsqit = TMV_ptr(rsq);
        const std::complex<double>* zit = TMV_cptr(z);
        const int zsize = z.size();
        for(int i=0;i<zsize;++i) {
            rsqit[i] = std::norm(zit[i]);
        }

        DVector zr = z.TMV_realPart();
        DVector zi = z.TMV_imagPart();
        for(int N=order+1,k=N*(N+1)/2;N<=order+2;++N) {
            psi.col(k).array() = zr.array() * psi.col(k-N).array();
            psi.col(k).array() += zi.array() * psi.col(k-N+1).array();
            psi.col(k+1).array() = zr.array() * psi.col(k-N+1).array();
            psi.col(k+1).array() -= zi.array() * psi.col(k-N).array();
            double sqrt_1_N = sqrt(1./N);
            psi.col(k) *= sqrt_1_N;
            psi.col(k+1) *= sqrt_1_N;
            k+=2;

            TMV_colRange(psi,k,k+N-1) = rsq.asDiagonal() * TMV_colRange(psi,k-2*N-1,k-N-2);
            TMV_colRange(psi,k,k+N-1) -=
                (N-1.) * TMV_colRange(psi,k-2*N-1,k-N-2);
            for(int m=N-2,p=N-1,q=1;m>=0;--p,++q,m-=2) {
                double pq = p*q;
                if (m==0) {
                    psi.col(k) /= sqrt(pq);
                    if (q > 1) psi.col(k) -= sqrt(1.-(N-1.)/pq)*psi.col(k+2-4*N);
                    ++k;
                } else {
                    TMV_colRange(psi,k,k+2) /= sqrt(pq);
                    if (q > 1)
                        TMV_colRange(psi,k,k+2) -= 
                            sqrt(1.-(N-1.)/pq)*TMV_colRange(psi,k+2-4*N,k+4-4*N);
                    k+=2;
                }
            }
        }
    }

#endif

    void setupGx(DMatrix& Gx, int order1, int order2)
    {
        Assert(int(Gx.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(Gx.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        Gx.setZero();
        for(int n=0,k=0;n<=order2;++n) for(int p=n,q=0;p>=q;--p,++q,++k) {
            double dp = p;
            double dq = q;
            // d(psi(x,y))/dx = psi(x,y) Gx
            // d/dx = 1/2(ap + aq - apt - aqt)
            // Gx( p-1,q , pq ) = 1/2 sqrt(p)
            // Gx( p,q-1 , pq ) = 1/2 sqrt(q)
            // Gx( p+1,q , pq ) = -1/2 sqrt(p+1)
            // Gx( p,q+1 , pq ) = -1/2 sqrt(q+1)
            if (p==q) {
                if (q>0 && n<=order1+1) Gx(k-n-2,k) = sqrt(dp)/2.;
                if (n+1<=order1) Gx(k+n+1,k) = -sqrt(dp+1.)/2.;
            } else if (p==q+1) {
                if (n<=order1+1) Gx(k-n,k) = sqrt(dp);
                if (q>0 && n<=order1+1) Gx(k-n-2,k) = sqrt(dq)/2.;
                if (n+1<=order1) Gx(k+n+1,k) = -sqrt(dp+1.)/2.;
                if (n+1<=order1) Gx(k+n+3,k) = -sqrt(dq+1.);
                if (q>0 && n<=order1+1) Gx(k-n-1,k+1) = sqrt(dq)/2.;
                if (n+1<=order1) Gx(k+n+2,k+1) = -sqrt(dp+1.)/2.;
            } else {
                if (n<=order1+1) Gx(k-n,k) = sqrt(dp)/2.;
                if (q>0 && n<=order1+1) Gx(k-n-2,k) = sqrt(dq)/2.;
                if (n+1<=order1) Gx(k+n+1,k) = -sqrt(dp+1.)/2.;
                if (n+1<=order1) Gx(k+n+3,k) = -sqrt(dq+1.)/2.;
                if (n<=order1+1) Gx(k-n+1,k+1) = sqrt(dp)/2.;
                if (q>0 && n<=order1+1) Gx(k-n-1,k+1) = sqrt(dq)/2.;
                if (n+1<=order1) Gx(k+n+2,k+1) = -sqrt(dp+1.)/2.;
                if (n+1<=order1) Gx(k+n+4,k+1) = -sqrt(dq+1.)/2.;
            }
            if (p > q) ++k;
        }
    }

    void setupGy(DMatrix& Gy, int order1, int order2)
    {
        Assert(int(Gy.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(Gy.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        Gy.setZero();
        for(int n=0,k=0;n<=order2;++n) for(int p=n,q=0;p>=q;--p,++q,++k) {
            double dp = p;
            double dq = q;
            // d(psi(x,y))/dx = psi(x,y) Gx
            // d/dy = 1/2 i (ap - aq + apt - aqt)
            // Gy( p-1,q , pq ) = 1/2 i sqrt(p)
            // Gy( p,q-1 , pq ) = -1/2 i sqrt(q)
            // Gy( p+1,q , pq ) = 1/2 i sqrt(p+1)
            // Gy( p,q+1 , pq ) = -1/2 i sqrt(q+1)
            if (p==q) {
                if (q>0 && n<=order1+1) Gy(k-n-1,k) = -sqrt(dp)/2.;
                if (n+1<=order1) Gy(k+n+2,k) = sqrt(dp+1.)/2.;
            } else if (p==q+1) {
                if (q>0 && n<=order1+1) Gy(k-n-1,k) = -sqrt(dq)/2.;
                if (n+1<=order1) Gy(k+n+2,k) = sqrt(dp+1.)/2.;
                if (n<=order1+1) Gy(k-n,k+1) = -sqrt(dp);
                if (q>0 && n<=order1+1) Gy(k-n-2,k+1) = sqrt(dq)/2.;
                if (n+1<=order1) Gy(k+n+1,k+1) = -sqrt(dp+1.)/2.;
                if (n+1<=order1) Gy(k+n+3,k+1) = sqrt(dq+1.);
            } else {
                if (n+1<=order1) Gy(k-n+1,k) = sqrt(dp)/2.;
                if (q>0 && n<=order1+1) Gy(k-n-1,k) = -sqrt(dq)/2.;
                if (n+1<=order1) Gy(k+n+2,k) = sqrt(dp+1.)/2.;
                if (n+1<=order1) Gy(k+n+4,k) = -sqrt(dq+1.)/2.;
                if (n<=order1+1) Gy(k-n,k+1) = -sqrt(dp)/2.;
                if (q>0 && n<=order1+1) Gy(k-n-2,k+1) = sqrt(dq)/2.;
                if (n+1<=order1) Gy(k+n+1,k+1) = -sqrt(dp+1.)/2.;
                if (n+1<=order1) Gy(k+n+3,k+1) = sqrt(dq+1.)/2.;
            }

            if (p > q) ++k;
        }
    }

    void setupGg1(DMatrix& Gg1, int order1, int order2)
    {
        Assert(int(Gg1.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(Gg1.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        Gg1.setZero();
        for(int n=0,k=0;n<=order2;++n) for(int p=n,q=0;p>=q;--p,++q,++k) {
            double dp = p;
            double dq = q;
            // d(psi(x,y))/dg1 = psi(x,y) Gg1
            // Gg1 is G_g1
            // d/dg1 = x d/dx - y d/dy
            //       = 1/2 (ap^2 + aq^2 - apt^2 - aqt^2)
            // Gg1( p-2,q , pq ) = 1/2 sqrt(p(p-1))
            // Gg1( p,q-2 , pq ) = 1/2 sqrt(q(q-1))
            // Gg1( p+2,q , pq ) = -1/2 sqrt((p+1)(p+2))
            // Gg1( p,q+2 , pq ) = -1/2 sqrt((q+1)(q+2))
            if (p==q) {
                if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dp*(dp-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
            } else if (p==q+1) {
                if (p>1 && n<=order1+2) Gg1(k-2*n-1,k) = sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg1(k+2*n+5,k) = -sqrt((dq+1.)*(dq+2.))/2.;
                if (p>1 && n<=order1+2) Gg1(k-2*n,k+1) = -sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg1(k-2*n-2,k+1) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+4,k+1) = -sqrt((dp+1)*(dp+2.))/2.;
                if (n+2<=order1) Gg1(k+2*n+6,k+1) = sqrt((dq+1)*(dq+2.))/2.;
            } else if (p==q+2) {
                if (n<=order1+2) Gg1(k-2*n+1,k) = sqrt(dp*(dp-1.));
                if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg1(k+2*n+7,k) = -sqrt((dq+1.)*(dq+2.));
                if (q>1 && n<=order1+2) Gg1(k-2*n-2,k+1) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+4,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
            } else {
                if (n<=order1+2) Gg1(k-2*n+1,k) = sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg1(k-2*n-3,k) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+3,k) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg1(k+2*n+7,k) = -sqrt((dq+1.)*(dq+2.))/2.;
                if (n<=order1+2) Gg1(k-2*n+2,k+1) = sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg1(k-2*n-2,k+1) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg1(k+2*n+4,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg1(k+2*n+8,k+1) = -sqrt((dq+1.)*(dq+2.))/2.;
            }

            if (p > q) ++k;
        }
    }

    void setupGg2(DMatrix& Gg2, int order1, int order2)
    {
        Assert(int(Gg2.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(Gg2.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        Gg2.setZero();
        for(int n=0,k=0;n<=order2;++n) for(int p=n,q=0;p>=q;--p,++q,++k) {
            double dp = p;
            double dq = q;
            // d(psi(x,y))/dg2 = psi(x,y) Gg2
            // Gg2 is G_g2
            // d/dg2 = y d/dx + x d/dy
            //       = 1/2 i (ap^2 - aq^2 + apt^2 - aqt^2)
            // Gg2( p-2,q , pq ) = 1/2 i sqrt(p(p-1))
            // Gg2( p,q-2 , pq ) = -1/2 i sqrt(q(q-1))
            // Gg2( p+2,q , pq ) = 1/2 i sqrt((p+1)(p+2))
            // Gg2( p,q+2 , pq ) = -1/2 i sqrt((q+1)(q+2))
            if (p==q) {
                if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dp*(dp-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
            } else if (p==q+1) {
                if (p>1 && n<=order1+2) Gg2(k-2*n,k) = -sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg2(k+2*n+6,k) = sqrt((dq+1.)*(dq+2.))/2.;
                if (p>1 && n<=order1+2) Gg2(k-2*n-1,k+1) = -sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg2(k-2*n-3,k+1) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+3,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg2(k+2*n+5,k+1) = sqrt((dq+1.)*(dq+2.))/2.;
            } else if (p==q+2) {
                if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
                if (n<=order2+2) Gg2(k-2*n+1,k+1) = -sqrt(dp*(dp-1.));
                if (q>1 && n<=order1+2) Gg2(k-2*n-3,k+1) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+3,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg2(k+2*n+7,k+1) = sqrt((dq+1.)*(dq+2.));
            } else {
                Gg2(k-2*n+2,k) = sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg2(k-2*n-2,k) = -sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+4,k) = sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg2(k+2*n+8,k) = -sqrt((dq+1.)*(dq+2.))/2.;
                if (n<=order2+2) Gg2(k-2*n+1,k+1) = -sqrt(dp*(dp-1.))/2.;
                if (q>1 && n<=order1+2) Gg2(k-2*n-3,k+1) = sqrt(dq*(dq-1.))/2.;
                if (n+2<=order1) Gg2(k+2*n+3,k+1) = -sqrt((dp+1.)*(dp+2.))/2.;
                if (n+2<=order1) Gg2(k+2*n+7,k+1) = sqrt((dq+1.)*(dq+2.))/2.;
            }

            if (p > q) ++k;
        }
    }

    void setupGmu(DMatrix& Gmu, int order1, int order2)
    {
        Assert(int(Gmu.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(Gmu.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        Gmu.setZero();
        for(int n=0,k=0;n<=order2;++n) for(int p=n,q=0;p>=q;--p,++q,++k) {
            double dp = p;
            double dq = q;
            // d(psi(x,y))/dmu = psi(x,y) Gmu
            // d/dmu = x d/dx + y d/dy
            //         ap aq - apt aqt - 1
            // Gmu( p-1,q-1 , pq ) = sqrt(pq)
            // Gmu( p+1,q+1 , pq ) = -sqrt((p+1)(q+1))
            // Gmu( pq , pq ) = -1
            if (q>0 && n<=order1+2) Gmu(k-2*n-1,k) = sqrt(dp*dq);
            if (n+2<=order1) Gmu(k+2*n+5,k) = -sqrt((dp+1.)*(dq+1.));
            if (n<=order1) Gmu(k,k) = -1.;
            if (p > q) {
                if (q>0 && n<=order1+2) Gmu(k-2*n,k+1) = sqrt(dp*dq);
                if (n+2<=order1) Gmu(k+2*n+6,k+1) = -sqrt((dp+1.)*(dq+1.));
                if (n<=order1) Gmu(k+1,k+1) = -1.;
            }

            if (p > q) ++k;
        }
    }

    void setupGth(DMatrix& Gth, int order1, int order2)
    {
        Assert(int(Gth.TMV_colsize()) == (order1+1)*(order1+2)/2);
        Assert(int(Gth.TMV_rowsize()) == (order2+1)*(order2+2)/2);

        Gth.setZero();
        for(int n=0,k=0;n<=order2;++n) for(int p=n,q=0;p>=q;--p,++q,++k) {
            double dp = p;
            double dq = q;
            // d(psi(x,y))/dtheta = psi(x,y) Gth
            // d/dtheta = x d/dy - y d/dx
            //          = m i
            // Gth( pq , pq ) = m i
            if (p>q && n<=order1) {
                Gth(k+1,k) = (dp-dq);
                Gth(k,k+1) = -(dp-dq);
            }

            if (p > q) ++k;
        }
    }

}}}}
