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

#include "lsst/meas/algorithms/shapelet/EllipseSolver.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"

//#define JTEST
#define MAX_X_DEV 4.

#ifdef JTEST
#include "TestHelper.h"
bool shouldShowTests = true;
bool shouldThrow = false;
std::string lastSuccess = "";
std::ostream* testout = &std::cout;
#endif

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    struct EllipseSolver3::ESImpl3
    {

        ESImpl3(
            const BVec& b0, int order,
            bool isFixedCen, bool isFixedGamma, bool isFixedMu);

        void calculateF(const DVector& x, DVector& f) const;
        void calculateJ(const DVector& x, const DVector& f, DMatrix& J) const;

        void doF3(const DVector& x, DVector& f) const;
        void doJ3(const DVector& x, const DVector& f, DMatrix& J) const;

        void doF2(const DVector& x, DVector& f) const;
        void doJ2(const DVector& x, const DVector& f, DMatrix& J) const;

        void doF1(const DVector& x, DVector& f) const;
        void doJ1(const DVector& x, const DVector& f, DMatrix& J) const;

        int b0order, bxorder, bxorderp2;
        int b0size, bxsize, bxsizep2;
        mutable BVec b0;
        mutable DVector bx;
        mutable DVector bxsave;
        mutable DMatrix Daug;
        mutable DMatrix Saug;
        mutable DMatrix Taug;
        mutable DMatrixView D;
        mutable DMatrixView S;
        mutable DMatrixView T;
        mutable DVector dbdE;
        mutable DVector Db0;
        mutable DVector SDb0;
        mutable DVector GDb0;
        mutable DVector GthDb0;
        bool fixcen, fixgam, fixmu, numeric_j, zerob11;
        DMatrix U;
        mutable DVector xinit;
        mutable DVector xx;
        mutable DVector ff;
        mutable DMatrix jj;
        mutable DVector x_short;
        mutable DVector f_short;
        mutable double fixuc,fixvc,fixg1,fixg2,fixm;
        DMatrix Gx;
        DMatrix Gy;
        DMatrix Gg1;
        DMatrix Gg2;
        DMatrix Gth;
        DMatrix Gmu;
    };

    EllipseSolver3::EllipseSolver3(
        const BVec& b0, int order,
        bool fixcen, bool fixgam, bool fixmu) :
        _pimpl(new ESImpl3(b0,order,fixcen,fixgam,fixmu))
    {}

    EllipseSolver3::~EllipseSolver3() 
    { delete _pimpl; }

    void EllipseSolver3::calculateF(const DVector& x, DVector& f) const
    { _pimpl->calculateF(x,f); }

    void EllipseSolver3::calculateJ(
        const DVector& x, const DVector& f, DMatrix& J) const
    { 
#ifdef ALWAYS_NUMERIC_J
        NLSolver::calculateJ(x,f,J);
#else
        if (_pimpl->numeric_j) NLSolver::calculateJ(x,f,J);
        else _pimpl->calculateJ(x,f,J); 
#endif
    }

    EllipseSolver3::ESImpl3::ESImpl3(
        const BVec& _b0, int _order,
        bool _fixcen, bool _fixgam, bool _fixmu) :
        b0order(_b0.getOrder()), bxorder(_order), bxorderp2(bxorder+2),
        b0size(_b0.size()), bxsize((bxorder+1)*(bxorder+2)/2),
        bxsizep2((bxorderp2+1)*(bxorderp2+2)/2),
        b0(_b0), bx(bxsize), bxsave(bxsize),
        Daug(bxsize,bxsizep2), Saug(bxsize,bxsizep2),
        Taug(bxsize,bxsizep2),
        D(TMV_colRange(Daug,0,bxsize)), S(TMV_colRange(Saug,0,bxsize)),
        T(TMV_colRange(Taug,0,bxsize)),
        dbdE(6), Db0(bxsize), SDb0(bxsize), GDb0(bxsizep2), GthDb0(bxsizep2),
        fixcen(_fixcen), fixgam(_fixgam), fixmu(_fixmu),
        numeric_j(false), zerob11(true),
        U((fixcen?0:2)+(fixgam?0:2)+(fixmu?0:1),5), 
        xinit(5), xx(5), ff(5), jj(5,5),
        x_short(U.TMV_colsize()), f_short(U.TMV_colsize()),
        Gx(bxsizep2,bxsize), Gy(bxsizep2,bxsize), 
        Gg1(bxsizep2,bxsize), Gg2(bxsizep2,bxsize), 
        Gth(bxsizep2,bxsize), Gmu(bxsizep2,bxsize)
        {
            //xdbg<<"EllipseSolver3: \n";
            //xdbg<<"b0 = "<<b0.vec()<<std::endl;
            //xdbg<<"b0.order, sigma = "<<b0.getOrder()<<"  "<<b0.getSigma()<<std::endl;
            //xdbg<<"order = "<<_order<<std::endl;
            //xdbg<<"fixcen,gam,mu = "<<fixcen<<"  "<<fixgam<<"  "<<fixmu<<std::endl;
            //xdbg<<"bxorder = "<<bxorder<<"  "<<bxsize<<std::endl;
            //xdbg<<"bxorderp2 = "<<bxorderp2<<"  "<<bxsizep2<<std::endl;
            U.setZero();
            xinit.setZero();

            int k=0;
            if (!fixcen) { U(k++,0) = 1.; U(k++,1) = 1.; }
            if (!fixgam) { U(k++,2) = 1.; U(k++,3) = 1.; }
            if (!fixmu) { U(k,4) = 1.; }

            Daug.setZero();
            Saug.setZero();
            Taug.setZero();

            setupGx(Gx,bxorderp2,bxorder);
            setupGy(Gy,bxorderp2,bxorder);
            setupGg1(Gg1,bxorderp2,bxorder);
            setupGg2(Gg2,bxorderp2,bxorder);
            setupGth(Gth,bxorderp2,bxorder);
            setupGmu(Gmu,bxorderp2,bxorder);

            if (fixcen) { 
                T.TMV_diag().TMV_setAllTo(1.);
            }
            if (fixgam) { 
                S.TMV_diag().TMV_setAllTo(1.);
            }
            if (fixmu) { 
                D.TMV_diag().TMV_setAllTo(1.);
            }
        }

    void EllipseSolver3::ESImpl3::doF3(const DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);

        //xdbg<<"Start doF3\n";
        //xdbg<<"x = "<<x<<std::endl;

#ifdef MAX_X_DEV
        if ((x-xinit).TMV_normInf() > MAX_X_DEV) {
            dbg<<"bad x-xinit: "<<(x-xinit)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            f = 2.* bx.TMV_subVector(1,6) / bx(0);
            return;
        }
#endif

        std::complex<double> zc(x[0],x[1]);
        std::complex<double> g(x[2],x[3]);
        double mu = x[4];

        if (norm(g) > 0.99) {
            dbg<<"bad gsq = "<<norm(g)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            f = 2.* bx.TMV_subVector(1,6) / bx(0);
            return;
        }

        bxsave = bx;
        bool bxset = false;

        //xdbg<<"Start with b0 = "<<b0.vec()<<std::endl;
        //xdbg<<"bx = "<<bx<<std::endl;
        //TODO: Switch to doing calculate??Transform calls for only the portion
        //      of the matrix that I acutally use.
        //      This has implications for doJ though...
        if (!fixmu) {
            //xdbg<<"calc Mu:\n";
            calculateMuTransform(mu,bxorder,Daug);
            bx = TMV_colRange(D,0,b0size) * b0.vec();
            bxset = true;
            //xdbg<<"mu = "<<mu<<" bx => "<<bx<<std::endl;
        } else if (fixm != 0.) {
            bx = TMV_colRange(D,0,b0size) * b0.vec();
            bxset = true;
        }
        if (!fixgam) {
            //xdbg<<"calc G:\n";
            calculateGTransform(g,bxorder,Saug);
            if (bxset) {
                bx = S * bx;
            } else {
                bx = TMV_colRange(S,0,b0size) * b0.vec();
                bxset = true;
            }
            //xdbg<<"g = "<<g<<" bx => "<<bx<<std::endl;
        } else if (fixg1 != 0. || fixg2 != 0.) {
            if (bxset) {
                bx = S * bx;
            } else {
                bx = TMV_colRange(S,0,b0size) * b0.vec();
                bxset = true;
            }
        }
        if (!fixcen) {
            //xdbg<<"calc Z:\n";
            calculateZTransform(zc,bxorder,Taug);
            //xdbg<<"after calc Z\n";
            if (bxset) {
                //xdbg<<"T = "<<T<<std::endl;
                bx = T * bx;
            } else {
                //xdbg<<"T = "<<TMV_colRange(T,0,b0size)<<std::endl;
                bx = TMV_colRange(T,0,b0size) * b0.vec();
                bxset = true;
            }
            //xdbg<<"zc = "<<zc<<" bx => "<<bx<<std::endl;
        } else if (fixuc != 0. || fixvc != 0.) {
            if (bxset) {
                bx = T * bx;
            } else {
                bx = TMV_colRange(T,0,b0size) * b0.vec();
                bxset = true;
            }
        }

        if (bx(0) <= 0.) {
            dbg<<"bad bx(0): "<<bx(0)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            xdbg<<"bx = "<<bx<<std::endl;
            f = 2.* bxsave.TMV_subVector(1,6) / bxsave(0);
            bx = bxsave;
            return;
        }

        f = bx.TMV_subVector(1,6) / bx(0);

        if (!zerob11) f(4) = 0.;
    }

    void EllipseSolver3::ESImpl3::doF2(const DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);

        //xdbg<<"Start doF3\n";
        //xdbg<<"x = "<<x<<std::endl;

#ifdef MAX_X_DEV
        if ((x-xinit).TMV_normInf() > MAX_X_DEV) {
            dbg<<"bad x-xinit: "<<(x-xinit)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            f = 2.* bx.TMV_subVector(1,6) / bx(0);
            return;
        }
#endif

        std::complex<double> zc(x[0],x[1]);
        double mu = x[4];

        bxsave = bx;
        bool bxset = false;

        //xdbg<<"Start with b0 = "<<b0.vec()<<std::endl;
        //xdbg<<"bx = "<<bx<<std::endl;
        if (!fixmu) {
            //xdbg<<"calc Mu:\n";
            calculateMuTransform(mu,bxorder,Daug);
            bx = TMV_colRange(D,0,b0size) * b0.vec();
            bxset = true;
            //xdbg<<"mu = "<<mu<<" bx => "<<bx<<std::endl;
        } else if (fixm != 0.) {
            bx = TMV_colRange(D,0,b0size) * b0.vec();
            bxset = true;
        }
        if (!fixcen) {
            //xdbg<<"calc Z:\n";
            calculateZTransform(zc,bxorder,Taug);
            //xdbg<<"after calc Z\n";
            if (bxset) {
                //xdbg<<"T = "<<T<<std::endl;
                bx = T * bx;
            } else {
                //xdbg<<"T = "<<TMV_colRange(T,0,b0size)<<std::endl;
                bx = TMV_colRange(T,0,b0size) * b0.vec();
                bxset = true;
            }
            //xdbg<<"zc = "<<zc<<" bx => "<<bx<<std::endl;
        } else if (fixuc != 0. || fixvc != 0.) {
            if (bxset) {
                bx = T * bx;
            } else {
                bx = TMV_colRange(T,0,b0size) * b0.vec();
                bxset = true;
            }
        }

        if (bx(0) <= 0.) {
            dbg<<"bad bx(0): "<<bx(0)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            xdbg<<"bx = "<<x<<std::endl;
            f = 2.* bxsave.TMV_subVector(1,6) / bxsave(0);
            bx = bxsave;
            return;
        }

        //xdbg<<"Done: bx = "<<bx<<std::endl;
        f = bx.TMV_subVector(1,6) / bx(0);

        if (!zerob11) f(4) = 0.;
    }

    void EllipseSolver3::ESImpl3::doF1(const DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);

        //xdbg<<"Start doF3\n";
        //xdbg<<"x = "<<x<<std::endl;

#ifdef MAX_X_DEV
        if ((x-xinit).TMV_normInf() > MAX_X_DEV) {
            dbg<<"bad x-xinit: "<<(x-xinit)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            f = 2.* bx.TMV_subVector(1,6) / bx(0);
            return;
        }
#endif

        std::complex<double> zc(x[0],x[1]);

        bxsave = bx;

        //xdbg<<"Start with b0 = "<<b0.vec()<<std::endl;
        //xdbg<<"bx = "<<bx<<std::endl;
        if (!fixcen) {
            //xdbg<<"calc Z:\n";
            calculateZTransform(zc,bxorder,Taug);
            //xdbg<<"after calc Z\n";
            bx = TMV_colRange(T,0,b0size) * b0.vec();
            //xdbg<<"zc = "<<zc<<" bx => "<<bx<<std::endl;
        } else if (fixuc != 0. || fixvc != 0.) {
            bx = TMV_colRange(T,0,b0size) * b0.vec();
        }

        if (bx(0) <= 0.) {
            dbg<<"bad bx(0): "<<bx(0)<<std::endl;
            xdbg<<"x = "<<x<<std::endl;
            xdbg<<"bx = "<<x<<std::endl;
            f = 2.* bxsave.TMV_subVector(1,6) / bxsave(0);
            bx = bxsave;
            return;
        }

        f = bx.TMV_subVector(1,6) / bx(0);

        if (!zerob11) f(4) = 0.;
    }

    void EllipseSolver3::ESImpl3::doJ3(
        const DVector& x, const DVector& f, DMatrix& J) const
    {
        // bx = T S D b0
        //
        // dbx/dzc = (dT/dzc) S D b0
        // dbx/dg  = T (dS/dg) D b0
        // dbx/dmu = T S (dD/dmu) b0
        //
        // f = bx(1:6) / bx(0)
        // df/dE = dbx/dE(1:6) / bx(0) - (1/bx(0)^2) (dbx(0)/dE) bx(1:6)
        //       = ( dbx/dE(1:6) - dbx/dE(0) f ) / bx(0)

        //xdbg<<"Start doJ3\n";
        //xdbg<<"x = "<<x<<std::endl;
        //xdbg<<"f = "<<f<<std::endl;

        Assert(x.size() == 5);
        Assert(J.TMV_rowsize() == 5);
        Assert(J.TMV_colsize() == 5);

        std::complex<double> zc(x[0],x[1]);
        std::complex<double> g(x[2],x[3]);
        double mu = x[4];
        double fact = 1./(1.-norm(g));

        if (!fixcen || !fixgam) {
            // Both cen and gam use this:
            Db0 = TMV_colRange(D,0,b0size) * b0.vec();
        }

        // Leave off one factor of bx(0) for now.
        // We apply it at the end to the whole matrix.
        if (!fixcen) {
            // dT/dx = T Gx;
            // dT/dy = T Gy;
            // db/dx = T Gx S D b0
            // db/dy = T Gy S D b0
            augmentZTransformCols(zc,bxorder,Taug);
            SDb0 = S * Db0;
            GDb0 = Gx * SDb0;
            dbdE = TMV_rowRange(Taug,0,6) * GDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(0) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(0) = "<<J.col(0)<<std::endl;
            GDb0 = Gy * SDb0;
            dbdE = TMV_rowRange(Taug,0,6) * GDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(1) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(1) = "<<J.col(1)<<std::endl;
        } else {
            J.col(0).setZero();
            //dbg<<"J(0) = "<<J.col(0)<<std::endl;
            J.col(1).setZero();
            //dbg<<"J(1) = "<<J.col(1)<<std::endl;
        }

        if (!fixgam) {
            // dS/dg1 = fact * S * (Gg1 + g2 * Gth);
            // dS/dg2 = fact * S * (Gg2 - g1 * Gth);
            augmentGTransformCols(g,bxorder,Saug);
            double g1 = real(g);
            double g2 = imag(g);
            GDb0 = Gg1 * Db0;
            GthDb0 = Gth * Db0;
            GDb0 += g2*GthDb0;
            SDb0 = fact * Saug * GDb0;
            dbdE = TMV_rowRange(T,0,6) * SDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(2) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(2) = "<<J.col(2)<<std::endl;
            GDb0 = Gg2 * Db0;
            GDb0 -= g1*GthDb0;
            SDb0 = fact * Saug * GDb0;
            dbdE = TMV_rowRange(T,0,6) * SDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(3) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(3) = "<<J.col(3)<<std::endl;
        } else {
            J.col(2).setZero();
            //dbg<<"J(2) = "<<J.col(2)<<std::endl;
            J.col(3).setZero();
            //dbg<<"J(3) = "<<J.col(3)<<std::endl;
        }

        if (!fixmu) {
            // dD/dmu = D Gmu
            augmentMuTransformCols(mu,bxorder,Daug);
            GDb0 = TMV_colRange(Gmu,0,b0size) * b0.vec();
            Db0 = Daug * GDb0;
            SDb0 = S * Db0; 
            dbdE = TMV_rowRange(T,0,6) * SDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(4) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(4) = "<<J.col(4)<<std::endl;
        } else {
            J.col(4).setZero();
            //dbg<<"J(4) = "<<J.col(4)<<std::endl;
        }
        J /= bx(0);
        if (!zerob11) J.row(4).setZero();
        //xdbg<<"J = "<<J<<std::endl;

#ifdef JTEST
        double dE = 1.e-6;
        dbg<<"dE = "<<dE<<std::endl;
        testout = dbgout;

        // This section was only used for debugging purposes, but
        // I'm leaving it in, since the derivations of Gx[i] are
        // helpful in understanding some of the code following this section.

        //
        // dD/dmu:
        //
        DMatrix Dbig(bxsizep2,bxsizep2,0.);
        DMatrix D0(bxsize,bxsize,0.);
        DMatrix D1(bxsize,bxsize,0.);
        DMatrix D2(bxsize,bxsize,0.);
        calculateMuTransform(mu,bxorderp2,Dbig);
        calculateMuTransform(mu,bxorder,D0);
        calculateMuTransform(mu-dE,bxorder,D1);
        calculateMuTransform(mu+dE,bxorder,D2);
        DMatrix Gmubig(bxsizep2,bxsize,0.);
        setupGmu(Gmubig,bxorderp2,bxorder);

        DMatrix dDdmu_num = (D2-D1)/(2.*dE);
        DMatrix d2D_num = (D2+D1-2.*D0)/(dE*dE);
        DMatrix dDdmu = TMV_rowRange(Dbig,0,bxsize) * Gmubig;
        dbg<<"Gmu = "<<Gmubig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"D = "<<Dbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dDdmu = "<<dDdmu.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dDdmu_num = "<<dDdmu_num.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Norm(dD/dmu - numeric dD/dmu) = "<<Norm(dDdmu-dDdmu_num)<<std::endl;
        dbg<<"dE*Norm(d2D_num) = "<<dE*Norm(d2D_num)<<std::endl;
        test(Norm(dDdmu-dDdmu_num) < 30.*dE*Norm(d2D_num),"dDdmu");

        //
        // dS/dg1
        //
        DMatrix Sbig(bxsizep2,bxsizep2,0.);
        DMatrix S0(bxsize,bxsize,0.);
        DMatrix S1(bxsize,bxsize,0.);
        DMatrix S2(bxsize,bxsize,0.);
        calculateGTransform(g,bxorderp2,Sbig);
        calculateGTransform(g,bxorder,S0);
        calculateGTransform(g-std::complex<double>(dE,0),bxorder,S1);
        calculateGTransform(g+std::complex<double>(dE,0),bxorder,S2);
        DMatrix Gg1big(bxsizep2,bxsize,0.);
        DMatrix Gthbig(bxsizep2,bxsize,0.);
        setupGg1(Gg1big,bxorderp2,bxorder);
        setupGth(Gthbig,bxorderp2,bxorder);

        DMatrix dSdg1_num = (S2-S1)/(2.*dE);
        DMatrix d2S_num = (S2+S1-2.*S0)/(dE*dE);
        DMatrix dSdg1 = fact * TMV_rowRange(Sbig,0,bxsize) * (Gg1big + imag(g) * Gthbig);
        dbg<<"Gg1 = "<<Gg1big.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Gth = "<<Gthbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"S = "<<Sbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dSdg1 = "<<dSdg1.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dSdg1_num = "<<dSdg1_num.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Norm(dS/dg1 - numeric dS/dg1) = "<<Norm(dSdg1-dSdg1_num)<<std::endl;
        dbg<<"dE*Norm(d2S_num) = "<<dE*Norm(d2S_num)<<std::endl;
        test(Norm(dSdg1-dSdg1_num) < 30.*dE*Norm(d2S_num),"dSdg1");

        //
        // dS/dg2
        //
        calculateGTransform(g-std::complex<double>(0,dE),bxorder,S1);
        calculateGTransform(g+std::complex<double>(0,dE),bxorder,S2);
        DMatrix Gg2big(bxsizep2,bxsize,0.);
        setupGg2(Gg2big,bxorderp2,bxorder);

        DMatrix dSdg2_num = (S2-S1)/(2.*dE);
        d2S_num = (S2+S1-2.*S0)/(dE*dE);
        DMatrix dSdg2 = fact * TMV_rowRange(Sbig,0,bxsize) * (Gg2big - real(g) * Gthbig);
        dbg<<"Gg2 = "<<Gg2big.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Gth = "<<Gthbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"S = "<<Sbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dSdg2 = "<<dSdg2.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dSdg2_num = "<<dSdg2_num.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Norm(dS/dg2 - numeric dS/dg2) = "<<Norm(dSdg2-dSdg2_num)<<std::endl;
        dbg<<"dE*Norm(d2S_num) = "<<dE*Norm(d2S_num)<<std::endl;
        test(Norm(dSdg2-dSdg2_num) < 30.*dE*Norm(d2S_num),"dSdg2");

        //
        // dT/dx
        //
        DMatrix Tbig(bxsizep2,bxsizep2,0.);
        DMatrix T0(bxsize,bxsize,0.);
        DMatrix T1(bxsize,bxsize,0.);
        DMatrix T2(bxsize,bxsize,0.);
        calculateZTransform(zc,bxorderp2,Tbig);
        calculateZTransform(zc,bxorder,T0);
        calculateZTransform(zc-std::complex<double>(dE,0),bxorder,T1);
        calculateZTransform(zc+std::complex<double>(dE,0),bxorder,T2);
        DMatrix Gxbig(bxsizep2,bxsize,0.);
        setupGx(Gxbig,bxorderp2,bxorder);

        DMatrix dTdx_num = (T2-T1)/(2.*dE);
        DMatrix d2T_num = (T2+T1-2.*T0)/(dE*dE);
        DMatrix dTdx = TMV_rowRange(Tbig,0,bxsize) * Gxbig;
        dbg<<"Gx = "<<Gxbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"T = "<<Tbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dTdx = "<<dTdx.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dTdx_num = "<<dTdx_num.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Norm(dT/dx - numeric dT/dx) = "<<Norm(dTdx-dTdx_num)<<std::endl;
        dbg<<"dE*Norm(d2T_num) = "<<dE*Norm(d2T_num)<<std::endl;
        test(Norm(dTdx-dTdx_num) < 30.*dE*Norm(d2T_num),"dTdx");

        //
        // dT/dy
        //
        calculateZTransform(zc-std::complex<double>(0,dE),bxorder,T1);
        calculateZTransform(zc+std::complex<double>(0,dE),bxorder,T2);
        DMatrix Gybig(bxsizep2,bxsize,0.);
        setupGy(Gybig,bxorderp2,bxorder,0.);

        DMatrix dTdy_num = (T2-T1)/(2.*dE);
        d2T_num = (T2+T1-2.*T0)/(dE*dE);
        DMatrix dTdy = TMV_rowRange(Tbig,0,bxsize) * Gybig;
        dbg<<"Gy = "<<Gybig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"T = "<<Tbig.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dTdy = "<<dTdy.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"dTdy_num = "<<dTdy_num.TMV_subMatrix(0,6,0,6)<<std::endl;
        dbg<<"Norm(dT/dy - numeric dT/dy) = "<<Norm(dTdy-dTdy_num)<<std::endl;
        dbg<<"dE*Norm(d2T_num) = "<<dE*Norm(d2T_num)<<std::endl;
        test(Norm(dTdy-dTdy_num) < 30.*dE*Norm(d2T_num),"dTdy");

        //
        // J 
        //

        DMatrix J_num(5,5,0.);
        DVector x0 = x;
        DVector f0(5);
        doF3(x0,f0);
        DVector bx0 = bx.TMV_subVector(0,6);
        dbg<<"x0 = "<<x0<<std::endl;
        dbg<<"b0 = "<<bx0<<std::endl;
        dbg<<"f0 = "<<f0<<std::endl;

        for(int k=0;k<5;k++) {
            dbg<<"k = "<<k<<std::endl;

            DVector f1(5);
            DVector x1 = x; x1(k) -= dE;
            doF3(x1,f1);
            DVector b1 = bx.TMV_subVector(0,6);
            dbg<<"x1 = "<<x1<<std::endl;
            dbg<<"b1 = "<<b1<<std::endl;
            dbg<<"f1 = "<<f1<<std::endl;

            DVector f2(5);
            DVector x2 = x; x2(k) += dE;
            doF3(x2,f2);
            DVector b2 = bx.TMV_subVector(0,6);
            dbg<<"x2 = "<<x2<<std::endl;
            dbg<<"b2 = "<<b2<<std::endl;
            dbg<<"f2 = "<<f2<<std::endl;

            DVector dbdE_num = (b2-b1)/(2.*dE);
            dbg<<"dbdE_num = "<<dbdE_num<<std::endl;
            DVector d2bdE = (b2+b1-2.*bx0)/(dE*dE);
            dbg<<"d2bdE = "<<d2bdE<<std::endl;

            DVector dbdE(6);
            double g1 = real(g);
            double g2 = imag(g);
            Db0 = TMV_colRange(D,0,b0size) * b0.vec();
            dbg<<"Db0 = "<<Db0<<std::endl;
            switch (k) {
              case 0: SDb0 = S * Db0;
                      dbg<<"SDb0 = "<<SDb0<<std::endl;
                      GDb0 = Gx * SDb0;
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      dbdE = TMV_rowRange(Tbig,0,6) * GDb0;
                      dbg<<"dbdx = "<<dbdE<<std::endl;
                      dbg<<"dbdx = dTdx S D b0 = "<<
                          TMV_rowRange(dTdx,0,6) * S * Db0<<std::endl;
                      break;
              case 1: SDb0 = S * Db0;
                      dbg<<"SDb0 = "<<SDb0<<std::endl;
                      GDb0 = Gy * SDb0;
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      dbdE = TMV_rowRange(Tbig,0,6) * GDb0;
                      dbg<<"dbdy = "<<dbdE<<std::endl;
                      dbg<<"dbdy = dTdy S D b0 = "<<
                          TMV_rowRange(dTdy,0,6) * S * Db0<<std::endl;
                      break;
              case 2: GDb0 = Gg1 * Db0;
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      GthDb0 = Gth * Db0;
                      dbg<<"GthDb0 = "<<GthDb0<<std::endl;
                      GDb0 += g2*GthDb0;
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      SDb0 = fact * TMV_rowRange(Sbig,0,bxsize) * GDb0;
                      dbg<<"SDb0 = "<<SDb0<<std::endl;
                      dbdE = TMV_rowRange(T,0,6) * SDb0;
                      dbg<<"dbdg1 = "<<dbdE<<std::endl;
                      dbg<<"dbdg1 = T dSdg1 D b0 = "<<
                          TMV_rowRange(T,0,6) * dSdg1 * Db0<<std::endl;
                      break;
              case 3: GDb0 = Gg2 * Db0;
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      GthDb0 = Gth * Db0;
                      dbg<<"GthDb0 = "<<GthDb0<<std::endl;
                      GDb0 -= g1*GthDb0;
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      SDb0 = fact * TMV_rowRange(Sbig,0,bxsize) * GDb0;
                      dbg<<"SDb0 = "<<SDb0<<std::endl;
                      dbdE = TMV_rowRange(T,0,6) * SDb0;
                      dbg<<"dbdg2 = "<<dbdE<<std::endl;
                      dbg<<"dbdg2 = T dSdg2 D b0 = "<<
                          TMV_rowRange(T,0,6) * dSdg2 * Db0<<std::endl;
                      break;
              case 4: GDb0 = TMV_colRange(Gmu,0,b0size) * b0.vec();
                      dbg<<"GDb0 = "<<GDb0<<std::endl;
                      Db0 = TMV_rowRange(Dbig,0,bxsize) * GDb0;
                      dbg<<"Db0 = "<<Db0<<std::endl;
                      SDb0 = S * Db0; 
                      dbg<<"SDb0 = "<<SDb0<<std::endl;
                      dbdE = TMV_rowRange(T,0,6) * SDb0;
                      dbg<<"dbdmu = "<<dbdE<<std::endl;
                      dbg<<"dbdmu = T S dDdmu b0 = "<<
                          TMV_rowRange(T,0,6) * S * TMV_colRange(dDdmu,0,b0size) * b0.vec()
                          <<std::endl;
                      break;
            }
            dbg<<"dbdE = "<<dbdE<<std::endl;
            dbg<<"Norm(dbdE-dbdE_num) = "<<Norm(dbdE-dbdE_num)<<std::endl;
            dbg<<"dE*Norm(d2bdE) = "<<dE*Norm(d2bdE)<<std::endl;
            test(Norm(dbdE-dbdE_num) < 30.*dE*Norm(d2bdE),"dbdE");

            DVector dfdE = dbdE.TMV_subVector(1,6) - dbdE(0) * f0;
            dfdE /= bx(0);
            dbg<<"dfdE from dbdE = "<<dfdE<<std::endl;
            dbg<<"dfdE from dbdE_num = "<<
                (dbdE_num.TMV_subVector(1,6) - dbdE_num(0) * f0)/bx(0)<<std::endl;;
            dbg<<"dfdE in J = "<<J.col(k)<<std::endl;
            J_num.col(k) = (f2-f1)/(2.*dE);
            dbg<<"dfdE_num = "<<J_num.col(k)<<std::endl;
            DVector d2f_num = (f2+f1-2.*f0)/(dE*dE);
            dbg<<"d2f_num = "<<d2f_num<<std::endl;
            dbg<<"Norm(dfdE_num-dfdE) ="<<Norm(J_num.col(k)-dfdE)<<std::endl;
            dbg<<"dE*Norm(d2f_num) ="<<dE*Norm(d2f_num)<<std::endl;
            test(Norm(J_num.col(k)-dfdE) < 30.*dE*Norm(d2f_num),"dfdE");
        }

        dbg<<"J_num = "<<J_num<<std::endl;
        dbg<<"analytic: "<<J<<std::endl;
        if (fixcen) TMV_colRange(J_num,0,2).setZero(); 
        if (fixgam) TMV_colRange(J_num,2,4).setZero(); 
        if (fixmu) J_num.col(4).setZero(); 
        dbg<<"J_num => "<<J_num<<std::endl;
        dbg<<"Norm(diff) = "<<Norm(J-J_num)<<std::endl;
        dbg<<"dE*Norm(J) = "<<dE*Norm(J_num)<<std::endl;
        test(Norm(J_num-J) < dE*Norm(J_num),"dfdE");
        doF3(x,f0); // return everything to original values
#endif
    }

    void EllipseSolver3::ESImpl3::doJ2(
        const DVector& x, const DVector& f, DMatrix& J) const
    {
        //xdbg<<"Start doJ2\n";
        //xdbg<<"x = "<<x<<std::endl;
        //xdbg<<"f = "<<f<<std::endl;

        Assert(x.size() == 5);
        Assert(J.TMV_rowsize() == 5);
        Assert(J.TMV_colsize() == 5);

        std::complex<double> zc(x[0],x[1]);
        double mu = x[4];

        J.setZero();
        if (!fixcen) {
            augmentZTransformCols(zc,bxorder,Taug);
            Db0 = TMV_colRange(D,0,b0size) * b0.vec();
            //dbg<<"Db0 = "<<Db0<<std::endl;
            GDb0 = Gx * Db0;
            dbdE = TMV_rowRange(Taug,0,6) * GDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(0) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(0) = "<<J.col(0)<<std::endl;
            GDb0 = Gy * Db0;
            dbdE = TMV_rowRange(Taug,0,6) * GDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(1) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(1) = "<<J.col(1)<<std::endl;
        } else {
            //dbg<<"J(0) = "<<J.col(0)<<std::endl;
            //dbg<<"J(1) = "<<J.col(1)<<std::endl;
        }

        if (!fixmu) {
            // dD/dmu = D Gmu
            augmentMuTransformCols(mu,bxorder,Daug);
            GDb0 = TMV_colRange(Gmu,0,b0size) * b0.vec();
            //dbg<<"GDb0 = "<<GDb0<<std::endl;
            Db0 = Daug * GDb0;
            //dbg<<"Db0 = "<<Db0<<std::endl;
            SDb0 = S * Db0; 
            //dbg<<"SDb0 = "<<SDb0<<std::endl;
            dbdE = TMV_rowRange(T,0,6) * SDb0;
            //dbg<<"dbdE = "<<dbdE<<std::endl;
            J.col(4) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
            //dbg<<"J(4) = "<<J.col(4)<<std::endl;
        } else {
            //dbg<<"J(4) = "<<J.col(4)<<std::endl;
        }
        //xdbg<<"bx(0) = "<<bx(0)<<std::endl;
        J /= bx(0);
        if (!zerob11) J.row(4).setZero();
        //xdbg<<"J = "<<J<<std::endl;
    }

    void EllipseSolver3::ESImpl3::doJ1(
        const DVector& x, const DVector& f, DMatrix& J) const
    {
        //xdbg<<"Start doJ1\n";
        //xdbg<<"x = "<<x<<std::endl;
        //xdbg<<"f = "<<f<<std::endl;

        Assert(x.size() == 5);
        Assert(J.TMV_rowsize() == 5);
        Assert(J.TMV_colsize() == 5);

        std::complex<double> zc(x[0],x[1]);

        J.setZero();

        augmentZTransformCols(zc,bxorder,Taug);
        GDb0 = TMV_colRange(Gx,0,b0size) * b0.vec();
        dbdE = TMV_rowRange(Taug,0,6) * GDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(0) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(0) = "<<J.col(0)<<std::endl;
        GDb0 = TMV_colRange(Gy,0,b0size) * b0.vec();
        dbdE = TMV_rowRange(Taug,0,6) * GDb0;
        //dbg<<"dbdE = "<<dbdE<<std::endl;
        J.col(1) = dbdE.TMV_subVector(1,6) - dbdE(0) * f;
        //dbg<<"J(1) = "<<J.col(1)<<std::endl;

        J /= bx(0);
        if (!zerob11) J.row(4).setZero();
        //xdbg<<"J = "<<J<<std::endl;
    }

    void EllipseSolver3::ESImpl3::calculateF(const DVector& x, DVector& f) const
    {
        Assert(x.size() == U.TMV_colsize());
        Assert(f.size() == U.TMV_colsize());

        Assert(xx.size() == U.TMV_rowsize());
        EIGEN_Transpose(xx) = EIGEN_Transpose(x)*U;
        if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
        if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
        if (fixmu) { xx[4] = fixm; }

        if (!fixgam || fixg1 != 0. || fixg2 != 0.)
            doF3(xx,ff);
        else if (!fixmu || fixm != 0.) 
            doF2(xx,ff);
        else
            doF1(xx,ff);

        f = U*ff;
    }

    void EllipseSolver3::ESImpl3::calculateJ(
        const DVector& x, const DVector& f, DMatrix& j) const
    {
        Assert(x.size() == U.TMV_colsize());
        Assert(f.size() == U.TMV_colsize());
        Assert(j.TMV_rowsize() == U.TMV_colsize());
        Assert(j.TMV_colsize() == U.TMV_colsize());

        EIGEN_Transpose(xx) = EIGEN_Transpose(x)*U;
        if (fixcen) { xx[0] = fixuc;  xx[1] = fixvc; }
        if (fixgam) { xx[2] = fixg1;  xx[3] = fixg2; }
        if (fixmu) { xx[4] = fixm; }

        if (!fixgam || fixg1 != 0. || fixg2 != 0.)
            doJ3(xx,ff,jj);
        else if (!fixmu || fixm != 0.) 
            doJ2(xx,ff,jj);
        else
            doJ1(xx,ff,jj);

        j = U*jj*U.transpose();
    }

    void EllipseSolver3::useNumericJ() { _pimpl->numeric_j = true; }

    void EllipseSolver3::dontZeroB11() { _pimpl->zerob11 = false; this->useSVD(); }

    void EllipseSolver3::getCovariance(DMatrix& cov) const 
    {
        DMatrix cov1(_pimpl->x_short.size(),_pimpl->x_short.size());
        NLSolver::getCovariance(cov1);
        cov = _pimpl->U.transpose()*cov1*_pimpl->U;
        dbg<<"getCovariance:\n";
        dbg<<"cov1 = "<<cov1<<std::endl;
        dbg<<"full cov = "<<cov<<std::endl;
    }

    void EllipseSolver3::getInverseCovariance(DMatrix& invcov) const 
    {
        DMatrix invcov1(_pimpl->x_short.size(),_pimpl->x_short.size());
        NLSolver::getInverseCovariance(invcov1);
        invcov = _pimpl->U.transpose()*invcov1*_pimpl->U;
        dbg<<"getInverseCovariance:\n";
        dbg<<"invcov1 = "<<invcov1<<std::endl;
        dbg<<"full invcov = "<<invcov<<std::endl;
    }

    void EllipseSolver3::callF(const DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);
        _pimpl->xinit = x;
        if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
        if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
        if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

        _pimpl->x_short = _pimpl->U * x;

        calculateF(_pimpl->x_short,_pimpl->f_short);

        EIGEN_Transpose(f) = EIGEN_Transpose(_pimpl->f_short)*_pimpl->U;
    }

    bool EllipseSolver3::solve(DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);
        _pimpl->xinit = x;
        _pimpl->bx.TMV_setAllTo(1.); // In case we bail out with f = 2*bx/b(0)
        if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
        if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
        if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }
        if (_pimpl->fixcen && (x[0] != 0. || x[1] != 0.)) { 
            calculateZTransform(
                std::complex<double>(x[0],x[1]),_pimpl->bxorder,_pimpl->Taug);
        }
        if (_pimpl->fixgam && (x[2] != 0. || x[3] != 0.)) { 
            calculateGTransform(
                std::complex<double>(x[2],x[3]),_pimpl->bxorder,_pimpl->Saug);
        }
        if (_pimpl->fixmu && x[4] != 0.) { 
            calculateMuTransform(x[4],_pimpl->bxorder,_pimpl->Daug);
        }

        _pimpl->x_short = _pimpl->U * x;

        bool ret;
        try {
            ret = NLSolver::solve(_pimpl->x_short,_pimpl->f_short);
        } catch (...) {
            xdbg<<"Caught exception during NLSolver::solve"<<std::endl;
            ret = false;
        }

        EIGEN_Transpose(x) = EIGEN_Transpose(_pimpl->x_short)*_pimpl->U;
        if (_pimpl->fixcen) { x[0] = _pimpl->fixuc; x[1] = _pimpl->fixvc; }
        if (_pimpl->fixgam) { x[2] = _pimpl->fixg1; x[3] = _pimpl->fixg2; }
        if (_pimpl->fixmu) { x[4] = _pimpl->fixm; }
        EIGEN_Transpose(f) = EIGEN_Transpose(_pimpl->f_short)*_pimpl->U;

        return ret;
    }

    bool EllipseSolver3::testJ(
        const DVector& x, DVector& f,
        std::ostream* os, double relerr) const 
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);
        _pimpl->xinit = x;
        if (_pimpl->fixcen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
        if (_pimpl->fixgam) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
        if (_pimpl->fixmu) { _pimpl->fixm = x[4]; }

        _pimpl->x_short = _pimpl->U * x;

        return NLSolver::testJ(_pimpl->x_short,_pimpl->f_short,os,relerr);
    }


}}}}
