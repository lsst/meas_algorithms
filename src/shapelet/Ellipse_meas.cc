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

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    bool Ellipse::doMeasure(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>* psf,
        int galOrder, int galOrder2, int maxm,
        double sigma, long& flag, double thresh, DMatrix* cov)
    {
        dbg<<"Start DoMeasure: galOrder = "<<galOrder<<", psf = "<<bool(psf)<<std::endl;
        dbg<<"fix = "<<_isFixedCen<<"  "<<_isFixedGamma<<"  "<<_isFixedMu<<std::endl;
        dbg<<"Thresh = "<<thresh<<std::endl;
        int npix = 0;
        for(size_t i=0;i<pix.size();++i) {
            xdbg<<"npix["<<i<<"] = "<<pix[i].size()<<std::endl;
            npix += pix[i].size();
        }

        int galSize = (galOrder+1)*(galOrder+2)/2;
        if (npix <= galSize) {
            dbg<<"Too few pixels ("<<npix<<") for given gal_order. \n";
            return false;
        }

        BVec b(galOrder,sigma);

        if (!doMeasureShapelet(pix,psf,b,galOrder,galOrder2,maxm)) {
            xdbg<<"Could not measure a shapelet vector.\n";
            return false;
        }
        if (b(0) <= 0) {
            xdbg<<"Bad flux in measured shapelet\n";
            return false;
        }
        xdbg<<"b = "<<b<<std::endl;

        return findRoundFrame(b,psf,galOrder2,thresh,flag,cov);
    }

    bool Ellipse::findRoundFrame(
        const BVec& b, bool psf, int galOrder2, double thresh,
        long& flag, DMatrix* cov)
    {
        DVector x(5);
        DVector f(5);

        x.setZero();

        EllipseSolver3 solver(b,galOrder2,_isFixedCen,_isFixedGamma,_isFixedMu);

#ifdef NOTHROW
        solver.noUseCholesky();
#endif
        double ftol = thresh*thresh;
        double gtol = thresh*ftol;
        solver.setTol(ftol,gtol);
        solver.setMinStep(gtol*thresh);
        solver.setOutput(*dbgout);
        if (XDEBUG) solver.useVerboseOutput();
        solver.setMinStep(1.e-6*gtol);
        solver.setDelta0(0.01);
        solver.setMaxIter(200);
        if (psf && !_isFixedMu) {
            solver.setDelta0(0.01);
            solver.useDogleg();
            solver.dontZeroB11();
            solver.useSVD();
        } else {
            solver.useHybrid();
        }
        if (solver.solve(x,f)) {
            dbg<<"Found good round frame:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            double f_normInf = f.TMV_normInf();
            if (psf && !_isFixedMu && !(f_normInf < solver.getFTol())) {
                xdbg<<"Oops, Local minimum, not real solution.\n";
                xdbg<<"f.norm = "<<f.norm()<<std::endl;
                xdbg<<"f.normInf = "<<f_normInf<<std::endl;
                xdbg<<"ftol = "<<solver.getFTol()<<std::endl;
                dbg<<"FLAG SHEAR_LOCAL_MIN\n";
                flag |= SHEAR_LOCAL_MIN;
                return false;
            }
        }  else {
            dbg<<"findRoundFrame solver failed:\n";
            dbg<<"x = "<<EIGEN_Transpose(x)<<std::endl;
            dbg<<"f = "<<EIGEN_Transpose(f)<<std::endl;
            //if (XDEBUG) if (!solver.testJ(x,f,dbgout,1.e-5)) exit(1);
            dbg<<"FLAG SHEAR_DIDNT_CONVERGE\n";
            flag |= SHEAR_DIDNT_CONVERGE;
            return false;
        }

        double sigma = b.getSigma();
        preShiftBy(std::complex<double>(x(0),x(1))*sigma,
                   std::complex<double>(x(2),x(3)),
                   x(4));

        dbg<<"ell => "<<*this<<std::endl;

        if (cov) {
            solver.useSVD();
            solver.getCovariance(*cov);
            xdbg<<"cov = "<<*cov<<std::endl;
        }

        return true;
    }


}}}}
