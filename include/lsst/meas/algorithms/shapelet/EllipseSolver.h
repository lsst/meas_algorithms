#ifndef MeasAlgoShapeletEllipseSolver_H
#define MeasAlgoShapeletEllipseSolver_H

#include "lsst/meas/algorithms/shapelet/Pixel.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
#include <vector>
#include "lsst/meas/algorithms/shapelet/NLSolver.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class EllipseSolver3 : public NLSolver
    {
    public :

        EllipseSolver3(
            const BVec& b0, int order,
            bool fixcen=false, bool fixgam=false, bool fixmu=false);
        ~EllipseSolver3();

        void calculateF(const DVector& x, DVector& f) const;
        void calculateJ(const DVector& x, const DVector& f, DMatrix& df) const;

        void useNumericJ();
        void dontZeroB11();
        void getCovariance(DMatrix& cov) const;
        void getInverseCovariance(DMatrix& invcov) const;

        // CallF takes x and f of length 5, rather than whatever shorter
        // length that F takex (depending on if things are fixed).
        void callF(const DVector& x, DVector& f) const;
        bool solve(DVector& x, DVector& f) const;
        bool testJ(const DVector& x, DVector& f,
                   std::ostream* os=0, double relerr=0.) const;

    private :

        struct ESImpl3;

        ESImpl3* _pimpl;
    };

}}}}

#endif
