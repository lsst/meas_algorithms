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

    class BaseEllipseSolver : public NLSolver
    {
    public :

        BaseEllipseSolver() {}
        virtual ~BaseEllipseSolver() {}
        virtual void useNumericJ() = 0;
        virtual const BVec& getB() const = 0;
        virtual void callF(const DVector& x, DVector& f) const = 0;
        virtual void getBCov(DMatrix& bcov) const = 0;
    };

    class EllipseSolver : public BaseEllipseSolver
    {
    public :

        EllipseSolver(
            const std::vector<PixelList>& pix, 
            int order, double sigma, 
            bool fixcen=false, bool fixgam=false, bool fixmu=false,
            bool useflux=false);
        EllipseSolver(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>& psf, double fp,
            int order, double sigma, 
            bool fixcen=false, bool fixgam=false, bool fixmu=false,
            bool useflux=false);
        ~EllipseSolver();

        void calculateF(const DVector& x, DVector& f) const;
        void calculateJ(const DVector& x, const DVector& f, DMatrix& df) const;

        void useNumericJ();
        const BVec& getB() const;
        void getBCov(DMatrix& bcov) const;
        void getCovariance(DMatrix& cov) const;
        void getInverseCovariance(DMatrix& invcov) const;

        // CallF takes x and f of length 5, rather than whatever shorter
        // length that F takex (depending on if things are fixed).
        void callF(const DVector& x, DVector& f) const;
        bool solve(DVector& x, DVector& f) const;
        bool testJ(const DVector& x, DVector& f,
                   std::ostream* os=0, double relerr=0.) const;
#ifdef USE_TMV
        void calculateNumericH(
            const DVector& x, const DVector& f, DSymMatrix& h) const;
#endif

    private :

        struct ESImpl;

        ESImpl* _pimpl;
    };

    // Use the integration method, rather than least-squares, to find b.
    class EllipseSolver2 : public BaseEllipseSolver
    {
    public :

        EllipseSolver2(
            const std::vector<PixelList>& pix,
            int order, double sigma, double pixscale,
            bool fixcen=false, bool fixgam=false, bool fixmu=false,
            bool useflux=false);
        EllipseSolver2(
            const std::vector<PixelList>& pix,
            const std::vector<BVec>& psf, double fp,
            int order, double sigma, double pixscale,
            bool fixcen=false, bool fixgam=false, bool fixmu=false,
            bool useflux=false);
        ~EllipseSolver2();

        void calculateF(const DVector& x, DVector& f) const;
        void calculateJ(const DVector& x, const DVector& f, DMatrix& df) const;

        void callF(const DVector& x, DVector& f) const;
        bool solve(DVector& x, DVector& f) const;
        bool testJ(const DVector& x, DVector& f,
                   std::ostream* os=0, double relerr=0.) const;

        void useNumericJ();
        const BVec& getB() const;
        void getBCov(DMatrix& bcov) const;

    private :

        struct ESImpl2;

        ESImpl2* _pimpl;
    };

}}}}

#endif
