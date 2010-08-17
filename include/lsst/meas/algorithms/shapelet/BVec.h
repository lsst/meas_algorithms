#ifndef MeasAlgoShapeletBVec_H
#define MeasAlgoShapeletBVec_H

#include "lsst/meas/algorithms/shapelet/MyMatrix.h"
#include <complex>
#include <vector>
#include "lsst/meas/algorithms/shapelet/dbg.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class BVec;

    class AssignableToBVec
    {
    public:
        virtual void AssignTo(BVec& bvec) const = 0;
        virtual int getOrder() const = 0;
        virtual double getSigma() const = 0;
        virtual ~AssignableToBVec() {}
    };

    class BVec : public AssignableToBVec
    {

    public :

        BVec(int order, double sigma) :
            _order(order), _sigma(sigma),
            _b((_order+1)*(_order+2)/2) 
        { _b.setZero(); }

#ifdef USE_TMV
        BVec(int order, double sigma, double* bvec) :
            _order(order), _sigma(sigma),
            _b((_order+1)*(_order+2)/2,bvec) 
        {}
#else
        BVec(int order, double sigma, double* bvec) :
            _order(order), _sigma(sigma),
            _b(DVector::Map(bvec,(_order+1)*(_order+1)/2))
        {} 
#endif

        BVec(int order, double sigma, const DVector& bvec) :
            _order(order), _sigma(sigma), _b(bvec)
        {}

        BVec(const BVec& rhs) :
            _order(rhs._order), _sigma(rhs._sigma), _b(rhs._b)
        {}

        virtual ~BVec() {}

        BVec& operator=(const AssignableToBVec& rhs);
        BVec& operator=(const BVec& rhs);

        void AssignTo(BVec& bvec) const;

        DVector& vec() { return _b; }
        const DVector& vec() const { return _b; }
        double& operator()(int i) { return _b(i); }
        double operator()(int i) const { return _b(i); }

        int getOrder() const { return _order; }
        double getSigma() const { return _sigma; }
        const DVector& getValues() const { return _b; }
        int size() const { return _b.size(); }

        void setSigma(double sigma) { _sigma = sigma; }

        // For this one v is allowed to be smaller than size, 
        // and the rest of the values are filled in as 0's.
        void setValues(const DVector& v);

        // For this one, the sizes must match.
        // b.setValues(v) is equivalent to b.vec() = v.
        template <typename V>
        void setValues(const V& v) 
        { 
            Assert(v.size() == size());
            _b = v; 
        }

        void normalize() { _b /= _b(0); }

    private :

        int _order;
        double _sigma;
        DVector _b;

    };

    void calculateZTransform(
        std::complex<double> z, int order, DMatrix& T);
    void applyZ(std::complex<double> z, BVec& b);

    void calculateMuTransform(double mu, int order, DMatrix& D);
    void augmentMuTransformRows(double mu, int order, DMatrix& D);
    void applyMu(double mu, BVec& b);

    void calculateGTransform(std::complex<double> g, int order, DMatrix& S);
    void augmentGTransformCols(std::complex<double> g, int order, DMatrix& S);
    void applyG(std::complex<double> g, BVec& b);

    void calculatePsfConvolve(
        const BVec& bpsf, int order, double sigma, DMatrix& C);
    void applyPsf(const BVec& bpsf, BVec& b);

}}}} // namespace shapelet

#endif
