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
        virtual void assignTo(BVec& bvec) const = 0;
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

        BVec(int order, double sigma, double* bvec) :
            _order(order), _sigma(sigma),
#ifdef USE_TMV
            _b((_order+1)*(_order+2)/2,bvec) 
#else
                _b(DVector::Map(bvec,(_order+1)*(_order+1)/2))
#endif
                {} 

        BVec(int order, double sigma, const DVector& bvec) :
            _order(order), _sigma(sigma), _b(bvec)
        {}

        BVec(const BVec& rhs) :
            _order(rhs._order), _sigma(rhs._sigma), _b(rhs._b)
        {}

        ~BVec() {}

        BVec& operator=(const AssignableToBVec& rhs);
        BVec& operator=(const BVec& rhs);

        void assignTo(BVec& bvec) const;

        DVector& vec() { return _b; }
        const DVector& vec() const { return _b; }
        double& operator()(int i) { return _b(i); }
        double operator()(int i) const { return _b(i); }

        int getOrder() const { return _order; }
        double getSigma() const { return _sigma; }
        const DVector& getValues() const { return _b; }
        size_t size() const { return _b.size(); }

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

        void conjugateSelf();

    private :

        int _order;
        double _sigma;
        DVector _b;

    };

    inline std::ostream& operator<<(std::ostream& os, const BVec& b)
    { os << b.getOrder()<<"  "<<b.getSigma()<<"  "<<b.vec(); return os; }

    // NB: All of the following calculate and augment functions assume that
    // the input matrix has already been zeroed before calling the function.
    // This is because the sparsity of the matrices maintains its structure
    // for different values of mu, g, theta, and z.  So it is faster to 
    // just overwrite the locations that need to be written and skip the zeroes
    // each time you call these functions.
    void calculateZTransform(
        std::complex<double> z, int order1, int order2, DMatrix& T);
    inline void calculateZTransform(std::complex<double> z, int order, DMatrix& T)
    { calculateZTransform(z,order,order,T); }
    void augmentZTransformCols(std::complex<double> z, int order, DMatrix& T);
    void applyZ(std::complex<double> z, BVec& b);

    void calculateMuTransform(double mu, int order1, int order2, DMatrix& D);
    inline void calculateMuTransform(double mu, int order, DMatrix& D)
    { calculateMuTransform(mu,order,order,D); }
    void augmentMuTransformRows(double mu, int order, DMatrix& D);
    void augmentMuTransformCols(double mu, int order, DMatrix& D);
    void applyMu(double mu, BVec& b);

    void calculateThetaTransform(
        double theta, int order1, int order2, DBandMatrix& R);
    inline void calculateThetaTransform(double theta, int order, DBandMatrix& R)
    { calculateThetaTransform(theta,order,order,R); }
    void applyTheta(double theta, BVec& b);

    void calculateGTransform(
        std::complex<double> g, int order1, int order2, DMatrix& S);
    inline void calculateGTransform(std::complex<double> g, int order, DMatrix& S)
    { calculateGTransform(g,order,order,S); }
    void augmentGTransformCols(std::complex<double> g, int order, DMatrix& S);
    void applyG(std::complex<double> g, BVec& b);

    void calculatePsfConvolve(
        const BVec& bpsf, int order1, int order2, double sigma, DMatrix& C);
    inline void calculatePsfConvolve(
        const BVec& bpsf, int order, double sigma, DMatrix& C)
    { calculatePsfConvolve(bpsf,order,order,sigma,C); }
    void applyPsf(const BVec& bpsf, BVec& b);

}}}} 

#endif
