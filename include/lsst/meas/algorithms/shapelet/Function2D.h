#ifndef MeasAlgoShapeletFunction2D_H
#define MeasAlgoShapeletFunction2D_H

#include <iostream>
#include <functional>
#include <memory>
#include "lsst/meas/algorithms/shapelet/Bounds.h"
#include "lsst/meas/algorithms/shapelet/MyMatrix.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {


    struct RangeException : public std::runtime_error
    {
    public :
        RangeException(const Position& p, const Bounds& b);
        ~RangeException() throw() {} 

        const Bounds& getBounds() const { return _b; }
        const Position& getPosition() const { return _p; }

    private : 
        Position _p;
        Bounds _b;

    };

    class Function2D 
    {
    public:

        Function2D() : _xOrder(0), _yOrder(0), _coeffs(new DMatrix(1,1)) 
        { (*_coeffs)(0,0) = 0.0; }

        Function2D(int xo, int yo) : 
            _xOrder(xo), _yOrder(yo), _coeffs(new DMatrix(xo+1,yo+1)) 
        { _coeffs->setZero(); }

        Function2D(int xo, int yo, const DMatrix& c) :
            _xOrder(xo),_yOrder(yo),_coeffs(new DMatrix(c)) {}

        Function2D(const Function2D& rhs) :
            _xOrder(rhs._xOrder), _yOrder(rhs._yOrder),
            _coeffs(new DMatrix(*rhs._coeffs)) {}

        virtual ~Function2D() {}

        virtual void write(std::ostream& fout) const =0;

        static std::auto_ptr<Function2D> read(std::istream& fin);

        virtual std::auto_ptr<Function2D> copy() const =0;

        // Adds a + bx + cy to function
        virtual void addLinear(double a, double b, double c) = 0;

        // Converts function to f(a+bx+cy,d+ex+fy)
        virtual void linearPreTransform(
            double a, double b, double c, double d, double e, double f) = 0;

        // returns new function x derivative of f
        virtual std::auto_ptr<Function2D> dFdX() const = 0;

        // returns new function y derivative of f
        virtual std::auto_ptr<Function2D> dFdY() const = 0;

        virtual std::auto_ptr<Function2D> conj() const;

        virtual void operator*=(double scale)
        { *_coeffs *= scale; }

        virtual double operator()(double x,double y) const;

        double operator()(const Position& p) const 
        { return operator()(p.getX(),p.getY()); }

        virtual void setTo(double value) 
        {
            if (_xOrder || _yOrder) {
                _xOrder = 0; _yOrder = 0; 
                _coeffs.reset(new DMatrix(1,1));
            }
            (*_coeffs)(0,0) = value;
        }

        bool isNonZero() const 
        { return (*_coeffs)(0,0) != 0.0 || _xOrder!=0 || _yOrder!=0; }

        int getXOrder() const { return _xOrder; }

        int getYOrder() const { return _yOrder; }

        const DMatrix& getCoeffs() const { return *_coeffs; }

        // Sets function to fit of f(pos_i) = v_i using i only if 
        // shouldUse[i] = true
        virtual void simpleFit(
            int order, const std::vector<Position>& pos, 
            const std::vector<double>& v, const std::vector<bool>& shouldUse,
            const std::vector<double>* sigList=0,
            double* chisqOut = 0, int* dofOut=0, DMatrix* cov=0);

        // Sets function to fit of f(pos_i) = v_i using i if fit is within 
        // nsig sigma.  *shouldUse is returned as list of good i's.
        virtual void outlierFit(
            int order,double nsig,
            const std::vector<Position>& pos, const std::vector<double>& v,
            std::vector<bool>* shouldUse,
            const std::vector<double>* sigList=0, 
            double* chisqout = 0, int* dofout = 0, DMatrix* cov=0);

        // Sets function to fit of f(pos_i) = v_i reducing the order as far
        // as possible keeping quality of fit the same as for maxOrder
        // equivProb = rejection percentile.  eg. 0.9 means a low order fit is
        // rejected if it is statistically rejected at the 90th percentile
        // Higher values of equivProb result in lower order fits.
        virtual void orderFit(
            int maxOrder, double equivProb,
            const std::vector<Position>& pos, const std::vector<double>& v,
            const std::vector<bool>& shouldUse, const std::vector<double>* sigList=0, 
            double *chisqout = 0, int *dofout = 0, DMatrix* cov=0);

        // Sets function to h(x,y) = a + b*f(x,y) + c*g(x,y)
        virtual void linearTransform(
            double a, double b, double c, 
            const Function2D& f, const Function2D& g);

        // Sets function to g(x,y) = a + b*f(x,y)
        virtual void linearTransform(double a, double b, const Function2D& f)
        { linearTransform(a,b,0.,f,f); }

        virtual void operator+=(const Function2D& rhs) = 0;

        virtual void setFunction(
            int _xOrder, int _yOrder, const DVector& fVect) = 0;

    protected:

        virtual DVector definePX(int order, double x) const = 0;
        virtual DVector definePY(int order, double y) const = 0;

        int _xOrder,_yOrder;
        std::auto_ptr<DMatrix> _coeffs;

    private:

        void doSimpleFit(
            int xOrder, int yOrder, 
            const std::vector<Position>& pos, const std::vector<double>& v,
            const std::vector<bool>& shouldUse, DVector *f, 
            const std::vector<double>* sigList=0, int *dof=0,
            DVector *diff=0, DMatrix* cov=0);
    };

    inline std::ostream& operator<<(std::ostream& fout, const Function2D& f)
    { f.write(fout); return fout; }

    inline std::istream& operator>>(
        std::istream& fin, std::auto_ptr<Function2D>& f)
    { f = Function2D::read(fin); return fin; }

    class Constant2D : public Function2D 
    {

    public:

        Constant2D() : Function2D() {}

        Constant2D(const Constant2D& rhs) : Function2D(rhs) {}

        Constant2D(double value) : Function2D(0,0)
        { Function2D::setTo(value); }

        Constant2D(std::istream& fin);

        virtual ~Constant2D() {}

        virtual double operator()(double ,double ) const 
        { return (*this->_coeffs)(0,0); }

        virtual void write(std::ostream& fout) const;

        virtual std::auto_ptr<Function2D> dFdX() const 
        { return std::auto_ptr<Function2D>(new Constant2D()); }

        virtual std::auto_ptr<Function2D> dFdY() const 
        { return std::auto_ptr<Function2D>(new Constant2D()); }

        std::auto_ptr<Function2D> copy() const 
        { return std::auto_ptr<Function2D>(new Constant2D(*this)); }

        virtual void addLinear(double a, double b, double c) 
        { 
            Assert(b == 0.);
            Assert(c == 0.);
            (*this->_coeffs)(0,0) += a; 
        }

        virtual void linearPreTransform(
            double , double , double , double , double , double )
        { Assert(false); return; }

        virtual void operator+=(const Function2D& rhs);

        virtual void setFunction(
            int /*xOrder*/, int /*yOrder*/, const DVector& fVect)
        { 
            Assert(_xOrder == 0);
            Assert(_yOrder == 0);
            Assert(fVect.size()==1);
            (*this->_coeffs)(0,0) = fVect(0); 
        }

    private:

        virtual DVector definePX(int order, double ) const
        { 
            Assert(order == 0);
            DVector temp(1);
            temp(0) = 1.0;
            return temp;
        }
        virtual DVector definePY(int order, double y) const
        { return definePX(order,y); }

        using Function2D::_coeffs;
        using Function2D::_xOrder;
        using Function2D::_yOrder;
    };

    class Polynomial2D : public Function2D 
    {

    public:

        Polynomial2D(double scale=1.) : Function2D(),_scale(scale) {}

        Polynomial2D(const Polynomial2D& rhs) :
            Function2D(rhs),_scale(rhs._scale) {}

        Polynomial2D(const DMatrix& a, double scale=1.) : 
            Function2D(a.TMV_colsize()-1,a.TMV_rowsize()-1,a), _scale(scale) {}

        Polynomial2D(std::istream& fin);

        Polynomial2D(int xo,int yo,double scale=1.) :
            Function2D(xo,yo), _scale(scale) {}

        virtual ~Polynomial2D() {}

        virtual void write(std::ostream& fout) const;

        virtual std::auto_ptr<Function2D> dFdX() const;

        virtual std::auto_ptr<Function2D> dFdY() const;

        std::auto_ptr<Function2D> copy() const 
        { return std::auto_ptr<Function2D>(new Polynomial2D(*this)); }

        virtual void addLinear(double a, double b, double c);

        virtual void linearPreTransform(
            double a, double b, double c, double d, double e, double f);

        virtual void operator+=(const Function2D& rhs);

        void makeProductOf(const Polynomial2D& f, const Polynomial2D& g);

        virtual void setFunction(
            int xOrder, int yOrder, const DVector& fVect);

    private:

        double _scale;

        virtual DVector definePX(int order, double x) const
        {
            DVector temp(order+1);
            temp(0) = 1.0;
            for(int i=1;i<=order;++i) temp(i) = temp(i-1)*x/_scale;
            return temp;
        }
        virtual DVector definePY(int order, double y) const
        { return definePX(order,y); }

        using Function2D::_coeffs;
        using Function2D::_xOrder;
        using Function2D::_yOrder;
    };

}}}}

#endif
