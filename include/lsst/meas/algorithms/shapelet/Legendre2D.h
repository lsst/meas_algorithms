#ifndef MeasAlgoShapeletLegendre2D_H
#define MeasAlgoShapeletLegendre2D_H

#include <iostream>
#include "lsst/meas/algorithms/shapelet/Function2D.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class Legendre2D : public Function2D
    {

    public:

        Legendre2D() {}

        Legendre2D(double xMin,double xMax,double yMin,double yMax) :
            Function2D(), _bounds(xMin,xMax,yMin,yMax) {}

        Legendre2D(const Bounds& b) : Function2D(), _bounds(b) {}

        Legendre2D(const Legendre2D& rhs) :
            Function2D(rhs), _bounds(rhs._bounds) {}

        Legendre2D(const Bounds& b, const DMatrix& a) : 
            Function2D(a.TMV_colsize()-1,a.TMV_rowsize()-1,a), _bounds(b) {}

        Legendre2D(int xo, int yo, const Bounds& b) :
            Function2D(xo,yo), _bounds(b) {}

        Legendre2D(std::istream& fin);

        virtual ~Legendre2D() {}

        virtual void write(std::ostream& fout) const;

        virtual std::auto_ptr<Function2D> dFdX() const;

        virtual std::auto_ptr<Function2D> dFdY() const;

        virtual std::auto_ptr<Function2D> copy() const
        { return std::auto_ptr<Function2D>(new Legendre2D(*this)); }

        virtual void addLinear(double a, double b, double c);

        virtual void linearPreTransform(
            double a, double b, double c, double d, double e, double f);

        virtual void operator+=(const Function2D& rhs);

        double getXMin() const {return _bounds.getXMin();}

        double getXMax() const {return _bounds.getXMax();}

        double getYMin() const {return _bounds.getYMin();}

        double getYMax() const {return _bounds.getYMax();}

        const Bounds& getBounds() const {return _bounds;}

        virtual void setFunction(
            int xorder, int yorder, const DVector& fvect);

    private:

        Bounds _bounds;
        using Function2D::_xOrder;
        using Function2D::_yOrder;
        using Function2D::_coeffs;

        DVector definePXY(
            int order,double xy,double min,double max) const;

        virtual DVector definePX(int order, double x) const;

        virtual DVector definePY(int order, double y) const;

    };

}}}}

#endif
