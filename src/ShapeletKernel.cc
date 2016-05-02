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
#include "assert.h"

#include "lsst/meas/algorithms/ShapeletKernel.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"
#include "lsst/afw/math/Function.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletKernelFunction : 
        public lsst::afw::math::Function2<double>
    {
    public :
        typedef lsst::afw::math::Function2<double> base;

        typedef std::shared_ptr<ShapeletKernelFunction> Ptr;
        typedef std::shared_ptr<const ShapeletKernelFunction> ConstPtr;

        ShapeletKernelFunction(Shapelet::ConstPtr shapelet) :
            base(shapelet->size()+1),
            _order(shapelet->getOrder()), _size(shapelet->size())
        {
            const Shapelet::ShapeletVector& vals = shapelet->getValues();
            std::vector<double> v(_size+1);
            for(int i=0;i<_size;++i) v[i] = vals(i);
            v[_size] = shapelet->getSigma();
            base::setParameters(v);
        }

        ShapeletKernelFunction(int order, int size) :
            base(size+1), _order(order), _size(size)
        {}

        base::Ptr clone() const
        { return base::Ptr(new ShapeletKernelFunction(*this)); }

        double operator()(double x, double y) const
        {
            Shapelet::ShapeletVector vals(_size);
            const std::vector<double>& v = base::getParameters();
            for(int i=0;i<_size;++i) vals(i) = v[i];
            double sigma = v[_size];
            Shapelet shapelet(_order,sigma,vals);
            return shapelet.evaluateAt(x,y);
        }

    private :
        const int _order;
        const int _size;
    };

    class ShapeletSpatialFunction : 
        public lsst::afw::math::Function2<double>
    {
    public :
        typedef lsst::afw::math::Function2<double> base;

        typedef std::shared_ptr<ShapeletSpatialFunction> Ptr;
        typedef std::shared_ptr<const ShapeletSpatialFunction> ConstPtr;

        ShapeletSpatialFunction(
            ShapeletInterpolation::ConstPtr interp, int i
        ) :
            base(interp->getFitSize()), _interp(interp), _i(i)
        {}

        base::Ptr clone() const
        { return base::Ptr(new ShapeletSpatialFunction(*this)); }

        double operator()(double x, double y) const
        { return _interp->interpolateSingleElement(x,y,_i); }

    private :
        const ShapeletInterpolation::ConstPtr _interp;
        const int _i;

    };

    // Use a radius of 5 sigma for width,height if they are not specified.
    template <class ShapeletPtr>
    inline int getImageSize(ShapeletPtr shapelet, int size)
    { return size == 0 ? int(ceil(shapelet->getSigma()*10.)) : size; }

    LocalShapeletKernel::LocalShapeletKernel(
        Shapelet::ConstPtr shapelet, const Wcs::ConstPtr& wcsPtr, const Extent& size) :
        base( getImageSize(shapelet,size.getX()),
              getImageSize(shapelet,size.getY()),
              ShapeletKernelFunction(shapelet) ),
        _shapelet(shapelet), _wcsPtr(wcsPtr->clone())
    {}

    LocalShapeletKernel::LocalShapeletKernel(
        Shapelet::ConstPtr shapelet, const Wcs::ConstPtr& wcsPtr) :
        base( getImageSize(shapelet,0),
              getImageSize(shapelet,0),
              ShapeletKernelFunction(shapelet) ),
        _shapelet(shapelet), _wcsPtr(wcsPtr->clone())
    {}

    double LocalShapeletKernel::computeImage(
        Image& image, bool doNormalize, double /*x*/, double /*y*/) const
    {
        //TODO This isn't quite right.  I need to account for the WCS.
        using shapelet::CDVector;
        using shapelet::DMatrix;
        using shapelet::DVector;
        using shapelet::makePsi;

        const int xCen = base::getCtrX();
        const int yCen = base::getCtrY();
        const int nX = base::getWidth();
        const int nY = base::getHeight();
        const int nPixels = nX * nY;
        const double sigma = _shapelet->getSigma();
        CDVector zList(nPixels);
        int k=0;
        for(int i=0;i<nX;++i) {
            for(int j=0;j<nY;++j) {
                std::complex<double> z(i-xCen,j-yCen);
                z /= sigma;
                zList[k++] = z;
            }
        }
        Assert(k == nPixels);
        DMatrix psi(nPixels,_shapelet->size());
        makePsi(psi,TMV_vview(zList),_shapelet->getOrder());
        DVector flux = psi * _shapelet->getValues();

        double sum = flux.TMV_sumElements();
        if (doNormalize) flux /= sum;

        k=0;
        for(int i=0;i<nX;++i) {
            for(int j=0;j<nY;++j) {
                image(i,j) = flux(k++);
            }
        }

        return sum;
    }

    std::vector<lsst::afw::math::Kernel::SpatialFunctionPtr> buildSetOfSpatialFunctions(
        ShapeletInterpolation::ConstPtr interp)
    {
        const int nParam = interp->getFitSize();
        typedef lsst::afw::math::Kernel::SpatialFunctionPtr SFPtr;
        std::vector<SFPtr> ret;
        ret.reserve(nParam);
        for(int i=0;i<nParam;++i) {
            ret.push_back(SFPtr(new ShapeletSpatialFunction(interp,i)));
        }
        return ret;
    }

    ShapeletKernel::ShapeletKernel(
        ShapeletInterpolation::ConstPtr interp, 
        const Wcs::ConstPtr& wcsPtr, const Extent& size
    ) :
        base( getImageSize(interp,size.getX()),
              getImageSize(interp,size.getY()),
              ShapeletKernelFunction(interp->getOrder(),interp->getSize()),
              buildSetOfSpatialFunctions(interp) ),
        _interp(interp), _wcsPtr(wcsPtr->clone())
    {}

    ShapeletKernel::ShapeletKernel(
        ShapeletInterpolation::ConstPtr interp, const Wcs::ConstPtr& wcsPtr
    ) :
        base( getImageSize(interp,0),
              getImageSize(interp,0),
              ShapeletKernelFunction(interp->getOrder(),interp->getSize()),
              buildSetOfSpatialFunctions(interp) ),
        _interp(interp), _wcsPtr(wcsPtr->clone())
    {}

    LocalShapeletKernel::ConstPtr ShapeletKernel::getLocalKernel(
        const Point& pos) const
    { 
        Shapelet::ConstPtr localShapelet(_interp->interpolate(pos));
        LocalShapeletKernel::ConstPtr localShapeletKernel(
            new LocalShapeletKernel(localShapelet, _wcsPtr));
        return localShapeletKernel;
    }

    double ShapeletKernel::computeImage(
        Image& image, bool doNormalize, double x, double y) const
    {
        Point point(x,y);
        return getLocalKernel(point)->computeImage(image,doNormalize);
    }

}}}
