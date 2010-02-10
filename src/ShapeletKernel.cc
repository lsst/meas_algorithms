
#include "lsst/meas/algorithms/ShapeletKernel.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"
#include "lsst/afw/math/Function.h"
#include "assert.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletKernelFunction : 
        public lsst::afw::math::Function2<double>
    {
    public :
        typedef lsst::afw::math::Function2<double> base;

        typedef boost::shared_ptr<ShapeletKernelFunction> Ptr;
        typedef boost::shared_ptr<const ShapeletKernelFunction> ConstPtr;

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

        typedef boost::shared_ptr<ShapeletSpatialFunction> Ptr;
        typedef boost::shared_ptr<const ShapeletSpatialFunction> ConstPtr;

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

    LocalShapeletKernel::LocalShapeletKernel(
        Shapelet::ConstPtr shapelet, Wcs::Ptr wcs) :
        // Use a radius of 5 sigma for width,height:
        base( int(ceil(shapelet->getSigma()*10.)),
              int(ceil(shapelet->getSigma()*10.)),
              ShapeletKernelFunction(shapelet) ),
        _shapelet(shapelet), _wcs(wcs)
    {}

    double LocalShapeletKernel::computeImage(
        Image& image, bool doNormalize, double /*x*/, double /*y*/) const
    {
        //TODO This isn't right.  I need to account for the WCS.
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

        if (doNormalize) flux /= flux.TMV_sumElements();

        k=0;
        for(int i=0;i<nX;++i) {
            for(int j=0;j<nY;++j) {
                image(i,j) = flux(k++);
            }
        }

        // FIXME: This is dumb.  The only thing the return value is used
        // for is to test that the image was normalized correctly.
        return flux.TMV_sumElements();
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
        ShapeletInterpolation::ConstPtr interp, Wcs::Ptr wcs
    ) :
        base( interp->getSigma()*5.,
              interp->getSigma()*5.,
              ShapeletKernelFunction(interp->getOrder(),interp->getSize()),
              buildSetOfSpatialFunctions(interp) ),
        _interp(interp), _wcs(wcs)
    {}

    LocalShapeletKernel::ConstPtr ShapeletKernel::getLocalKernel(
        const PointD& pos) const
    { 
        Shapelet::ConstPtr localShapelet(_interp->interpolate(pos));
        LocalShapeletKernel::ConstPtr localShapeletKernel(
            new LocalShapeletKernel(localShapelet,_wcs));
        return localShapeletKernel;
    }

    double ShapeletKernel::computeImage(
        Image& image, bool doNormalize, double x, double y) const
    {
        PointD point(x,y);
        return getLocalKernel(point)->computeImage(image,doNormalize);
    }

}}}
