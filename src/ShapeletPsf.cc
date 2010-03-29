
#include "lsst/meas/algorithms/ShapeletPsf.h"
#include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
#include "lsst/meas/algorithms/ShapeletInterpolation.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
#include "assert.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class MeanSizeVisitor : 
        public lsst::afw::math::CandidateVisitor 
    {
    public :
        typedef lsst::afw::math::SpatialCellCandidate SpatialCellCandidate;

        MeanSizeVisitor() : _sum(0.0), _n(0) {}

        void reset() { _sum = 0.0; _n = 0; }

        void processCandidate(SpatialCellCandidate* cand) 
        {
            ShapeletPsfCandidate* psfCand = 
                dynamic_cast<ShapeletPsfCandidate*>(cand);
            Assert(psfCand);
            _sum += psfCand->getSize();
            ++_n;
        }

        double getMean() { return _sum / _n; }

    private :
        double _sum;
        int _n;
    };

    class ShapeletPsfVisitor : 
        public lsst::afw::math::CandidateVisitor 
    {
    public :
        typedef lsst::afw::math::SpatialCellCandidate SpatialCellCandidate;
        typedef lsst::afw::image::Image<double> Image;
        typedef lsst::afw::image::Wcs Wcs;

        ShapeletPsfVisitor(
            int order, double sigma, double aperture,
            Image::ConstPtr image, Wcs::Ptr wcs, Image::ConstPtr weightImage
        ) :
            _order(order), _sigma(sigma), _aperture(aperture),
            _image(image), _wcs(wcs), _weightImage(weightImage)
        {}

        void reset() {}

        void processCandidate(SpatialCellCandidate* cand) 
        {
            using lsst::afw::detection::Source;
            using lsst::afw::image::PointD;

            ShapeletPsfCandidate* psfCand = 
                dynamic_cast<ShapeletPsfCandidate*>(cand);
            Assert(psfCand);
            Shapelet::Ptr shape(new Shapelet(_order,_sigma));
            const Source& source = *(psfCand->getSource());
            PointD pos(psfCand->getX(),psfCand->getY());

            // Convert the aperture to pixels.
            // pixelScale is arcsec/pixel
            double pixelScale = sqrt(getJacobian(*_wcs,pos).determinant());
            double pixelAperture = _aperture / pixelScale;

            if (!shape->measureFromImage(
                    source,pos,false,true,pixelAperture,
                    _image,_wcs,_weightImage)) {
                psfCand->setBad();
            }
            psfCand->setShapelet(shape);
        }

    private :
        int _order;
        double _sigma;
        double _aperture;
        Image::ConstPtr _image;
        Wcs::Ptr _wcs;
        Image::ConstPtr _weightImage;
    };

    class ShapeletPsfImpl 
    {
    public :
        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::math::SpatialCellSet SpatialCellSet;
        typedef lsst::afw::image::Image<double> Image;
        typedef lsst::afw::image::Wcs Wcs;
        typedef lsst::afw::image::PointD PointD;

        ShapeletPsfImpl(
            const Policy& policy,
            SpatialCellSet::Ptr cellSet,
            Image::ConstPtr image,
            Wcs::Ptr wcs,
            Image::ConstPtr weightImage
        ) : 
            _cellSet(cellSet), 
            _interp(new ShapeletInterpolation(policy)),
            _wcs(wcs)
        {
            const int order = policy.getInt("shapeletOrder");

            // Note: This aperture is in arcsec.  Will need to convert to
            // pixels for each star.
            const double aperture = policy.getInt("psfAperture");
             
            // First find the mean size.
            // WARNING: if we stop using the shapelet sigma for the
            // size measurement in the StarFinder, then we should 
            // add a step here to measure sigma for each star before
            // taking the mean.
            MeanSizeVisitor visitor1;
            _cellSet->visitCandidates(&visitor1);
            double sigma = visitor1.getMean();
            
            // ShapeletPsfVisitor visits each candidate and measures the
            // shapelet decomposition.
            ShapeletPsfVisitor visitor2(
                order, sigma, aperture, image, _wcs, weightImage);
            _cellSet->visitCandidates(&visitor2);

            // Resort the Spatial cell, since the ratings have changed.
            // TODO: This needs the trunk version of afw.  Not working yet..
            //_cellSet->sortCandidates();
            
            // Finally do the interpolation with a FittedShapelet object.
            _interp->calculate(_cellSet,image,wcs,weightImage);
        }

        // Default destructor, copy constructor and op= do the right thing.
        
        LocalShapeletKernel::Ptr getLocalKernel(const PointD& pos, double color, int width, int height)
        { 
            return LocalShapeletKernel::Ptr(
                new LocalShapeletKernel(_interp->interpolate(pos),_wcs,width,height));
        }

        ShapeletKernel::Ptr getKernel(double color, int width, int height)
        { return ShapeletKernel::Ptr(new ShapeletKernel(_interp,_wcs,width,height)); }

        const SpatialCellSet& getCellSet() const
        { return *_cellSet; }

    private :
        SpatialCellSet::Ptr _cellSet;
        ShapeletInterpolation::Ptr _interp;
        Wcs::Ptr _wcs;
    };

    ShapeletPsf::ShapeletPsf(
        const Policy& policy, SpatialCellSet::Ptr cellSet,
        Image::ConstPtr image, Wcs::Ptr wcs, Image::ConstPtr weightImage
    ) :
        pImpl(new ShapeletPsfImpl(policy,cellSet,image,wcs,weightImage))
    {}

    ShapeletPsf::~ShapeletPsf()
    { delete pImpl; pImpl = 0; }

    ShapeletPsf::ShapeletPsf(const ShapeletPsf& rhs)
    { 
        if(pImpl) delete pImpl; 
        pImpl = new ShapeletPsfImpl(*rhs.pImpl);
    }

    ShapeletPsf& ShapeletPsf::operator=(const ShapeletPsf& rhs)
    {
        *pImpl = *rhs.pImpl;
        return *this;
    }

    LocalShapeletKernel::Ptr ShapeletPsf::getLocalKernel(
        const PointD& pos, double color, int width, int height) const
    { return pImpl->getLocalKernel(pos,color,width,height); }

    ShapeletKernel::Ptr ShapeletPsf::getKernel(double color, int width, int height) const
    { return pImpl->getKernel(color,width,height); }

    const lsst::afw::math::SpatialCellSet& ShapeletPsf::getCellSet() const
    { return pImpl->getCellSet(); }

}}}


