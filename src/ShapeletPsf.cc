
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
        typedef lsst::afw::image::MaskedImage<double> MaskedImage;
        typedef lsst::afw::image::Wcs Wcs;

        ShapeletPsfVisitor(
            int order, double sigma, double aperture,
            const MaskedImage& image, const Wcs& wcs
        ) :
            _order(order), _sigma(sigma), _aperture(aperture),
            _image(image), _wcs(wcs)
        {}

        void reset() {}

        void processCandidate(SpatialCellCandidate* cand) 
        {
            using lsst::afw::detection::Source;
            using lsst::afw::geom::PointD;

            ShapeletPsfCandidate* psfCand = 
                dynamic_cast<ShapeletPsfCandidate*>(cand);
            Assert(psfCand);
            Shapelet::Ptr shape(new Shapelet(_order,_sigma));
            const Source& source = *(psfCand->getSource());
            PointD pos = lsst::afw::geom::makePointD(psfCand->getX(),psfCand->getY());

            // Convert the aperture to pixels.
            // pixelScale is arcsec/pixel
            double pixelScale = sqrt(getJacobian(_wcs,pos).determinant());
            double pixelAperture = _aperture / pixelScale;

            if (!shape->measureFromImage(
                    source,pos,false,true,pixelAperture,_image,_wcs)) {
                psfCand->setBad();
            }
            psfCand->setShapelet(shape);
        }

    private :
        int _order;
        double _sigma;
        double _aperture;
        const MaskedImage& _image;
        const Wcs& _wcs;
    };

    class ShapeletPsfImpl 
    {
    public :
        typedef ShapeletPsf::Policy Policy;
        typedef ShapeletPsf::SpatialCellSet SpatialCellSet;
        typedef ShapeletPsf::MaskedImage MaskedImage;
        typedef ShapeletPsf::Wcs Wcs;
        typedef ShapeletPsf::Point Point;
        typedef ShapeletPsf::Extent Extent;
        typedef ShapeletPsf::Color Color;
        typedef ShapeletPsf::Kernel Kernel;
        typedef ShapeletPsf::Ptr Ptr;
        typedef ShapeletPsf::ConstPtr ConstPtr;

        ShapeletPsfImpl(
            const Policy& policy,
            const SpatialCellSet& cellSet,
            const MaskedImage& image,
            const Wcs& wcs
        ) : 
            _cellSet(new SpatialCellSet(cellSet)), 
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
                order, sigma, aperture, image, _wcs);
            _cellSet->visitCandidates(&visitor2);

            // Resort the Spatial cell, since the ratings have changed.
            //
            // Note: this currently requires the trunk version of afw, which in turn
            // requires the trunk version of daf_base.  So if you're not well-versed in 
            // eups and setup, here is the procedure to make this work:
            //
            // Start in the svn root directory.  Then:
            // 
            // cd daf/base/trunk
            // setup -r .
            // scons
            // cd ../../../afw/trunk
            // setup -r .
            // setup -r ../../daf/base/trunk
            // scons
            // cd ../../meas/algorithms/tickets/1158
            // setup -r .
            // setup -r ../../../../daf/base/trunk
            // setup -r ../../../../afw/trunk
            // scons -c
            // scons
            //
            _cellSet->sortCandidates();
            
            // Finally do the interpolation with a FittedShapelet object.
            _interp->calculate(_cellSet,image,wcs);
        }

        // Default destructor, copy constructor and op= do the right thing.
        
        Kernel::ConstPtr getLocalKernel(
            const Color& color, const Point& pos) const
        {
            // TODO: For now we ignore the color argument.
            // This functionality needs to be added!
            return LocalShapeletKernel::ConstPtr(
                new LocalShapeletKernel(_interp->interpolate(pos),_wcs));
        }

        Kernel::Ptr getLocalKernel(const Color& color, const Point& pos)
        {
            return LocalShapeletKernel::Ptr(
                new LocalShapeletKernel(_interp->interpolate(pos),_wcs));
        }

        Kernel::ConstPtr getKernel(const Color& color) const
        {
            return ShapeletKernel::ConstPtr(
                new ShapeletKernel(_interp,_wcs)); 
        }

        Kernel::Ptr getKernel(const Color& color)
        {
            return ShapeletKernel::Ptr(
                new ShapeletKernel(_interp,_wcs)); 
        }

        const SpatialCellSet& getCellSet() const
        { return *_cellSet; }

    private :
        SpatialCellSet::Ptr _cellSet;
        ShapeletInterpolation::Ptr _interp;
        const Wcs& _wcs;
    };

    ShapeletPsf::ShapeletPsf(
        const Policy& policy,
        const SpatialCellSet& cellSet,
        const MaskedImage& image,
        const Wcs& wcs
    ) :
        pImpl(new ShapeletPsfImpl(policy,cellSet,image,wcs))
    {}

    ShapeletPsf::~ShapeletPsf()
    { delete pImpl; pImpl = 0; }

    ShapeletPsf::ShapeletPsf(const ShapeletPsf& rhs)
    { 
        if(pImpl) delete pImpl; 
        pImpl = new ShapeletPsfImpl(*rhs.pImpl);
    }

    ShapeletPsf::Kernel::ConstPtr ShapeletPsf::doGetLocalKernel(
        const ShapeletPsf::Point& pos,
        const ShapeletPsf::Color& color
    ) const
    { return pImpl->getLocalKernel(color,pos); }

    ShapeletPsf::Kernel::Ptr ShapeletPsf::doGetLocalKernel(
        const ShapeletPsf::Point& pos,
        const ShapeletPsf::Color& color
    )
    { return pImpl->getLocalKernel(color,pos); }

    ShapeletPsf::Kernel::ConstPtr ShapeletPsf::doGetKernel(
        const ShapeletPsf::Color& color
    ) const
    { return pImpl->getKernel(color); }

    ShapeletPsf::Kernel::Ptr ShapeletPsf::doGetKernel(
        const ShapeletPsf::Color& color
    )
    { return pImpl->getKernel(color); }

    const lsst::afw::math::SpatialCellSet& ShapeletPsf::getCellSet() const
    { return pImpl->getCellSet(); }


}}}
