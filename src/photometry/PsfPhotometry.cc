// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief A class that knows how to calculate fluxes using the PSF photometry algorithm
 * @ingroup meas/algorithms
 */
template<typename ExposureT>
class PsfPhotometer : public Algorithm<afwDet::Photometry, ExposureT>
{
public:
    typedef Algorithm<afwDet::Photometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<PsfPhotometer> Ptr;
    typedef boost::shared_ptr<PsfPhotometer const> ConstPtr;

    /// Ctor
    PsfPhotometer() : AlgorithmT() {}

    virtual std::string getName() const { return "PSF"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<PsfPhotometer<ExposureT> >();
    }

    virtual PTR(afwDet::Photometry) measureNull(void) const {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<afwDet::Photometry>(NaN, NaN);
    }

    virtual void configure(lsst::pex::policy::Policy const&) {}
    virtual PTR(afwDet::Photometry) measureOne(ExposurePatch<ExposureT> const&, afwDet::Source const&) const;
};

namespace {
/**
 * Accumulate sum(x) and sum(x**2)
 */
template<typename T>
struct getSum2 {
    getSum2() : sum(0.0), sum2(0.0) {}
    
    getSum2& operator+(T x) {
        sum += x;
        sum2 += x*x;
        
        return *this;
    }
    
    double sum;                         // \sum_i(x_i)
    double sum2;                        // \sum_i(x_i^2)
};

template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage, ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _sumVar(0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset() {}        
    void reset(afwDetection::Footprint const& foot) {
        _sumVar = _sum = 0.0;

        afwGeom::BoxI const& bbox(foot.getBBox());
        _x0 = bbox.getMinX();
        _y0 = bbox.getMinY();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size "
                                             "for %d x %d weight image") %
                               bbox.getMinX() % bbox.getMinY() % bbox.getMaxX() % bbox.getMaxY() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    
    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        typename MaskedImageT::Variance::Pixel vval = iloc.variance(0, 0);
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
        _sumVar += wval*wval*vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

    /// Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }
private:
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    double _sumVar;
    int _x0, _y0;                                     // the origin of the current Footprint
};

}
    
/************************************************************************************************************/
/**
 * Calculate the desired psf flux
 */
template<typename ExposureT>
PTR(afwDet::Photometry) PsfPhotometer<ExposureT>::measureOne(ExposurePatch<ExposureT> const& patch,
                                                             afwDet::Source const& source) const
{
    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;

    CONST_PTR(ExposureT) exposure = patch.getExposure();
    CONST_PTR(afwDet::Peak) peak = patch.getPeak();
    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = peak->getFx();   ///< object's column position
    double const ycen = peak->getFy();   ///< object's row position
    
    // BBox for data image
    afwGeom::BoxI imageBBox(mimage.getBBox(afwImage::PARENT));
    
    afwDetection::Psf::ConstPtr psf = exposure->getPsf();
    if (!psf) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, "No PSF provided for PSF photometry");
    }

    afwDetection::Psf::Image::Ptr wimage;

    try {
        wimage = psf->computeImage(afwGeom::PointD(xcen, ycen));
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Computing PSF at (%.3f, %.3f)") % xcen % ycen).str());
        throw e;
    }
        
    FootprintWeightFlux<MaskedImageT, afwDetection::Psf::Image> wfluxFunctor(mimage, wimage);
    // Build a rectangular Footprint corresponding to wimage
    afwDetection::Footprint foot(wimage->getBBox(afwImage::PARENT), imageBBox);
    wfluxFunctor.apply(foot);
    
    getSum2<afwDetection::Psf::Pixel> sum;
    sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);
    
    double flux = wfluxFunctor.getSum()*sum.sum/sum.sum2;
    double fluxErr = ::sqrt(wfluxFunctor.getSumVar())*::fabs(sum.sum)/sum.sum2;
    return boost::make_shared<afwDet::Photometry>(flux, fluxErr);
}

// Declare the existence of a "PSF" algorithm to MeasurePhotometry
DECLARE_ALGORITHM(PsfPhotometer, afwDet::Photometry);

}}}
