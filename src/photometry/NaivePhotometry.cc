// -*- LSST-C++ -*-

#if 0 && defined(__ICC)
#pragma warning (push)
#pragma warning (disable: 21)           // type qualifiers are meaningless in this declaration
#pragma warning disable: 68)            // integer conversion resulted in a change of sign
#pragma warning (disable: 279)          // controlling expression is constant
#pragma warning (disable: 304)          // access control not specified ("public" by default)
#pragma warning (disable: 444)          // destructor for base class ... is not virtual
//#pragma warning (pop)
#endif

#include <limits>
#include <numeric>
#include "Eigen/LU"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace detection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {
/**
 * Implement "Naive" photometry.
 */
class NaivePhotometry : public lsst::afw::detection::Photometry
{
public:
    typedef boost::shared_ptr<NaivePhotometry> Ptr;
    typedef boost::shared_ptr<NaivePhotometry const> ConstPtr;

    /// Ctor
    NaivePhotometry(double flux, float fluxErr=-1) {
        init();                         // This allocates space for fields added by defineSchema
        set<FLUX>(flux);                // ... if you don't, these set calls will fail an assertion
        set<FLUX_ERR>(fluxErr);         // the type of the value must match the schema
    }

    /// Add desired fields to the schema
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Photometry::defineSchema(schema);
    }

    template<typename ImageT>
    static Photometry::Ptr doMeasure(typename ImageT::ConstPtr im, detection::Peak const&);
};

template <typename MaskedImageT>
class FootprintFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = 0.0;
    }
    void reset(detection::Footprint const&) {}        

    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int,                                   ///< column-position of pixel
                    int                                    ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);
        _sum += val;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    double _sum;
};

template <typename MaskedImageT, typename WeightImageT>
class FootprintWeightFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    FootprintWeightFlux(MaskedImageT const& mimage,          ///< The image the source lives in
                        typename WeightImageT::Ptr wimage    ///< The weight image
                       ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                           _wimage(wimage),
                           _sum(0), _x0(0), _y0(0) {}
    
    /// @brief Reset everything for a new Footprint
    void reset(detection::Footprint const& foot) {
        _sum = 0.0;
        
        afwImage::BBox const& bbox(foot.getBBox());
        _x0 = bbox.getX0();
        _y0 = bbox.getY0();

        if (bbox.getDimensions() != _wimage->getDimensions()) {
            throw LSST_EXCEPT(pexExceptions::LengthErrorException,
                              (boost::format("Footprint at %d,%d -- %d,%d is wrong size for "
                                             "%d x %d weight image") %
                               bbox.getX0() % bbox.getY0() % bbox.getX1() % bbox.getY1() %
                               _wimage->getWidth() % _wimage->getHeight()).str());
        }
    }
    void reset() {}
    
    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        typename WeightImageT::Pixel wval = (*_wimage)(x - _x0, y - _y0);
        _sum += wval*ival;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    typename WeightImageT::Ptr const& _wimage;        // The weight image
    double _sum;                                      // our desired sum
    int _x0, _y0;                                     // the origin of the current Footprint
};

            
/*****************************************************************************************************/
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

/************************************************************************************************************/
/**
 * Process the image; calculate values
 */
template<typename MaskedImageT>
detection::Photometry::Ptr
NaivePhotometry::doMeasure(typename MaskedImageT::ConstPtr img, detection::Peak const& peak)
{
    double const xcen = peak.getFx();   ///< object's column position
    double const ycen = peak.getFy();   ///< object's row position
    PSF const *psf = NULL;              ///< image's PSF

    int const ixcen = afwImage::positionToIndex(xcen);
    int const iycen = afwImage::positionToIndex(ycen);
    
    afwImage::BBox imageBBox(afwImage::PointI(img->getX0(), img->getY0()),
                             img->getWidth(), img->getHeight()); // BBox for data image

    /* ******************************************************* */
    // Aperture photometry
    double aperFlux = std::numeric_limits<double>::quiet_NaN();
    {
        FootprintFlux<MaskedImageT> fluxFunctor(*img);
        
        double const radius = 10;       // == getRadius()
        detection::Footprint const foot(afwImage::BCircle(afwImage::PointI(ixcen, iycen), radius),
                                        imageBBox);
        fluxFunctor.apply(foot);
        aperFlux = fluxFunctor.getSum();
    }

    /* ******************************************************** */
    // Weighted aperture photometry, using a PSF weight --- i.e. a PSF flux
    double psfFlux = std::numeric_limits<double>::quiet_NaN();
    if (psf) {
        PSF::Image::Ptr wimage = psf->getImage(xcen, ycen);
        
        FootprintWeightFlux<MaskedImageT, PSF::Image> wfluxFunctor(*img, wimage);
        
        // Build a rectangular Footprint corresponding to wimage
        detection::Footprint foot(afwImage::BBox(afwImage::PointI(0, 0),
                                              psf->getWidth(), psf->getHeight()), imageBBox);
        foot.shift(ixcen - psf->getWidth()/2, iycen - psf->getHeight()/2);
        
        wfluxFunctor.apply(foot);
        
        getSum2<PSF::Pixel> sum;
        sum = std::accumulate(wimage->begin(true), wimage->end(true), sum);
        psfFlux = wfluxFunctor.getSum()/sum.sum2;
    }

    return boost::make_shared<NaivePhotometry>(aperFlux, psfFlux);
}

//
// Explicit instantiations
//
// We need to make an instance here so as to register it with MeasurePhotometry
//
// \cond
#define INSTANTIATE(TYPE) \
    NewMeasurePhotometry<afwImage::MaskedImage<TYPE> >::declare("NAIVE", \
                                                    &NaivePhotometry::doMeasure<afwImage::MaskedImage<TYPE> >)

volatile bool isInstance[] = {
    INSTANTIATE(float),
    INSTANTIATE(double)
};

// \endcond

}}}}
