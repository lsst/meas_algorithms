/// \file

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/Shape.h"
#include "lsst/meas/algorithms/Photometry.h"

namespace lsst { namespace meas { namespace algorithms {
namespace pexLogging = lsst::pex::logging;
namespace image = lsst::afw::image;
namespace detection = lsst::afw::detection;
/************************************************************************************************************/
/**
 * \brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintCentroid : public detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintCentroid(MaskedImageT const& mimage ///< The image the source lives in
                              ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                                  _n(0), _sum(0), _sumx(0), _sumy(0),
                                  _min( std::numeric_limits<double>::max()), _xmin(0), _ymin(0),
                                  _max(-std::numeric_limits<double>::max()), _xmax(0), _ymax(0),
                                  _bits(0)
        {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _n = 0;
        _sum = _sumx = _sumy = 0.0;
        _min =  std::numeric_limits<double>::max();
        _xmin = _ymin = 0;
        _max = -std::numeric_limits<double>::max();
        _xmax = _ymax = 0;
        _bits = 0x0;
    }

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

        _n++;
        _sum += val;
        _sumx += lsst::afw::image::indexToPosition(x)*val;
        _sumy += lsst::afw::image::indexToPosition(y)*val;
        _bits |= loc.mask(0, 0);

        if (val < _min) {
            _min = val;
            _xmin = x;
            _ymin = y;
        }
        if (val > _max) {
            _max = val;
            _xmax = x;
            _ymax = y;
        }
    }

    /// Return the number of pixels
    int getN() const { return _n; }
    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the Footprint's column centroid
    double getX() const { return _sumx/_sum; }
    /// Return the Footprint's row centroid
    double getY() const { return _sumy/_sum; }
    /// Return the Footprint's peak pixel
    detection::Peak getPeak(bool isNegative) const {
        return isNegative ? detection::Peak(_xmin, _ymin) : detection::Peak(_xmax, _ymax);
    }
    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    int _n;
    double _sum, _sumx, _sumy;
    double _min;
    int _xmin, _ymin;
    double _max;
    int _xmax, _ymax;
    typename MaskedImageT::Mask::Pixel _bits;
};


template <typename MaskedImageT>
class FootprintFlux : public detection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : detection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0)
        {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _sum = 0.0;
    }

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);
        _sum += val;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    double _sum;
};

/************************************************************************************************************/
/**
 * Use *this to measure the Footprint foot, setting fields in src
 */
template<typename MaskedImageT>
void MeasureSources<MaskedImageT>::apply(
        lsst::afw::detection::Source::Ptr src,       ///< the Source to receive results
        lsst::afw::detection::Footprint const& foot  ///< Footprint to measure
                                                                   ) {
    float background = 0;               // background level to subtract XXX
    MaskedImageT const& mimage = getExposure().getMaskedImage();
    PSF::ConstPtr psf = getPsf();

    bool const isNegative = (src->getFlagForDetection() & Flags::DETECT_NEGATIVE);
    //
    // Measure some properties of the Footprint
    //
    FootprintCentroid<MaskedImageT> centroidFunctor(mimage);
    centroidFunctor.apply(foot);

    detection::Peak const& peak = centroidFunctor.getPeak(isNegative);
    //
    // Check for bits set in the Footprint
    //
    if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::EDGE);
    }
    if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::INTERP);
    }
    if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::SATUR);
    }
    //
    // Now run measure objects code (but not for edge objects)
    //
    typename MaskedImageT::Mask &mask = *mimage.getMask();
    if (mask(peak.getIx() - mask.getX0(), peak.getIy() - mask.getY0(),
             MaskedImageT::Mask::getMaskPlane("EDGE"))) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::EDGE);
        return;
    }
    //
    // Centroids
    //
    try {
        Centroid cen = getMeasureCentroid()->apply(*mimage.getImage(),
                                                   peak.getIx(), peak.getIy(), psf.get(), background);
        
        src->setXAstrom(cen.getX());
        src->setXAstromErr(cen.getXErr());
        src->setYAstrom(cen.getY());
        src->setYAstromErr(cen.getYErr());
    } catch (lsst::pex::exceptions::LengthErrorException const& e) {
        src->setXAstrom(peak.getIx());
        src->setYAstrom(peak.getIy());
        src->setFlagForDetection(src->getFlagForDetection() | (Flags::EDGE | Flags::PEAKCENTER));

        return;
    } catch (lsst::pex::exceptions::RuntimeErrorException const& e) {
        src->setXAstrom(peak.getIx());
        src->setYAstrom(peak.getIy());
        src->setFlagForDetection(src->getFlagForDetection() | Flags::PEAKCENTER);

        return;
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Centroiding at (%d, %d)") % peak.getIx() % peak.getIy()).str());
        throw e;
    }
    //
    // Shapes
    //
    try {
        Shape shape = getMeasureShape()->apply(mimage,
                                               src->getXAstrom(), src->getYAstrom(), psf.get(), background);
        
        src->setIxx(shape.getMxx());       // <xx>
        src->setIxxErr(shape.getMxxErr()); // sqrt(Var<xx>)
        src->setIxy(shape.getMxy());       // <xy>
        src->setIxyErr(shape.getMxyErr()); // sign(Covar(x, y))*sqrt(|Covar(x, y)|))        
        src->setIyy(shape.getMyy());       // <yy>
        src->setIyyErr(shape.getMyyErr()); // sqrt(Var<yy>)

        src->setFlagForDetection(src->getFlagForDetection() | shape.getFlags());
    } catch (lsst::pex::exceptions::DomainErrorException const& e) {
        getLog().log(pexLogging::Log::INFO, boost::format("Measuring Shape at (%.3f,%.3f): %s") %
                     src->getXAstrom() % src->getYAstrom() % e.what());
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Measuring Shape at (%.3f, %.3f)") %
                            src->getXAstrom() % src->getYAstrom()).str());
        throw e;
    }


    //
    // Photometry
    //
    try {
        Photometry photometry = getMeasurePhotometry()->apply(mimage, src->getXAstrom(), src->getYAstrom(),
                                                              psf.get(), background);
        
        src->setApFlux(photometry.getApFlux());
        //src->setApFluxErr(photometry.getApMagErr());
        src->setPsfFlux(photometry.getPsfFlux());
        //src->setPsfFluxErr(photometry.getApPsfErr());
        
    } catch (lsst::pex::exceptions::DomainErrorException const& e) {
        getLog().log(pexLogging::Log::INFO, boost::format("Measuring Photometry at (%.3f,%.3f): %s") %
                     src->getXAstrom() % src->getYAstrom() % e.what());
    } catch (lsst::pex::exceptions::Exception & e) {
        LSST_EXCEPT_ADD(e, (boost::format("Measuring Photometry at (%.3f, %.3f)") %
                            src->getXAstrom() % src->getYAstrom()).str());
        throw e;
    }
    
    //
    // Check for bits set near the centroid
    //
    {
        image::PointI llc(image::positionToIndex(src->getXAstrom()) - 1,
                          image::positionToIndex(src->getYAstrom()) - 1);
        detection::Footprint const middle(image::BBox(llc, 3, 3)); // 3x3 square centred at the the centroid
        centroidFunctor.apply(middle);
        if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
            src->setFlagForDetection(src->getFlagForDetection() | Flags::INTERP_CENTER);
        }
        if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
            src->setFlagForDetection(src->getFlagForDetection() | Flags::SATUR_CENTER);
        }
    }
}

//
// Explicit instantiations
//
// \cond
template class MeasureSources<lsst::afw::image::Exposure<float> >;
// \endcond
}}}
