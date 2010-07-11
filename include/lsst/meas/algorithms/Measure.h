// -*- LSST-C++ -*-
#if !defined(LSST_DETECTION_MEASURE_H)
#define LSST_DETECTION_MEASURE_H

//!
// Measure properties of an image selected by a Footprint
//
#include <list>
#include <cmath>
#include "lsst/base.h"
#include "boost/cstdint.hpp"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"

namespace lsst {
namespace pex {
    namespace policy {
        class Policy;
    }
}
namespace afw {
    namespace detection {
        class Psf;
    }
}
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * Here's the object that remembers and can execute our choice of astrometric algorithms
 */
template<typename ImageT>
class NewMeasureAstrometry :
        public lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Astrometry,
                                                     ImageT, lsst::afw::detection::Peak> {
public:
    typedef PTR(NewMeasureAstrometry) Ptr;

    NewMeasureAstrometry(typename ImageT::ConstPtr im,
                         CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()
                        ) :
        lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Astrometry,
                                              ImageT, lsst::afw::detection::Peak>(im, policy) {}
};
    
/**
 * Here's the object that remembers and can execute our choice of photometric algorithms
 */
template<typename ImageT>
class NewMeasurePhotometry :
        public lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Photometry,
                                                     ImageT, lsst::afw::detection::Peak> {
public:
    typedef PTR(NewMeasurePhotometry) Ptr;
    
    NewMeasurePhotometry(typename ImageT::ConstPtr im,
                         CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()
                        ) :
        lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Photometry,
                                              ImageT, lsst::afw::detection::Peak>(im, policy) {}
};

/************************************************************************************************************/
/**
 * A class to provide a set of flags describing our processing
 *
 * This would be in MeasureSources, but that's inconvenient as it's templated
 */
struct Flags {
    enum {
        EDGE                      = 0x1,    ///< source is in region labelled EDGE
        SHAPE_SHIFT               = 0x2,    ///< centroid shifted while estimating adaptive moments
        SHAPE_MAXITER             = 0x4,    ///< too many iterations for adaptive moments
        SHAPE_UNWEIGHTED          = 0x8,    ///< "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_PSF      = 0x10,   ///< the PSF's "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_BAD      = 0x20,   ///< even the unweighted moments were bad
        PEAKCENTER                = 0x40,   ///< given centre is position of peak pixel
        BINNED1                   = 0x80,   ///< source was found in 1x1 binned image
        INTERP                    = 0x100,  ///< source's footprint includes interpolated pixels
        INTERP_CENTER             = 0x200,  ///< source's centre is close to interpolated pixels
        SATUR                     = 0x400,  ///< source's footprint includes saturated pixels
        SATUR_CENTER              = 0x800,  ///< source's centre is close to saturated pixels
        DETECT_NEGATIVE           = 0x1000, ///< source was detected as being significantly negative
        STAR                      = 0x2000, ///< source is thought to be point-like
        /// Should this this object be ignored in essentially all analyses?
        BAD                       = EDGE|INTERP_CENTER|SATUR_CENTER
    };
};

template<typename ExposureT>
class MeasureSources {
public:
    typedef PTR(MeasureSources) Ptr;
    typedef CONST_PTR(MeasureSources) ConstPtr;

    typedef typename ExposureT::MaskedImageT MaskedImageT;

    MeasureSources(ExposureT const& exposure,                   ///< Exposure wherein Sources dwell
                   lsst::pex::policy::Policy const& policy,     ///< Policy to describe processing
                   CONST_PTR(lsst::afw::detection::Psf) psf     ///< image's PSF \todo Cf #645
                  ) :
        _exposure(exposure), _policy( policy), _psf(psf),
        _moLog(lsst::pex::logging::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                                       lsst::pex::logging::Log::INFO)) {

        typename MaskedImageT::Ptr mi(new MaskedImageT(exposure.getMaskedImage()));
        
        _measureAstrom =
            boost::make_shared<NewMeasureAstrometry<MaskedImageT> >(mi, _policy.getPolicy("astrometry"));
        
        _measurePhotom =
            boost::make_shared<NewMeasurePhotometry<MaskedImageT> >(mi, _policy.getPolicy("photometry"));
#if 0
        _measureShape =
            boost::make_shared<NewMeasureShape<MaskedImageT> >(mi, _policy.getPolicy("shape"));
#endif
    }
    
    virtual ~MeasureSources() {
    }
    
    virtual void apply(lsst::afw::detection::Source::Ptr src,   ///< the Source to receive results
                       lsst::afw::detection::Footprint const& foot  ///< Footprint to measure
                      );
    
    /// Return the Exposure
    ExposureT const& getExposure() const { return _exposure; }
    /// Return the Policy used to describe processing
    lsst::pex::policy::Policy const& getPolicy() const { return _policy; }
    /// Return the PSF
    CONST_PTR(lsst::afw::detection::Psf) getPsf() const { return _psf; }
    /// Return the log
    lsst::pex::logging::Log &getLog() const { return *_moLog; }
    /// return the astrometric measurer
    typename NewMeasureAstrometry<MaskedImageT>::Ptr getMeasureAstrom() const { return _measureAstrom; }
    /// return the photometric measurer
    typename NewMeasurePhotometry<MaskedImageT>::Ptr getMeasurePhotom() const { return _measurePhotom; }

private:
    ExposureT const _exposure;              // Exposure wherein Sources dwell
    lsst::pex::policy::Policy const _policy;// Policy to describe processing
    CONST_PTR(lsst::afw::detection::Psf) _psf; // image's PSF \todo Cf #645

    PTR(lsst::pex::logging::Log) _moLog; // log for measureObjects
    /*
     * Objects that know how to measure the object's properties
     */
    typename NewMeasureAstrometry<MaskedImageT>::Ptr _measureAstrom;
    typename NewMeasurePhotometry<MaskedImageT>::Ptr _measurePhotom;
};

/**
 * A function to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
typename MeasureSources<ExposureT>::Ptr makeMeasureSources(
        ExposureT const& exposure,
        lsst::pex::policy::Policy const& policy,
        CONST_PTR(lsst::afw::detection::Psf) psf = CONST_PTR(lsst::afw::detection::Psf)()
    ) {
    return typename MeasureSources<ExposureT>::Ptr(new MeasureSources<ExposureT>(exposure, policy, psf));
}
            
}}}
#endif // LSST_DETECTION_MEASURE_H
