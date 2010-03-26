// -*- LSST-C++ -*-
#if !defined(LSST_DETECTION_MEASURE_H)
#define LSST_DETECTION_MEASURE_H

//!
// Measure properties of an image selected by a Footprint
//
#include <list>
#include <cmath>
#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/Shape.h"
#include "lsst/meas/algorithms/Photometry.h"

namespace lsst {
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
    typedef boost::shared_ptr<NewMeasureAstrometry> Ptr;

    NewMeasureAstrometry(typename ImageT::ConstPtr im) :
        lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Astrometry,
                                              ImageT, lsst::afw::detection::Peak>(im) {}
};
    
/**
 * Here's the object that remembers and can execute our choice of photometric algorithms
 */
template<typename ImageT>
class NewMeasurePhotometry :
        public lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Photometry,
                                                     ImageT, lsst::afw::detection::Peak> {
public:
    typedef boost::shared_ptr<NewMeasurePhotometry> Ptr;
    
    NewMeasurePhotometry(typename ImageT::ConstPtr im) :
        lsst::afw::detection::MeasureQuantity<lsst::afw::detection::Photometry,
                                              ImageT, lsst::afw::detection::Peak>(im) {}
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
        /// Should this this object be ignored in essentially all analyses?
        BAD                       = EDGE|INTERP_CENTER|SATUR_CENTER
    };
};

template<typename ExposureT>
class MeasureSources {
public:
    typedef boost::shared_ptr<MeasureSources> Ptr;
    typedef boost::shared_ptr<MeasureSources const> ConstPtr;

    typedef typename ExposureT::MaskedImageT MaskedImageT;

    MeasureSources(ExposureT const& exposure,                   ///< Exposure wherein Sources dwell
                   lsst::pex::policy::Policy const& policy,     ///< Policy to describe processing
                   PSF::ConstPtr psf                            ///< image's PSF \todo Cf #645
                  ) :
        _exposure(exposure), _policy( policy), _psf(psf),
        _moLog(lsst::pex::logging::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                                       lsst::pex::logging::Log::INFO)) {
        //
        // lookup algorithms in Policy
        //
        _mCentroid = 
            createMeasureCentroid<typename MaskedImageT::Image>(_policy.getString("centroidAlgorithm"), exposure.getMaskedImage().getImage());
        {
            typename MaskedImageT::Ptr mi(new MaskedImageT(exposure.getMaskedImage()));

            _mShape = 
                createMeasureShape<MaskedImageT>(_policy.getString("shapeAlgorithm"), mi);
            _mPhotometry = createMeasurePhotometry<MaskedImageT>(_policy.getString("photometryAlgorithm"), mi);
            _mPhotometry->setRadius(_policy.getDouble("apRadius"));
        }
        /****************************************************************************************************/
        {
            typename MaskedImageT::Ptr mi(new MaskedImageT(exposure.getMaskedImage()));

            typename NewMeasureAstrometry<MaskedImageT>::Ptr measureAstro =
                boost::make_shared<NewMeasureAstrometry<MaskedImageT> >(mi);
#if 0
            measureAstro->addAlgorithm("naive");
#endif
            
            // Create our photometric measuring object based on argv
            NewMeasurePhotometry<MaskedImageT> *measurePhoto = new NewMeasurePhotometry<MaskedImageT>(mi);
#if 0    
            for (int i = 1; i != argc; ++i) {
                measurePhoto->addAlgorithm(argv[i]);
            }
#endif
        }
    }
    
    virtual ~MeasureSources() {}
    
    virtual void apply(lsst::afw::detection::Source::Ptr src,   ///< the Source to receive results
                       lsst::afw::detection::Footprint const& foot  ///< Footprint to measure
                      );
    
    /// Return the Exposure
    ExposureT const& getExposure() const { return _exposure; }
    /// Return the Policy used to describe processing
    lsst::pex::policy::Policy const& getPolicy() const { return _policy; }
    /// Return the PSF
    PSF::ConstPtr getPsf() const { return _psf; }
    /// Return the log
    lsst::pex::logging::Log &getLog() const { return *_moLog; }
    /// return the centroid measurer
    MeasureCentroid<typename MaskedImageT::Image>* getMeasureCentroid() const { return _mCentroid; }
    /// return the shape measurer
    MeasureShape<MaskedImageT>* getMeasureShape() const { return _mShape; }
    MeasurePhotometry<MaskedImageT>* getMeasurePhotometry() const { return _mPhotometry; }

private:
    ExposureT const _exposure;              // Exposure wherein Sources dwell
    lsst::pex::policy::Policy const _policy;// Policy to describe processing
    PSF::ConstPtr _psf;                     // image's PSF \todo Cf #645

    boost::shared_ptr<lsst::pex::logging::Log> _moLog; // log for measureObjects
    /*
     * Objects that know how to measure the object's properties
     */
    MeasureCentroid<typename MaskedImageT::Image> * _mCentroid;
    MeasureShape<MaskedImageT> * _mShape;
    MeasurePhotometry<MaskedImageT> * _mPhotometry;
};

/**
 * A function to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
typename MeasureSources<ExposureT>::Ptr makeMeasureSources(
        ExposureT const& exposure,
        lsst::pex::policy::Policy const& policy,
        PSF::ConstPtr psf = PSF::ConstPtr()
    ) {
    return typename MeasureSources<ExposureT>::Ptr(new MeasureSources<ExposureT>(exposure, policy, psf));
}
            
}}}
#endif // LSST_DETECTION_MEASURE_H
