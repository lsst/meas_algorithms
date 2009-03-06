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
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/Shape.h"
#include "lsst/meas/algorithms/Photometry.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * A class to provide a set of flags describing our processing
 *
 * This would be in MeasureSources, but that's inconvenient as it's templated
 */
struct Flags {
    enum {
        EDGE                      = 0x1,    ///< object is in region labelled EDGE
        SHAPE_SHIFT               = 0x2,    ///< centroid shifted while estimating adaptive moments
        SHAPE_MAXITER             = 0x4,    ///< too many iterations for adaptive moments
        SHAPE_UNWEIGHTED          = 0x8,    ///< "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_PSF      = 0x10,   ///< the PSF's "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_BAD      = 0x20,   ///< even the unweighted moments were bad
        PEAKCENTER                = 0x40,   ///< given centre is position of peak pixel
        BINNED1                   = 0x80,   ///< object was found in 1x1 binned image
        INTERP                    = 0x100,  ///< object's footprint includes interpolated pixels
        INTERP_CENTER             = 0x200,  ///< object's centre is close to interpolated pixels
        SATUR                     = 0x400,  ///< object's footprint includes saturated pixels
        SATUR_CENTER              = 0x800,  ///< object's centre is close to saturated pixels
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
        _exposure(exposure), _policy(policy), _psf(psf),
        _moLog(lsst::pex::logging::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                                       lsst::pex::logging::Log::INFO)) {
        //
        // lookup algorithms in Policy
        //
        _mCentroid = 
            createMeasureCentroid<typename MaskedImageT::Image>(_policy.getString("measureObjects.centroidAlgorithm"));
        _mShape = createMeasureShape<MaskedImageT>(_policy.getString("measureObjects.shapeAlgorithm"));
        _mPhotometry = createMeasurePhotometry<MaskedImageT>(_policy.getString("measureObjects.photometryAlgorithm"),
                                                             _policy.getDouble("measureObjects.apRadius"));
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
    measureCentroid<typename MaskedImageT::Image>* getMeasureCentroid() const { return _mCentroid; }
    /// return the shape measurer
    measureShape<MaskedImageT>* getMeasureShape() const { return _mShape; }
    measurePhotometry<MaskedImageT>* getMeasurePhotometry() const { return _mPhotometry; }

private:
    ExposureT const _exposure;              // Exposure wherein Sources dwell
    lsst::pex::policy::Policy const _policy;// Policy to describe processing
    PSF::ConstPtr _psf;                     // image's PSF \todo Cf #645

    boost::shared_ptr<lsst::pex::logging::Log> _moLog; // log for measureObjects
    /*
     * Objects that know how to measure the object's properties
     */
    measureCentroid<typename MaskedImageT::Image> * _mCentroid;
    measureShape<MaskedImageT> * _mShape;
    measurePhotometry<MaskedImageT> * _mPhotometry;
};

/**
 * A function to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
typename MeasureSources<ExposureT>::Ptr makeMeasureSources(
        ExposureT const& exposure,
        lsst::pex::policy::Policy const& policy,
        PSF::ConstPtr psf=PSF::ConstPtr()
    ) {
    return typename MeasureSources<ExposureT>::Ptr(new MeasureSources<ExposureT>(exposure, policy, psf));
}

}}}
#endif // LSST_DETECTION_MEASURE_H
