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
 
#if !defined(LSST_MEAS_ALGORITHMS_MEASURE_H)
#define LSST_MEAS_ALGORITHMS_MEASURE_H

//!
// Measure properties of an image selected by a Footprint
//
#include <list>
#include <cmath>
#include "lsst/base.h"
#include "boost/cstdint.hpp"
#include "boost/type_traits.hpp"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/meas/algorithms/MeasureQuantity.h"

namespace {
    template<typename T>
    struct ElementTypeNoCV {            // Remove top-level cv-qualifiers
        typedef typename boost::remove_cv<T>::type type;
    };

    template<typename T>
    struct ElementTypeNoCV<boost::shared_ptr<T> > { // Return the type of a shared_ptr with cv-qualifiers removed
        typedef typename boost::remove_cv<T>::type type;
    };
}

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
/*
 * Defining e.g. makeMeasureAstrometry's a little tricky, as the whole point of helper functions a l\'a make_pair
 * is to allow the compiler to deduce the template arguments, and the constructor for e.g. MeasureAstrometry accepts
 * a shared_ptr<Foo>, not a Foo & and the deduction fails
 *
 * In C++ we work around this by templating on the shared pointer and using ElementTypeNoCV to deduce the desired
 * underlying type.  Unfortunately, swig doesn't understand this and generates unusable bindings.  It's probable
 * that this could be solved in pure swig, but it's simpler to provide a different version of makeMeasureAstrometry
 * that it can handle directly.  Swig handles argument checking differently from C++, and the simpler version works
 * just fine from python
 */
#if defined(DOXYGEN) || defined(SWIGPYTHON)
#define MAKE_MEASURE_ALGORITHM(ALGORITHM) \
    template<typename ImageT> typename boost::shared_ptr<Measure##ALGORITHM<ImageT> > \
    makeMeasure##ALGORITHM(typename ImageT::ConstPtr im=typename ImageT::ConstPtr(), \
                          CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) \
    { \
        return boost::make_shared<Measure##ALGORITHM<ImageT> >(im, policy); \
    }
#else
#define MAKE_MEASURE_ALGORITHM(ALGORITHM) \
    template<typename ImageConstPtr> \
    typename boost::shared_ptr<Measure##ALGORITHM<typename ElementTypeNoCV<ImageConstPtr>::type> > \
    makeMeasure##ALGORITHM(ImageConstPtr im=ImageConstPtr(),                       /**< "Image::ConstPtr" */ \
        CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)() /**< policy to configure with */ \
                         )                                              \
    {                                                                   \
        typedef typename ElementTypeNoCV<ImageConstPtr>::type Image;    \
                                                                        \
        typename Image::ConstPtr cim = im; /* Ensure that we have a ConstPtr even if passed a non-const Ptr */ \
        return boost::make_shared<Measure##ALGORITHM<Image> >(cim, policy); \
    }
#endif

/************************************************************************************************************/
/**
 * Here's the object that remembers and can execute our choice of astrometric algorithms
 */
template<typename ImageT>
class MeasureAstrometry :
        public MeasureQuantity<lsst::afw::detection::Astrometry, ImageT> {
public:
    typedef PTR(MeasureAstrometry) Ptr;

    MeasureAstrometry(typename ImageT::ConstPtr im,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()
                     ) :
        MeasureQuantity<lsst::afw::detection::Astrometry, ImageT>(im, policy) {}
};

MAKE_MEASURE_ALGORITHM(Astrometry)
    
/**
 * Here's the object that remembers and can execute our choice of photometric algorithms
 */
template<typename ImageT>
class MeasurePhotometry :
        public MeasureQuantity<lsst::afw::detection::Photometry, ImageT> {
public:
    typedef PTR(MeasurePhotometry) Ptr;
    
    MeasurePhotometry(typename ImageT::ConstPtr im,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()
                     ) :
        MeasureQuantity<lsst::afw::detection::Photometry, ImageT>(im, policy) {}
};

MAKE_MEASURE_ALGORITHM(Photometry)

/**
 * Here's the object that remembers and can execute our choice of shape measurement algorithms
 */
template<typename ImageT>
class MeasureShape :
        public MeasureQuantity<lsst::afw::detection::Shape, ImageT> {
public:
    typedef PTR(MeasureShape) Ptr;

    MeasureShape(typename ImageT::ConstPtr im,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()
                     ) :
        MeasureQuantity<lsst::afw::detection::Shape, ImageT>(im, policy) {}
};

MAKE_MEASURE_ALGORITHM(Shape)
    
/************************************************************************************************************/
/**
 * A class to provide a set of flags describing our processing
 *
 * This would be in MeasureSources, but that's inconvenient as it's templated
 */
struct Flags {
    enum {
        EDGE                            = 0x000001, ///< source is in region labelled EDGE
        SHAPE_SHIFT                     = 0x000002, ///< centroid shifted while estimating adaptive moments
        SHAPE_MAXITER                   = 0x000004, ///< too many iterations for adaptive moments
        SHAPE_UNWEIGHTED                = 0x000008, ///< "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_PSF            = 0x000010, ///< the PSF's "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_BAD            = 0x000020, ///< even the unweighted moments were bad
        PEAKCENTER                      = 0x000040, ///< given centre is position of peak pixel
        BINNED1                         = 0x000080, ///< source was found in 1x1 binned image
        INTERP                          = 0x000100, ///< source's footprint includes interpolated pixels
        INTERP_CENTER                   = 0x000200, ///< source's centre is close to interpolated pixels
        SATUR                           = 0x000400, ///< source's footprint includes saturated pixels
        SATUR_CENTER                    = 0x000800, ///< source's centre is close to saturated pixels
        DETECT_NEGATIVE                 = 0x001000, ///< source was detected as being significantly negative
        STAR                            = 0x002000, ///< source is thought to be point-like
        PSFSTAR                         = 0x004000, ///< source was used in PSF determination

        PHOTOM_NO_PSF                   = 0x008000, ///< NO Psf provided to photometry algorithm
        PHOTOM_NO_PEAK                  = 0x010000, ///< NO Peak provided to photometry algorithm
        PHOTOM_NO_SOURCE                = 0x020000, ///< NO source provided to photometry algorithm
        PHOTOM_NO_FOOTPRINT             = 0x040000, ///< NO FOOTPRINT provided to photometry algorithm

        SHAPELET_PHOTOM_NO_BASIS        = 0x080000, ///< ShapeletModelPhotometry configure without a basis
        SHAPELET_PHOTOM_BAD_MOMENTS     = 0x100000, ///< input moments are too large or not finite
        SHAPELET_PHOTOM_INVERSION_FAIL  = 0x200000, ///< ShapeletModelPhotometry failed
        SHAPELET_PHOTOM_INVERSION_UNSAFE= 0x400000, ///< ShapeletModelPhotometry should not be trusted
        SHAPELET_PHOTOM_GALAXY_FAIL     = 0x800000, ///< ShapeletModelPhotometry only fit a point source model
        ///ShapeletModelPhtoometry should be ignored in essentially all analyses
        SHAPELET_PHOTOM_BAD = PHOTOM_NO_PSF | PHOTOM_NO_SOURCE | PHOTOM_NO_FOOTPRINT | 
                SHAPELET_PHOTOM_NO_BASIS | SHAPELET_PHOTOM_BAD_MOMENTS | 
                SHAPELET_PHOTOM_INVERSION_FAIL | SHAPELET_PHOTOM_INVERSION_UNSAFE,
        /// Should this this object be ignored in essentially all analyses?
        BAD                       = EDGE|INTERP_CENTER|SATUR_CENTER

    };
};

template<typename ExposureT>
class MeasureSources {
public:
    typedef PTR(MeasureSources) Ptr;
    typedef CONST_PTR(MeasureSources) ConstPtr;

    MeasureSources(typename ExposureT::ConstPtr exposure,   ///< Exposure wherein Sources dwell
                   lsst::pex::policy::Policy const& policy ///< Policy to describe processing
                  ) :
        _exposure(exposure), _policy( policy),
        _moLog(lsst::pex::logging::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                                       lsst::pex::logging::Log::INFO)) {

        lsst::pex::policy::DefaultPolicyFile dictFile(
            "meas_algorithms", "MeasureSourcesDictionary.paf", "policy");
        CONST_PTR(lsst::pex::policy::Policy) dictPtr(
            lsst::pex::policy::Policy::createPolicy(
                dictFile, dictFile.getRepositoryPath()));

        lsst::pex::policy::DefaultPolicyFile defaultsFile(
            "meas_algorithms", "MeasureSourcesDefaults.paf", "policy");
        CONST_PTR(lsst::pex::policy::Policy) defaultsPtr(
            lsst::pex::policy::Policy::createPolicy(
                defaultsFile, defaultsFile.getRepositoryPath()));

        _policy.mergeDefaults(*defaultsPtr);
        _policy.mergeDefaults(*dictPtr);
        
        if (_policy.isPolicy("astrometry")) {
            _measureAstrom =
                boost::make_shared<MeasureAstrometry<ExposureT> >(exposure, _policy.getPolicy("astrometry"));
        }
        
        if (_policy.isPolicy("photometry")) {
            _measurePhotom =
                boost::make_shared<MeasurePhotometry<ExposureT> >(exposure, _policy.getPolicy("photometry"));
        }

        if (_policy.isPolicy("shape")) {
            _measureShape =
                boost::make_shared<MeasureShape<ExposureT> >(exposure, _policy.getPolicy("shape"));
        }
    }
    
    virtual ~MeasureSources() {
    }
    
    virtual void apply(PTR(lsst::afw::detection::Source) src,   ///< the Source to receive results
                       CONST_PTR(lsst::afw::detection::Footprint) foot=
                             PTR(lsst::afw::detection::Footprint)() ///< Footprint to measure
                      );
    
    /// Return the Exposure
    typename ExposureT::ConstPtr getExposure() const { return _exposure; }
    /// Set the Exposure
    void setExposure(typename ExposureT::ConstPtr exposure) {
        _exposure = exposure;

        if (_measureAstrom) {
            _measureAstrom->setImage(exposure);
        }
        if (_measurePhotom) {
            _measurePhotom->setImage(exposure);
        }
        if (_measureShape) {
            _measureShape->setImage(exposure);
        }
    }
    /// Return the Policy used to describe processing
    lsst::pex::policy::Policy const& getPolicy() const { return _policy; }
    /// Return the log
    lsst::pex::logging::Log &getLog() const { return *_moLog; }
    /// return the astrometric measurer
    typename MeasureAstrometry<ExposureT>::Ptr getMeasureAstrom() const { return _measureAstrom; }
    /// return the photometric measurer
    typename MeasurePhotometry<ExposureT>::Ptr getMeasurePhotom() const { return _measurePhotom; }
    /// return the shape measurer
    typename MeasureShape<ExposureT>::Ptr getMeasureShape() const { return _measureShape; }

private:
    typename ExposureT::ConstPtr _exposure;    // Exposure wherein Sources dwell
    lsst::pex::policy::Policy _policy;   // Policy to describe processing

    PTR(lsst::pex::logging::Log) _moLog; // log for measureObjects
    /*
     * Objects that know how to measure the object's properties
     */
    typename MeasureAstrometry<ExposureT>::Ptr _measureAstrom;
    typename MeasurePhotometry<ExposureT>::Ptr _measurePhotom;
    typename MeasureShape<ExposureT>::Ptr      _measureShape;
};

/**
 * A function to return a MeasureSources of the correct type (cf. std::make_pair)
 */
#if defined(DOXYGEN) || defined(SWIGPYTHON)
    template<typename Exposure>
    typename MeasureSources<Exposure>::Ptr makeMeasureSources(
        typename Exposure::ConstPtr exposure,
        lsst::pex::policy::Policy const& policy
                                                             )
    {
        return typename MeasureSources<Exposure>::Ptr(new MeasureSources<Exposure>(exposure, policy));
    }
#else    
template<typename ExposureConstPtr>
typename MeasureSources<typename ElementTypeNoCV<ExposureConstPtr>::type>::Ptr makeMeasureSources(
        ExposureConstPtr exposure,
        lsst::pex::policy::Policy const& policy
    ) {
    typedef typename ElementTypeNoCV<ExposureConstPtr>::type Exposure;

    typename Exposure::ConstPtr cexposure = exposure; // Ensure that we have a ConstPtr even if passed a non-const Ptr
    return typename MeasureSources<Exposure>::Ptr(new MeasureSources<Exposure>(cexposure, policy));
}
#endif
       
}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
