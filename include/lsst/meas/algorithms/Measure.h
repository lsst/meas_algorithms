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
#include "lsst/meas/algorithms/Flags.h"
#include "lsst/meas/algorithms/detail/Measure.h"

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

namespace pexPolicy = lsst::pex::policy;
namespace pexLog = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;

template<typename ExposureT>
class MeasureSources {
public:
    typedef PTR(MeasureSources) Ptr;
    typedef CONST_PTR(MeasureSources) ConstPtr;
    typedef MeasureAstrometry<ExposureT> MeasureAstrometryT;
    typedef MeasureShape<ExposureT> MeasureShapeT;
    typedef MeasurePhotometry<ExposureT> MeasurePhotometryT;

    MeasureSources(pexPolicy::Policy const& policy ///< Policy to describe processing
                  ) :
        _policy( policy),
        _moLog(pexLog::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                                       pexLog::Log::INFO)) {

        pexPolicy::DefaultPolicyFile dictFile(
            "meas_algorithms", "MeasureSourcesDictionary.paf", "policy");
        CONST_PTR(pexPolicy::Policy) dictPtr(
            pexPolicy::Policy::createPolicy(
                dictFile, dictFile.getRepositoryPath()));

        pexPolicy::DefaultPolicyFile defaultsFile(
            "meas_algorithms", "MeasureSourcesDefaults.paf", "policy");
        CONST_PTR(pexPolicy::Policy) defaultsPtr(
            pexPolicy::Policy::createPolicy(
                defaultsFile, defaultsFile.getRepositoryPath()));

        _policy.mergeDefaults(*defaultsPtr);
        _policy.mergeDefaults(*dictPtr);
        
        if (_policy.isPolicy("astrometry")) {
            _measureAstrom = boost::make_shared<MeasureAstrometryT>(_policy.getPolicy("astrometry"));
        }
        
        if (_policy.isPolicy("photometry")) {
            _measurePhotom = boost::make_shared<MeasurePhotometryT>(_policy.getPolicy("photometry"));
        }

        if (_policy.isPolicy("shape")) {
            _measureShape = boost::make_shared<MeasureShapeT>(_policy.getPolicy("shape"));
        }
    }
    
    virtual ~MeasureSources() {
    }
    
    
    /// Return the Policy used to describe processing
    pexPolicy::Policy const& getPolicy() const { return _policy; }
    /// Return the log
    pexLog::Log &getLog() const { return *_moLog; }
    /// return the astrometric measurer
    typename MeasureAstrometryT::Ptr getMeasureAstrom() const { return _measureAstrom; }
    /// return the photometric measurer
    typename MeasurePhotometryT::Ptr getMeasurePhotom() const { return _measurePhotom; }
    /// return the shape measurer
    typename MeasureShapeT::Ptr getMeasureShape() const { return _measureShape; }

    virtual void measure(afwDet::Source& source, ExposureT const& exp) {
        CONST_PTR(afwImage::Wcs) wcs = exp.getWcs();
        _measure<SingleMeasurer<ExposureT> >(source, source, *wcs, exp, _policy);
    }
    virtual void measure(afwDet::Source& target, afwDet::Source const& source,
                         afwImage::Wcs const& wcs, ExposureT const& exp) {
        _measure<SingleMeasurer<ExposureT> >(target, source, wcs, exp, _policy);
    }
    virtual void measureGroup(afwDet::Source& target, afwDet::Source const& source,
                              afwImage::Wcs const& wcs, ExposureGroup<ExposureT>& group) {
        _measure<GroupMeasurer<ExposureT> >(target, source, wcs, group, _policy);
    }
    virtual void measureGroups(std::vector<afwDet::Source>& target, afwDet::Source const& source,
                               afwImage::Wcs const& wcs, std::vector<ExposureGroup<ExposureT> >& groups) {
        _measure<GroupsMeasurer<ExposureT> >(target, source, wcs, groups, _policy);
    }

private:
    pexPolicy::Policy _policy;   // Policy to describe processing

    PTR(pexLog::Log) _moLog; // log for measureObjects
    /*
     * Objects that know how to measure the object's properties
     */
    typename MeasureAstrometry<ExposureT>::Ptr _measureAstrom;
    typename MeasurePhotometry<ExposureT>::Ptr _measurePhotom;
    typename MeasureShape<ExposureT>::Ptr      _measureShape;

    /// Common driver function for measure, measureGroup, measureGroups
    template<class Measurer>
    void _measure(typename Measurer::SourceContainerT& target, afwDet::Source const& source, 
                  afwImage::Wcs const& wcs, typename Measurer::ExposureContainerT& exp,
                  pexPolicy::Policy const& policy) {
        typedef typename ExposureT::MaskedImageT MaskedImageT;

        Measurer::footprints(exp, source, wcs);
        Measurer::check(exp, target);

        // Centroids
        if (!getMeasureAstrom()) {
            Measurer::nullAstrom(target, source, exp);
        } else {
            PTR(afwDet::Astrometry) astrom = Measurer::measure(*getMeasureAstrom(), exp, source);
            // XXX record in target with setAstrom
            Measurer::extract(target, *astrom, policy, AstrometryExtractor());
            Measurer::astrom(target, source, exp);
        }

        // Shapes
        if (getMeasureShape()) {
            PTR(afwDet::Shape) shapes = Measurer::measure(*getMeasureShape(), exp, source);
            // XXX record in target with setShapes
            Measurer::extract(target, *shapes, policy, ShapeExtractor());
        }

        // Photometry
        if (getMeasurePhotom()) {
            PTR(afwDet::Photometry) phot = Measurer::measure(*getMeasurePhotom(), exp, source);
            // XXX record in target with setPhotometry
            Measurer::extract(target, *phot, policy, ApPhotExtractor());
            Measurer::extract(target, *phot, policy, PsfPhotExtractor());
            Measurer::extract(target, *phot, policy, ModelPhotExtractor());
            Measurer::extract(target, *phot, policy, InstPhotExtractor());
            Measurer::photom(target, *phot, policy);
        }
    }
};



/**
 * Factory functions to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
PTR(MeasureSources<ExposureT>) makeMeasureSources(ExposureT const& exp,
                                                  pexPolicy::Policy const& policy) {
    return boost::make_shared<MeasureSources<ExposureT> >(policy);
}

template<typename ExposureT>
PTR(MeasureSources<ExposureT>) makeMeasureSources(ExposureGroup<ExposureT> const& group,
                                                  pexPolicy::Policy const& policy) {
    return boost::make_shared<MeasureSources<ExposureT> >(policy);
}
template<typename ExposureT>
PTR(MeasureSources<ExposureT>) makeMeasureSources(std::vector<ExposureGroup<ExposureT> > const& groups,
                                                  pexPolicy::Policy const& policy) {
    return boost::make_shared<MeasureSources<ExposureT> >(policy);
}

       
}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
