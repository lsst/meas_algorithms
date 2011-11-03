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
#include "boost/noncopyable.hpp"
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

/// High-level class to perform source measurement
///
/// Iterates over the various measurement types (Astrometry, Shape, Photometry).
template<typename ExposureT>
class MeasureSources : private boost::noncopyable {
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
            _measureAstrom = boost::make_shared<MeasureAstrometryT>(*_policy.getPolicy("astrometry"));
        }
        
        if (_policy.isPolicy("photometry")) {
            _measurePhotom = boost::make_shared<MeasurePhotometryT>(*_policy.getPolicy("photometry"));
        }

        if (_policy.isPolicy("shape")) {
            _measureShape = boost::make_shared<MeasureShapeT>(*_policy.getPolicy("shape"));
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

    virtual void measure(afwDet::Source& target, CONST_PTR(ExposureT) exp) {
        CONST_PTR(afwImage::Wcs) wcs = exp->getWcs();
        CONST_PTR(afwDet::Footprint) foot = target.getFootprint();
        bool negative = target.getFlagForDetection() & Flags::DETECT_NEGATIVE;
        // Get highest peak
        afwDet::Footprint::PeakList const& peakList = foot->getPeaks();
        if (peakList.size() == 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, 
                              (boost::format("No peak for source %d") % target.getId()).str());
        }
        PTR(afwDet::Peak) peak = peakList[0];
        for (size_t i = 1; i < peakList.size(); ++i) {
            float value = peakList[i]->getPeakValue();
            if (negative) {
                value *= -1;
            }
            if (value > peak->getPeakValue()) {
                peak = peakList[i];
            }
        }
        afwGeom::Point2D center(peak->getFx(), peak->getFy());
        ExposurePatch<ExposureT> patch(exp, foot, center);
        _measure<detail::SingleMeasurer<ExposureT> >(target, target, *wcs, patch);
    }
    virtual void measure(afwDet::Source& target, CONST_PTR(ExposureT) exp, afwGeom::Point2D const& center) {
        CONST_PTR(afwImage::Wcs) wcs = exp->getWcs();
        ExposurePatch<ExposureT> patch(exp, target.getFootprint(), center);
        _measure<detail::SingleMeasurer<ExposureT> >(target, target, *wcs, patch);
    }
    virtual void measure(afwDet::Source& target, afwDet::Source const& source,
                         afwImage::Wcs const& wcs, CONST_PTR(ExposureT) exp) {
        std::vector<CONST_PTR(ExposureT)> exposures(1);
        exposures[0] = exp;
        measure(target, source, wcs, exposures);
    }
    virtual void measure(afwDet::Source& target, afwDet::Source const& source,
                         afwImage::Wcs const& wcs, std::vector<CONST_PTR(ExposureT)> const& exposures) {
        size_t size = exposures.size();
        std::vector<PTR(ExposurePatch<ExposureT>)> patches(size);
        afwGeom::Point2D center(source.getXAstrom(), source.getYAstrom());
        for (size_t i = 0; i != size; ++i) {
            CONST_PTR(afwImage::Wcs) expWcs = exposures[i]->getWcs();
            if (!expWcs) {
                throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, 
                                  (boost::format("No WCS for exposure %d") % i).str());
            }
            afwDet::Footprint::Ptr foot = source.getFootprint()->transform(wcs, *expWcs,
                                                                           exposures[i]->getBBox());
            patches[i] = makeExposurePatch(exposures[i], foot, center, wcs);
        }
        _measure<detail::MultipleMeasurer<ExposureT> >(target, const_cast<afwDet::Source&>(source), 
                                                       wcs, patches);
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

    /// Common driver function for measure
    template<class Measurer>
    void _measure(afwDet::Source& target,
                  afwDet::Source& source, 
                  afwImage::Wcs const& wcs,
                  typename Measurer::ExposureContainerT& patches
        );
#if 0
{
        //std::cerr << "Measuring source " << source.getId() << std::endl;

        Measurer::footprints(exp, source, wcs);
        Measurer::check(exp, source);

        // Centroids
        if (!getMeasureAstrom()) {
            Measurer::nullAstrom(target, source, exp);
        } else {
            PTR(MeasureQuantity<afwDet::Astrometry, ExposureT>) meas = getMeasureAstrom();
            PTR(afwDet::Astrometry) astrom = 
                Measurer::template measure<afwDet::Astrometry>(meas, exp, source);
            Measurer::template set<detail::AstrometrySetter>(target, astrom);
            Measurer::template extract<detail::AstrometryExtractor>(target, _policy);
            Measurer::astrom(target, source, exp);
        }

        // Shapes
        if (getMeasureShape()) {
            PTR(afwDet::Shape) shapes = 
                Measurer::template measure<afwDet::Shape>(getMeasureShape(), exp, source);
            // XXX record in target with setShapes
            Measurer::template set<detail::ShapeSetter>(target, shapes);
            Measurer::template extract<detail::ShapeExtractor>(target, _policy);
        }

        // Photometry
        if (getMeasurePhotom()) {
            PTR(afwDet::Photometry) phot = 
                Measurer::template measure<afwDet::Photometry>(getMeasurePhotom(), exp, source);
            // XXX record in target with setPhotometry
            Measurer::template set<detail::PhotometrySetter>(target, phot);
            Measurer::template extract<detail::ApPhotExtractor>(target, _policy);
            Measurer::template extract<detail::PsfPhotExtractor>(target, _policy);
            Measurer::template extract<detail::ModelPhotExtractor>(target, _policy);
            Measurer::template extract<detail::InstPhotExtractor>(target, _policy);
            Measurer::photom(target, *phot, _policy);
        }

        Measurer::flags(target, exp);
    }
#endif
};



/**
 * Factory functions to return a MeasureSources of the correct type (cf. std::make_pair)
 */
template<typename ExposureT>
PTR(MeasureSources<ExposureT>) makeMeasureSources(ExposureT const& exp,
                                                  pexPolicy::Policy const& policy) {
    return boost::make_shared<MeasureSources<ExposureT> >(policy);
}
       
}}}
#endif // LSST_MEAS_ALGORITHMS_MEASURE_H
