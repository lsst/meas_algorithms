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

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/Measure.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

typedef std::pair<afw::table::KeyTuple<afw::table::Flux>,ScaledFlux::KeyTuple> ScaledFluxKeys;

template <typename T>
inline typename afw::table::Field<T>::Value get(
    afw::table::SourceRecord const & record,
    afw::table::Key<T> const & key,
    typename afw::table::Field<T>::Value default_ = false
) {
    return key.isValid() ? record.get(key) : default_;
}


void applyApCorr(
    afw::table::SourceRecord & source,
    ScaledFluxKeys const & keys,
    double apCorr, double apCorrErr, bool apCorrBad
) {
    if (apCorrBad) {
        if (keys.second.badCorrectionFlag.isValid()) {
            source.set(keys.second.badCorrectionFlag, true);
        }
        if (keys.first.flag.isValid()) {
            source.set(keys.first.flag, true);
        }
    }
    double flux = source.get(keys.first.meas);
    double fluxErr = source.get(keys.first.err);
    source.set(keys.first.meas, flux * apCorr);
    fluxErr *= apCorr;
    fluxErr *= fluxErr;
    flux *= apCorrErr;
    flux *= flux;
    source.set(keys.first.err, std::sqrt(fluxErr + flux));
}

} // anonymous

class MeasureSources::FluxCorrectionImpl {
public:
    ScaledFluxKeys canonical;
    std::vector<ScaledFluxKeys> others;
};

MeasureSourcesBuilder & MeasureSourcesBuilder::addAlgorithm(AlgorithmControl const & algorithmControl) {
    PTR(AlgorithmControl) p = algorithmControl.clone();
    p->name = _prefix + p->name;
    _ctrls.insert(p);
    return *this;
}

MeasureSourcesBuilder & MeasureSourcesBuilder::setCentroider(CentroidControl const & centroidControl) {
    PTR(CentroidControl) p = centroidControl.clone();
    p->name = _prefix + p->name;
    _centroider = p;
    return *this;
}

MeasureSources MeasureSourcesBuilder::build(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const & metadata,
    std::string const & canonicalFlux,
    int canonicalFluxIndex
) const {
    MeasureSources r;
    r._log.reset(pex::logging::Log::getDefaultLog().createChildLog(
                     "meas.algorithms.MeasureSource",
                     pex::logging::Log::INFO
                 ));
    AlgorithmControlMap ctrlMap;
    if (_centroider) {
        r._badCentroidKey = schema.addField<afw::table::Flag>(
            _prefix + "flags.badcentroid",
            "the centroid algorithm used to feed centers to other algorithms failed"
        );
        r._centroider = _centroider->makeAlgorithm(schema, metadata, ctrlMap);
        r._algorithms.push_back(r._centroider);
        ctrlMap[_centroider->name] = _centroider;
    }
    for (ControlSet::const_iterator i = _ctrls.begin(); i != _ctrls.end(); ++i) {
        r._algorithms.push_back((**i).makeAlgorithm(schema, metadata, ctrlMap));
        ctrlMap[(**i).name] = *i;
    }
    r._fluxCorrectionImpl.reset(new MeasureSources::FluxCorrectionImpl());
    for (
        MeasureSources::AlgorithmList::const_iterator i = r._algorithms.begin();
        i != r._algorithms.end();
        ++i
    ) {
        CONST_PTR(ScaledFlux) asScaledFlux = boost::dynamic_pointer_cast<ScaledFlux const>(*i);
        if ((**i).getControl().name == canonicalFlux) {
            if (!asScaledFlux) {
                throw LSST_EXCEPT(
                    pex::exceptions::InvalidParameterException,
                    (boost::format("Canonical flux (%s) is not an instance of ScaledFlux")
                     % canonicalFlux).str()
                );
            }
            int fluxCount = asScaledFlux->getFluxCount();
            if (canonicalFluxIndex < 0 || canonicalFluxIndex >= fluxCount) {
                throw LSST_EXCEPT(
                    pex::exceptions::InvalidParameterException,
                    (boost::format("Invalid index (%d) for canonical flux (must be between 0 and %d)")
                     % canonicalFluxIndex % (fluxCount - 1)).str()
                );
            }
            for (int i = 0; i < fluxCount; ++i) {
                if (i == canonicalFluxIndex) {
                    r._fluxCorrectionImpl->canonical =
                        ScaledFluxKeys(
                            asScaledFlux->getFluxKeys(i),
                            asScaledFlux->getFluxCorrectionKeys(i)
                        );
                } else {
                    r._fluxCorrectionImpl->others.push_back(
                        ScaledFluxKeys(
                            asScaledFlux->getFluxKeys(i),
                            asScaledFlux->getFluxCorrectionKeys(i)
                        )
                    );
                }
            }
        } else if (asScaledFlux) {
            int fluxCount = asScaledFlux->getFluxCount();
            for (int i = 0; i < fluxCount; ++i) {
                r._fluxCorrectionImpl->others.push_back(
                    ScaledFluxKeys(
                        asScaledFlux->getFluxKeys(i),
                        asScaledFlux->getFluxCorrectionKeys(i)
                    )
                );
            }
        }
    }
    return r;
}

MeasureSources::~MeasureSources() {}

void MeasureSources::correctFluxes(
    afw::table::SourceRecord & source,
    double apCorr,
    double apCorrErr,
    bool apCorrBad
) {
    if (!_fluxCorrectionImpl->canonical.first.meas.isValid()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot correct fluxes without defining a canonical flux measurement."
        );
    }
    bool badCanonicalFlux = get(source, _fluxCorrectionImpl->canonical.first.flag, false);
    bool badCanonicalPsfFactor = get(source, _fluxCorrectionImpl->canonical.second.psfFactorFlag, false);
    double canonicalPsfFactor = get(source, _fluxCorrectionImpl->canonical.second.psfFactor, 1.0);
    applyApCorr(source, _fluxCorrectionImpl->canonical, apCorr, apCorrErr, apCorrBad);
    for (
        std::vector<ScaledFluxKeys>::const_iterator i = _fluxCorrectionImpl->others.begin();
        i != _fluxCorrectionImpl->others.end();
        ++i
    ) {
        assert(i->first.meas.isValid());
        assert(i->first.err.isValid());
        if (badCanonicalFlux || badCanonicalPsfFactor) {
            // We can't do the correction because the canonical measurement failed.
            // We set "*.flags.badcorr" and the overall "*.flags" if we can.
            if (i->second.badCorrectionFlag.isValid()) {
                source.set(i->second.badCorrectionFlag, true);
            }
            if (i->first.flag.isValid()) {
                source.set(i->first.flag, true);
            }
        } else {
            if (i->second.psfFactorFlag.isValid() && source.get(i->second.psfFactorFlag)) {
                // Bad PSF factor for this measurement.  Mark it as failed overall if we can, then
                // let the uncorrected flux stand.
                if (i->first.flag.isValid()) {
                    source.set(i->first.flag, true);
                }
            } else {
                double psfFactor = get(source, i->second.psfFactor, 1.0);
                double scaling = canonicalPsfFactor / psfFactor;
                source[i->first.meas] *= scaling;
                source[i->first.err] *= scaling;
            }
            applyApCorr(source, *i, apCorr, apCorrErr, apCorrBad);
        }
    }
}

namespace {

/// Apply algorithm to measure source on image
///
/// Common code for executing an algorithm and catching exceptions.
template<typename PixelT>
void applyAlgorithm(
    Algorithm const& algorithm,         ///< Algorithm to apply
    afw::table::SourceRecord & source,  ///< Source to measure
    afw::image::Exposure<PixelT> const & exposure, ///< Exposure on which to measure
    afw::geom::Point2D const& center,              ///< Center to use
    PTR(pex::logging::Log) log,                     ///< Log for errors
    afw::table::SourceRecord const * reference = NULL, ///< Reference source for forced measurement.
    afw::geom::AffineTransform const * refToMeas = NULL ///< Transform from reference to measurement frame.
) {
    try {
        if (reference) {
            assert(refToMeas);
            algorithm.applyForced(source, exposure, center, *reference, *refToMeas);
        } else {
            algorithm.apply(source, exposure, center);
        }
    } catch (pex::exceptions::Exception const& e) {
        // Swallow all exceptions, because one bad measurement shouldn't affect all others
        log->log(pex::logging::Log::DEBUG, boost::format("Measuring %s on source %d at (%f,%f): %s") %
                 algorithm.getControl().name % source.getId() % center.getX() % center.getY() % e.what());
    } catch (...) {
        log->log(pex::logging::Log::WARN, 
                 boost::format("Measuring %s on source %d at (%f,%f): Unknown non-LSST exception.") %
                 algorithm.getControl().name % source.getId() % center.getX() % center.getY());
    }
}

} // namespace

template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::geom::Point2D const & center,
    bool refineCenter
) const {
    afw::geom::Point2D c(center);
    for (
        AlgorithmList::const_iterator i = _algorithms.begin();
        i != _algorithms.end();
        ++i
    ) {
        applyAlgorithm(**i, source, exposure, c, _log);
        if (refineCenter && *i == _centroider) { // should only match the first alg, but test is cheap
            if (source.get(_centroider->getKeys().flag)) {
                source.set(_badCentroidKey, true);
            } else {
                c = source.get(_centroider->getKeys().meas);
            }
        }
    }
}

template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure
) const {
    CONST_PTR(afw::detection::Footprint) foot = source.getFootprint();
    if (!foot) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                          (boost::format("No footprint for source %d") % source.getId()).str());
    }
    afw::detection::Footprint::PeakList const& peakList = foot->getPeaks();
    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    PTR(afw::detection::Peak) peak = peakList[0];
    // set the initial centroid in the patch using the peak, then refine it.
    afw::geom::Point2D center(peak->getFx(), peak->getFy());
    apply(source, exposure, center, true);
}


template <typename PixelT>
void MeasureSources::apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const & exposure,
    afw::table::SourceRecord const& reference,
    CONST_PTR(afw::image::Wcs) referenceWcs
) const {
    if (!referenceWcs) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException,
                          "referenceWcs argument is not optional");
    }

    CONST_PTR(afw::image::Wcs) wcs = exposure.getWcs();

    // Create a transformed footprint for the forced source.
    CONST_PTR(afw::detection::Footprint) refFoot = reference.getFootprint();
    if (!refFoot) {
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, 
                          (boost::format("No footprint for reference %d") % reference.getId()).str());
    }
    source.setFootprint(refFoot->transform(*referenceWcs, *wcs, exposure.getBBox()));

    // Compute the local transform from the reference frame to the measurement frame.
    afw::geom::AffineTransform refToSky = referenceWcs->linearizePixelToSky(reference.getCoord());
    afw::geom::AffineTransform skyToMeas = wcs->linearizeSkyToPixel(reference.getCoord());
    afw::geom::Point2D center = skyToMeas(reference.getCoord().getPosition(afw::geom::degrees));
    afw::geom::AffineTransform refToMeas = skyToMeas * refToSky;

    source.setCoord(reference.getCoord());
    for (
        AlgorithmList::const_iterator i = _algorithms.begin();
        i != _algorithms.end();
        ++i
    ) {
        applyAlgorithm(**i, source, exposure, center, _log, &reference, &refToMeas);
    }
}

#define INSTANTIATE(TYPE) \
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<TYPE> const &) const; \
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<TYPE> const &, \
                                    afw::geom::Point2D const &, bool) const; \
template void MeasureSources::applyWithCoord(afw::table::SourceRecord &, \
                                             afw::image::Exposure<TYPE> const &) const; \
template void MeasureSources::applyWithPixel(afw::table::SourceRecord &, \
                                             afw::image::Exposure<TYPE> const &) const; \
template void MeasureSources::apply(afw::table::SourceRecord &, afw::image::Exposure<TYPE> const &, \
                                    afw::table::SourceRecord const&, CONST_PTR(afw::image::Wcs)) const;

INSTANTIATE(float);
INSTANTIATE(double);

#undef INSTANTIATE

}}} // namespace lsst::meas::algorithms
