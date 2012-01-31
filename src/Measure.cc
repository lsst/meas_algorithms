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
 
/// \file

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/Measure.h"

namespace pexPolicy = lsst::pex::policy;
namespace pexExceptions = lsst::pex::exceptions;    
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;


namespace lsst {
namespace meas {
namespace algorithms {

namespace {

/************************************************************************************************************/
/**
 * @brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintCentroid : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintCentroid(MaskedImageT const& mimage ///< The image the source lives in
                              ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                                  _n(0), _sum(0), _sumx(0), _sumy(0),
                                  _min( std::numeric_limits<double>::max()), _xmin(0), _ymin(0),
                                  _max(-std::numeric_limits<double>::max()), _xmax(0), _ymax(0),
                                  _bits(0) {}

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
    virtual void reset(afwDetection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

        _n++;
        _sum += val;
        _sumx += afwImage::indexToPosition(x)*val;
        _sumy += afwImage::indexToPosition(y)*val;
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
    PTR(afwDetection::Peak) makePeak(bool isNegative) const {
        return boost::make_shared<afwDetection::Peak>(isNegative ? afwDetection::Peak(_xmin, _ymin) :
                                                      afwDetection::Peak(_xmax, _ymax));
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

/************************************************************************************************************/
/**
 * @brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintBits : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintBits(MaskedImageT const& mimage ///< The image the source lives in
                              ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                                  _bits(0) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _bits = 0x0;
    }
    virtual void reset(afwDetection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        _bits |= loc.mask(0, 0);
    }

    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    typename MaskedImageT::Mask::Pixel _bits;
};


/// Check footprint for bad pixels
template<typename ExposureT>
void checkFootprint(
    afw::table::SourceRecord & source,
    MeasureSourcesFlags::KeyArray const & keys,
    ExposurePatch<ExposureT> const & patch,                  // Patch to check
    typename ExposureT::MaskedImageT::Mask::Pixel const bits // Bits in footprint
) {
    
    // Check for bits set in the Footprint
    if (bits & ExposureT::MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
        source.set(keys[MeasureSourcesFlags::EDGE], true);
    }
    if (bits & ExposureT::MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        source.set(keys[MeasureSourcesFlags::INTERPOLATED], true);
    }
    if (bits & ExposureT::MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        source.set(keys[MeasureSourcesFlags::SATURATED], true);
    }

    // Check for bits set near the centroid
    afwGeom::Point2D const & center = patch.getCenter(); // Center in appropriate coordinates
    afwGeom::Point2I llc(afwImage::positionToIndex(center.getX()) - 1,
                         afwImage::positionToIndex(center.getY()) - 1);
    afwDetection::Footprint const middle(afwGeom::BoxI(llc, afwGeom::ExtentI(3))); // central 3x3
    
    FootprintBits<typename ExposureT::MaskedImageT> bitsFunctor(patch.getExposure()->getMaskedImage());
    bitsFunctor.apply(middle);
    if (bitsFunctor.getBits() & ExposureT::MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        source.set(keys[MeasureSourcesFlags::INTERPOLATED_CENTER], true);
    }
    if (bitsFunctor.getBits() & ExposureT::MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        source.set(keys[MeasureSourcesFlags::SATURATED_CENTER], true);
    }
}

} // anonymous namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MeasureSources implementation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename ExposureT>
MeasureSources<ExposureT>::MeasureSources(pexPolicy::Policy const& policy) :
    _policy(policy),
    _log(pexLogging::Log::getDefaultLog().createChildLog("meas.algorithms.measureSource",
                                                         pexLogging::Log::INFO)),
    _schema(afw::table::SourceTable::makeMinimalSchema())
{
    _extendednessKey = _schema.addField<double>("extendedness", "FIXME! NEVER DOCUMENTED!");
    _flagKeys[MeasureSourcesFlags::EDGE] = _schema.addField<afw::table::Flag>(
        "flags.meas.edge", "source is in region labeled EDGE"
    );
    _flagKeys[MeasureSourcesFlags::INTERPOLATED] = _schema.addField<afw::table::Flag>(
        "flags.meas.interpolated.any", "source's footprint includes interpolated pixels"
    );
    _flagKeys[MeasureSourcesFlags::INTERPOLATED_CENTER] = _schema.addField<afw::table::Flag>(
        "flags.meas.interpolated.center", "source's center is close to interpolated pixels"
    );
    _flagKeys[MeasureSourcesFlags::SATURATED] = _schema.addField<afw::table::Flag>(
        "flags.meas.saturated.any", "source's footprint includes saturated pixels"
    );
    _flagKeys[MeasureSourcesFlags::SATURATED_CENTER] = _schema.addField<afw::table::Flag>(
        "flags.meas.saturated.center", "source's center is close to saturated pixels"
    );
    _flagKeys[MeasureSourcesFlags::PEAKCENTER] = _schema.addField<afw::table::Flag>(
        "flags.meas.peakcenter", "given center is position of peak pixel"
    );
}

template <typename ExposureT>
PTR(afw::table::SourceTable) MeasureSources<ExposureT>::makeTable() const {
    PTR(afw::table::SourceTable) table = afw::table::SourceTable::make(_schema);
    if (_policy.isString("source.centroid")) {
        table->defineCentroid(_policy.getString("source.centroid"));
    }
    if (_policy.isString("source.shape")) {
        table->defineShape(_policy.getString("source.shape"));
    }
    if (_policy.isString("source.apFlux")) {
        table->defineApFlux(_policy.getString("source.apFlux"));
    }
    if (_policy.isString("source.modelFlux")) {
        table->defineModelFlux(_policy.getString("source.modelFlux"));
    }
    if (_policy.isString("source.psfFlux")) {
        table->definePsfFlux(_policy.getString("source.psfFlux"));
    }
    if (_policy.isString("source.instFlux")) {
        table->defineInstFlux(_policy.getString("source.instFlux"));
    }
    return table;
}

template<typename ExposureT>
void MeasureSources<ExposureT>::apply(
    afw::table::SourceRecord & source,
    CONST_PTR(ExposureT) exp,
    afw::geom::Point2D const & center
) const {
    CONST_PTR(afwImage::Wcs) wcs = exp->getWcs();
    CONST_PTR(afwDetection::Footprint) foot = source.getFootprint();
    ExposurePatch<ExposureT> patch(exp, foot, center);
    FootprintCentroid<typename ExposureT::MaskedImageT> centroider(patch.getExposure()->getMaskedImage());
    centroider.apply(*patch.getFootprint());
    checkFootprint(source, _flagKeys, patch, centroider.getBits());
    // FIXME: should we be setting the patch's center here, or just going with the user input?
    _apply(source, patch);
}

template<typename ExposureT>
void MeasureSources<ExposureT>::apply(afw::table::SourceRecord & source, CONST_PTR(ExposureT) exp) const {
    CONST_PTR(afwImage::Wcs) wcs = exp->getWcs();
    CONST_PTR(afwDetection::Footprint) foot = source.getFootprint();
#if 0 // FIXME: where does/should the DETECT_NEGATIVE flag get set?
    bool negative = source.getFlagForDetection() & Flags::DETECT_NEGATIVE;
#endif
    // Get highest peak; 
    afwDetection::Footprint::PeakList const& peakList = foot->getPeaks();
    if (peakList.size() == 0) {
        throw LSST_EXCEPT(pexExceptions::RuntimeErrorException, 
                          (boost::format("No peak for source %d") % source.getId()).str());
    }
    PTR(afwDetection::Peak) peak = peakList[0];
    for (size_t i = 1; i < peakList.size(); ++i) {
        float value = peakList[i]->getPeakValue();
#if 0 // FIXME (see above)
        if (negative) {
            value *= -1;
        }
#endif
        if (value > peak->getPeakValue()) {
            peak = peakList[i];
        }
    }
    // set the initial centroid in the patch using the peak, then refine it.
    afwGeom::Point2D center(peak->getFx(), peak->getFy());
    ExposurePatch<ExposureT> patch(exp, foot, center);
    FootprintCentroid<typename ExposureT::MaskedImageT> centroider(patch.getExposure()->getMaskedImage());
    centroider.apply(*patch.getFootprint());
    patch.setCenter(afw::geom::Point2D(centroider.getX(), centroider.getY()));
    checkFootprint(source, _flagKeys, patch, centroider.getBits());
    _apply(source, patch);
}

template<typename ExposureT>
void MeasureSources<ExposureT>::_apply(
    afw::table::SourceRecord & source,
    ExposurePatch<ExposureT> & patch
) const {
    CONST_PTR(afw::table::SourceTable) table = source.getTable();
    bool gotCentroidSlot = false;
    for (
        typename std::list<PTR(Algorithm<ExposureT>)>::const_iterator i = _algorithms.begin();
        i != _algorithms.end();
        ++i
    ) {
        try {
            (**i).apply(source, patch);
        } catch (pex::exceptions::Exception const & e) {
            // Swallow all exceptions, because one bad measurement shouldn't affect all others
            _log->log(pex::logging::Log::DEBUG, boost::format("Measuring %s at (%d,%d): %s") %
                      (**i).getName() % source.getX() % source.getY() % e.what());
        } catch (...) {
            _log->log(pex::logging::Log::WARN, 
                      boost::format("Measuring %s at (%d,%d): Unknown non-LSST exception.") %
                      (**i).getName() % source.getX() % source.getY());
        }
        if (!gotCentroidSlot && table->getCentroidFlagKey().isValid() && source.getCentroidFlag()) {
            patch.setCenter(source.getCentroid());
            gotCentroidSlot = true;
        }
    }

    // FIXME: should be handled more gracefully; previously the policy defaults were paf dict,
    // but now they're in the pure-Python config...so if we're configured from c++ with
    // an incomplete policy, we don't have these params.  Best just not to fill the field, and
    // fix it better when we finish killing policy.
    if (_policy.isPolicy("classification")) {
        // Add some star/galaxy information.  The "extendedness" parameter is supposed to be the
        // probability of being extended
        std::vector<float> fac(3);// Fiddle factors for star/galaxy separation
        fac[0] = _policy.getDouble("classification.sg_fac1");
        fac[1] = _policy.getDouble("classification.sg_fac2");
        fac[2] = _policy.getDouble("classification.sg_fac3");
        bool const isStar = ((fac[0]*source.getInstFlux() + fac[1]*source.getInstFluxErr()) <
                             (source.getPsfFlux() + fac[2]*source.getPsfFluxErr()) ? 0.0 : 1.0);
        source[_extendednessKey] = (isStar ? 0.0 : 1.0);
    }
}

//
// Explicit instantiations
//
// \cond

#define INSTANTIATE(PIXEL) \
    template class MeasureSources<afwImage::Exposure<PIXEL> >; \
    template PTR(MeasureSources<afwImage::Exposure<PIXEL> >) \
        makeMeasureSources(afwImage::Exposure<PIXEL> const&, pexPolicy::Policy const&);

INSTANTIATE(float);
INSTANTIATE(double);

// \endcond
}}}
