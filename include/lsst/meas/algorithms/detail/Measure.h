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
 
#if !defined(LSST_MEAS_ALGORITHMS_DETAIL_MEASURE_H)
#define LSST_MEAS_ALGORITHMS_DETAIL_MEASURE_H

#include <list>
#include <cmath>
#include "lsst/base.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Astrometry.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/Shape.h"
#include "lsst/meas/algorithms/MeasureQuantity.h"


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


namespace {
    /// Extractors to call the right extraction method 
    struct ApPhotExtractor {
        static std::string name = "source.apFlux";
        static void extract(afwDetection::Source source, afwDetection::Photometry phot) {
            source.extractApPhot(phot);
        }
    };
    struct PsfPhotExtractor {
        static std::string name = "source.psfFlux";
        static void extract(afwDetection::Source source, afwDetection::Photometry phot) {
            source.extractPsfPhot(phot);
        }
    };
    struct ModelPhotExtractor {
        static std::string name = "source.modelFlux";
        static void extract(afwDetection::Source source, afwDetection::Photometry phot) {
            source.extractModelPhot(phot);
        }
    };
    struct InstPhotExtractor {
        static std::string name = "source.instFlux";
        static void extract(afwDetection::Source source, afwDetection::Photometry phot) {
            source.extractInstPhot(phot);
        }
    };
    struct AstrometryExtractor {
        static std::string name = "source.astrom";
        static void extract(afwDetection::Source source, afwDetection::Astrometry astrom) {
            source.extractAstrometry(astrom);
        }
    };
    struct ShapeExtractor {
        static std::string name = "source.shape";
        static void extract(afwDetection::Source source, afwDetection::Shape shape) {
            source.extractShape(shape);
        }
    };

    /// Templated function to extract the correct measurement
    template<class MeasurementT, class Extractor>
    void extractMeasurements(afwDetection::Source& source,
                             afwDetection::Measurement<MeasurementT> const& measurements,
                             pexPolicy::Policy const& policy) {
        std::string const& name = Extractor::name;
        if (policy.isString(name)) {
            std::string const& alg = policy.getString(name);
            if (alg != "NONE") {
                afwDetection::Measurement<MeasurementT>::TPtr meas = measurements->find(alg);
                if (!meas) {
                    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                                      (boost::format("Can't find measurement from algorithm %s") % alg).str());
                }
                Extractor::extract(source, meas);
            }
        }
    }


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


/// How to make a peak
    template<typename ExposureT>
    struct SinglePeakMaker {
        static CONST_PTR(afwDetection::Peak) makePeak(ExposureT const& exp,
                                                      afwDet::Source const& source,
                                                      FootprintCentroid const& centroid) {
            return centroid.makePeak(source->getFlagForDetection() & Flags::DETECT_NEGATIVE);
        }
    };
    template<typename ExposureT>
    struct GroupPeakMaker {
        static CONST_PTR(afwDetection::Peak) makePeak(ExposureT const& exp,
                                                      afwDet::Source const& source,
                                                      FootprintCentroid const& centroid) {
            afwGeom::Point2D pix = exp.getWcs()->skyToPixel(source.getRaDec());
            return afwDet::Peak(pix.getX(), pix.getY());
        }
    };


    template<typename ExposureT, class PeakMaker>
    void checkPixels(ExposurePatch<ExposureT>& patch,
                     afwDet::Source const& source) {
        CONST_PTR(ExposureT) exp = patch.getExposure();

        FootprintCentroid centroidFunctor(exp->getMaskedImage());
        centroidFunctor.apply(*foot);

        CONST_PTR(afwDetection::Peak) peak = PeakMaker<ExposureT>::makePeak(exp, source, centroidFunctor);
        patch.setPeak(peak);

        // Check for bits set in the Footprint
        if (centroidFunctor.getBits() & exp->getMaskedImage()->getPlaneBitMask("EDGE")) {
            orFlag(ExposurePatch::EDGE);
        }
        if (centroidFunctor.getBits() & exp->getMaskedImage()->getPlaneBitMask("INTRP")) {
            orFlag(ExposurePatch::INTERP);
        }
        if (centroidFunctor.getBits() & exp->getMaskedImage()->getPlaneBitMask("SAT")) {
            orFlag(ExposurePatch::SAT);
        }

        // Check for bits set near the centroid
        afwGeom::Point2I llc(afwImage::positionToIndex(peak->xf) - 1,
                             afwImage::positionToIndex(peak->yf) - 1);
        afwDetection::Footprint const middle(afwGeom::BoxI(llc, afwGeom::ExtentI(3))); // central 3x3
        centroidFunctor.apply(middle);
        if (centroidFunctor.getBits() & exp->getMaskedImage()->getPlaneBitMask("INTRP")) {
            orFlag(ExposurePatch::INTERP_CENTER);
        }
        if (centroidFunctor.getBits() & exp->getMaskedImage()->getPlaneBitMask("SAT")) {
            orFlag(ExposurePatch::SAT_CENTER);
        }
    }


/// Measurers for use with _measure
///
/// These make up for the lack of support for passing method names: we simply
/// template on the Measurer class, and call a class (static) method of a
/// standard name.  The optimised compiled code should end up calling the
/// appropriate method name directly.
template<typename ExposureT>
struct OneMeasurer {
    typedef ExposureT ExposureContainerT;
    typedef afwDet::Source SourceContainerT;
    static void footprints(ExposurePatch<ExposureT>& exp,
                           afwDet::Source const& source,
                           afwImage::Wcs const& sourceWcs) {
        exp.setFootprint(source.getFootprint());
    }
    static void check(afwDet::Source &source, ExposurePatch<ExposureT>& exp) {
        return checkPixels(exp, source);
    }
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(ExposurePatch<ExposureT> const& exp,
                                     afwDet::Source const& source) {
        return alg->measureOne(exp, source);
    }
    template<typename MeasurementT, class Extractor>
    static void extract(afwDet::Source& target, MeasurementT const& meas, pexPolicy const& policy) {
        extractMeasurements<MeasurementT, Extractor>(target, meas, policy);
    }
    static void nullAstrom(afwDet::Source& target, afwDet::Source const& source,
                           ExposurePatch<ExposureT> const& exp) {
        CONST_PTR(afwDet::Peak) peak = exp.getPeak();
        target->setXAstrom(peak->getIx());
        target->setYAstrom(peak->getIy());
        target->setFlagForDetection(target->getFlagForDetection() | Flags::PEAKCENTER);
    }
    static void astrom(afwDet::Source& target, afwDet::Source const& source,
                       ExposurePatch<ExposureT> const& exp) {
        if (lsst::utils::isnan(target.getXAstrom()) || lsst::utils::isnan(target.getYAstrom())) {
            CONST_PTR(afwDet::Peak) peak = exp.getPeak();
            target->setXAstrom(peak->getFx());
            target->setYAstrom(peak->getFy());
            target->setFlagForDetection(src->getFlagForDetection() | Flags::PEAKCENTER);
        }
    }
    static void photom(afwDet::Source& target, afwDetection::Photometry const& phot,
                       pexPolicy::Policy const& policy) {
        // Set photometry flags
        boost::int64_t flag = target->getFlagForDetection();
        for (afwDetection::Measurement<afwDetection::Photometry>::const_iterator i = photom->begin();
            i != fluxes->end(); ++i) {
            flag |= (*i)->getFlag();
        }
        target->setFlagForDetection(flag);

        // Add some star/galaxy information.  The "extendedness" parameter is supposed to be the
        // probability of being extended
        std::vector<float> fac(3);// Fiddle factors for star/galaxy separation
        fac[0] = getNumeric(policy, "classification.sg_fac1");
        fac[1] = getNumeric(policy, "classification.sg_fac2");
        fac[2] = getNumeric(policy, "classification.sg_fac3");

        bool const isStar = ((fac[0]*target->getInstFlux() + fac[1]*target->getInstFluxErr()) <
                             (target->getPsfFlux() + fac[2]*target->getPsfFluxErr()) ? 0.0 : 1.0);
#if 0
        target->setExtendedness(isStar ? 0.0 : 1.0);
#else
        target->setApDia(isStar ? 0.0 : 1.0);
#endif
    }
};
template<typename ExposureT>
struct GroupMeasurer {
    typedef ExposureGroup<ExposureT> ExposureContainerT;
    typedef afwDet::Source SourceContainerT;
    static void footprints(ExposureGroup<ExposureT>& group,
                           afwDet::Source const& source,
                           afwImage::Wcs const& sourceWcs) {
        group.setFootprints(source.getFootprint(), sourceWcs);
    }
    static void check(ExposureGroup<ExposureT>& group,
                      afwDet::Source const& source) {
        for (typename ExposureGroup<ExposureT>::iterator iter = group.begin(); iter != group.end(); ++iter) {
            checkPixels(*iter, source);
        }
    }
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                     ExposureGroup<ExposureT> const& group,
                                     afwDet::Source const& source) {
        return alg->measureGroup(exp, source);
    }
    static void nullAstrom(afwDet::Source& target, afwDet::Source const& source,
                           ExposureGroup<ExposureT> const& group) {
        target->setXAstrom(source->setXAstrom());
        target->setYAstrom(source->setYAstrom());
        target->setFlagForDetection(target->getFlagForDetection() | Flags::PEAKCENTER);
    }
    static void astrom(afwDet::Source& target, afwDet::Source const& source,
                       ExposureGroup<ExposureT> const& exp) {
        if (lsst::utils::isnan(target.getXAstrom()) || lsst::utils::isnan(target.getYAstrom())) {
            target->setXAstrom(source->getXAstrom());
            target->setYAstrom(source->getYAstrom());
            target->setFlagForDetection(src->getFlagForDetection() | Flags::PEAKCENTER);
        }
    }
    static void photom(afwDet::Source& target, afwDetection::Photometry const& phot,
                       pexPolicy::Policy const& policy) {
        SingleMeasurer<MeasurementT, ExposureT>::photFlags(target, phot, policy);
    }
};
template<typename ExposureT>
struct GroupsMeasurer {
    typedef std::vector<ExposureGroup<ExposureT> > ExposureContainerT;
    typedef std::vector<afwDet::Source> SourceContainerT;
    static void footprints(std::vector<ExposureGroup<ExposureT> >& groups,
                           afwDet::Source const& source,
                           afwImage::Wcs const& sourceWcs) {
        for (typename std::vector<ExposureGroup<ExposureT> >::iterator iter = groups.begin(); 
             iter != groups.end(); ++ iter) {
            iter->setFootprints(source.getFootprint(), sourceWcs);
        }
    }
    static void check(std::vector<ExposureGroup<ExposureT> >& groups,
                      std::vector<afwDet::Source> const& sources) {
        typename std::vector<afwDet::Source>::iterator source = sources.begin();
        for (typename std::vector<ExposureGroup<ExposureT> >::iterator group = groups.begin();
             group != groups.end(); ++group, ++source) {
            checkPixels(*group, *source);
        }
    }
    template<typename MeasurementT>
    static PTR(MeasurementT) measure(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                     std::vector<ExposureGroup<ExposureT> > const& groups,
                                     afwDet::Source const& source) {
        return alg->measureGroups(exp, source);
    }
    static void nullAstrom(std::vector<afwDet::Source>& target, afwDet::Source const& source,
                           std::vector<ExposureGroup<ExposureT> > const& groups) {
        typename std::vector<afwDet::Source>::iterator target = targets.begin();
        for (typename std::vector<ExposureGroup<ExposureT> >::iterator group = groups.begin();
             group != groups.end(); ++group, ++source) {
            GroupMeasurer::nullAstrom(*target, source, *group);
        }
    }
    static void astrom(std::vector<afwDet::Source>& targets, afwDet::Source const& source,
                       std::vector<ExposureGroup<ExposureT> > const& groups) {
        typename std::vector<afwDet::Source>::iterator target = targets.begin();
        for (typename std::vector<ExposureGroup<ExposureT> >::iterator group = groups.begin();
             group != groups.end(); ++group, ++source) {
            GroupMeasurer::astrom(*target, source, *group);
        }
    }
    static void photom(std::vector<afwDet::Source>& targets, afwDetection::Photometry const& phots,
                       pexPolicy::Policy const& policy) {
        typename afwDetection::Measurement<afwDetection::Photometry>::iterator phot = phots.begin();
        for (typename std::vector<afwDet::Source>::iterator target = targets.begin(); 
             target != targets.end(); ++target, ++phot) {
            GroupMeasurer<ExposureT>::photFlags(target, phot, policy);
        }
    }
};




} // anonymous namespace

        

/// Common driver function for measureOne, measureGroup, measureGroups
template<typename ExposureT, typename Measurer>
virtual void _measure(Measurer::SourceContainerT& target, afwDet::Source const& src, 
                      afwImage::Wcs const& wcs, Measurer::ExposureContainerT const& exp,
                      pexPolicy::Policy const& policy)
{
    typedef typename ExposureT::MaskedImageT MaskedImageT;

    Measurer::footprints(exp, source, wcs);
    Measurer::check(exp, target);

    // Centroids
    if (!getMeasureAstrom()) {
        Measurer::nullAstrom(target, exp);
    } else {
        PTR(afwDet::Astrometry) astrom = Measurer::measure<afwDetection::Astrometry>(exp, source);
        Measurer::extract<afwDetection::Astrometry, AstrometryExtractor>(target, astrom, policy);
        Measurer::astrom(target, src, exp);
    }

    // Shapes
    if (getMeasureShape()) {
        PTR(afwDet::Shape) shapes = Measurer::measure<afwDetection::Shape>(exp, source);
        Measurer::extract<afwDetection::Shape, ShapeExtractor>(target, shapes, policy);
    }

    // Photometry
    if (getMeasurePhotom()) {
        PTR(afwDet::Photometry) phot = Measurer::measure<afwDetection::Photometry>(exp, source);
        Measurer::extract<afwDetection::Photometry, ApPhotExtractor>(target, phot, policy);
        Measurer::extract<afwDetection::Photometry, PsfPhotExtractor>(target, phot, policy);
        Measurer::extract<afwDetection::Photometry, ModelPhotExtractor>(target, phot, policy);
        Measurer::extract<afwDetection::Photometry, InstPhotExtractor>(target, phot, policy);
        Measurer::photom(target, phot, policy);
    }
}

}}} // lsst::meas::algorithms
#endif
