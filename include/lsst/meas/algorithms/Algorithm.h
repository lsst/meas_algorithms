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
 
#if !defined(LSST_MEAS_ALGORITHMS_ALGORITHM_H)
#define LSST_MEAS_ALGORITHMS_ALGORITHM_H

#include "boost/noncopyable.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/policy.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Filter.h"

namespace pexLog = lsst::pex::logging;
namespace pexPolicy = lsst::pex::policy;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

namespace lsst { namespace meas { namespace algorithms {

/// A convenience container for the exposure, peak and footprint that will be measured.
///
/// This is more useful than a std::pair or similar.
template<typename ExposureT>
class ExposurePatch : private boost::noncopyable {
public:
    typedef unsigned char FlagT;
    typedef PTR(ExposurePatch) Ptr;
    typedef CONST_PTR(ExposurePatch) ConstPtr;

    /// Flag values, indicating which measurement is bad
    enum { NONE           = 0x00,     /// None bad
           EDGE           = 0x01,     /// Footprint overlaps an edge
           INTERP         = 0x02,     /// Footprint includes interpolated pixels
           INTERP_CENTER  = 0x04,     /// Peak pixel is interpolated
           SAT            = 0x08,     /// Footprint includes saturated pixels
           SAT_CENTER     = 0x10,     /// Peak pixel is saturated
           ASTROMETRY     = 0x20,     /// Bad astrometry
           SHAPE          = 0x40,     /// Bad shapes
           PHOTOMETRY     = 0x80,     /// Bad photometry
           ALL            = 0xFF      /// All are bad
    };

    /// Constructor
    explicit ExposurePatch(CONST_PTR(ExposureT) exp, 
                           CONST_PTR(afwDet::Footprint) foot=CONST_PTR(afwDet::Footprint)(),
                           CONST_PTR(afwDet::Peak) peak=CONST_PTR(afwDet::Peak)(),
                           FlagT flags=NONE) :
        _exp(exp), _foot(foot), _peak(peak), _flags(flags) {}
    explicit ExposurePatch(CONST_PTR(ExposureT) exp, 
                           CONST_PTR(afwDet::Peak) peak,
                           CONST_PTR(afwDet::Footprint) foot=CONST_PTR(afwDet::Footprint)(),
                           FlagT flags=NONE) :
        _exp(exp), _foot(foot), _peak(peak), _flags(flags) {}

    /// Accessors
    CONST_PTR(ExposureT) getExposure() const { return _exp; }
    CONST_PTR(afwDet::Footprint) getFootprint() const { return _foot; }
    CONST_PTR(afwDet::Peak) getPeak() const { return _peak; }
    bool getFlags() const { return _flags; }

    /// Modifiers
    void setExposure(CONST_PTR(ExposureT) exp) { _exp = exp; }
    void setFootprint(CONST_PTR(afwDet::Footprint) foot) { _foot = foot; }
    void setPeak(CONST_PTR(afwDet::Peak) peak) { _peak = peak; }
    void setFlags(FlagT flags) { _flags = flags; }
    void orFlag(FlagT flags) { _flags |= flags; }

private:
    CONST_PTR(ExposureT) _exp;          // Exposure to be measured
    CONST_PTR(afwDet::Footprint) _foot; // Footprint to be measured, or NULL
    CONST_PTR(afwDet::Peak) _peak;      // Peak being measured, or NULL
    FlagT _flags;                       // Flags indicating which measurement is bad
};

/// Factory functions for ExposurePatch
template<typename ExposureT>
PTR(ExposurePatch<ExposureT>) makeExposurePatch(
    CONST_PTR(ExposureT) exp,
    CONST_PTR(afwDet::Peak) peak,
    CONST_PTR(afwDet::Footprint) foot=CONST_PTR(afwDet::Footprint)()) {
    return boost::make_shared<ExposurePatch<ExposureT> >(exp, peak, foot);
}    
template<typename ExposureT>
PTR(ExposurePatch<ExposureT>) makeExposurePatch(
    CONST_PTR(ExposureT) exp,
    CONST_PTR(afwDet::Footprint) foot=CONST_PTR(afwDet::Footprint)(),
    CONST_PTR(afwDet::Peak) peak=CONST_PTR(afwDet::Peak)()) {
    return boost::make_shared<ExposurePatch<ExposureT> >(exp, foot, peak);
}

/// A group of exposures to be measured.
///
/// The idea behind a "group" is that they share some quality in common, so that
/// algorithms can assume that some characteristics don't change across the
/// group.  Here we assume that the "group" is defined by the filter (e.g.,
/// galaxy shapes shouldn't change much for images taken with the same filter),
/// so it ensures that all inputs have the same filter.  We could make this
/// class more general by templating on the quality that is held constant.
///
/// Could make this inherit from std::vector, but I've read that std::vector
/// isn't supposed to be used as a base class (no virtual destructor).
template<typename ExposureT>
class ExposureGroup : private boost::noncopyable {
public:
    typedef ExposurePatch<ExposureT> PatchT;
    typedef std::vector<PTR(PatchT)> PatchSetT;
    typedef typename PatchSetT::iterator iterator;
    typedef typename PatchSetT::const_iterator const_iterator;
    typedef PTR(ExposureGroup) Ptr;
    typedef CONST_PTR(ExposureGroup) ConstPtr;

    /// Constructor
    ///
    /// Don't want a constructor like ExposureGroup(PatchSetT& patches) because
    /// what happens if the exposures don't all have the same filter?  We should
    /// avoid throwing an exception in the constructor.
    ExposureGroup() : _patches(), _filter() {}

    /// Iterators
    ///
    /// Access to the exposures is through the iterators.
    iterator begin() { return _patches.begin(); }
    const_iterator begin() const { return _patches.begin(); }
    iterator end() { return _patches.end(); }
    const_iterator end() const { return _patches.end(); }

    /// Accessors
    size_t size() const { return _patches.size(); }
    PTR(PatchT) operator[](size_t index) const { return _patches[index]; }
    afwImage::Filter getFilter() const { return _filter; }

    /// Append to the list of exposure patches.
    void append(PTR(PatchT) patch) {
        afwImage::Filter const& patchFilter = patch->getExposure()->getFilter();
        if (_filter.getId() == afwImage::Filter::UNKNOWN) {
            _filter = afwImage::Filter(patchFilter.getId());
        } else if (_filter.getId() != patchFilter.getId()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                              (boost::format("Exposure filter (%d) doesn't match group (%d)") %
                               patchFilter.getId() % _filter.getId()).str());
        }
        _patches.push_back(patch);
    }

    /// Set footprints for all exposures
    void setFootprints(afwDet::Footprint const& foot, // Original (source) footprint
                       afwImage::Wcs const& wcs       // Original (source) WCS
        ) {
        for (iterator iter = begin(); iter != end(); ++iter) {
            PTR(ExposurePatch<ExposureT>) patch = *iter;
            CONST_PTR(ExposureT) exp = patch->getExposure();
            PTR(afwDet::Footprint) targetFoot = foot.transform(wcs, *exp->getWcs(), exp->getBBox());
            patch->setFootprint(targetFoot);
        }
    }

private:
    PatchSetT _patches;                 // Exposure patches of interest
    afwImage::Filter _filter;           // Common filter for patches
};

/// Factory function for ExposureGroup
template<typename ExposureT>
PTR(ExposureGroup<ExposureT>) makeExposureGroup(PTR(ExposurePatch<ExposureT>) patch) {
    PTR(ExposureGroup<ExposureT>) group = boost::make_shared<ExposureGroup<ExposureT> >();
    group->append(patch);
    return group;
}    

/// Base class for algorithms for measuring MeasurementT (e.g., Photometry)
template<typename MeasurementT, typename ExposureT>
class Algorithm : private boost::noncopyable {
public:
    typedef ExposurePatch<ExposureT> PatchT;
    typedef ExposureGroup<ExposureT> GroupT;
    typedef std::vector<PTR(ExposureGroup<ExposureT>)> GroupSetT;

    /// Constructor
    Algorithm() {}

    /// Destructor
    virtual ~Algorithm() {}

    /// Measure a single value from a single image.
    ///
    /// Returns leaf of MeasurementT (single measurement).
    ///
    /// Pure-virtual, so subclasses MUST define: it is the essence of the
    /// measurement, as the other measure functions can (but need not) be
    /// defined from it.
    virtual PTR(MeasurementT) measureOne(PatchT const&, afwDet::Source const&) const = 0;
    
    /// Measure a values from a group of images.
    ///
    /// Returns composite of MeasurementT (one measurement per exposure).
    ///
    /// Because it is a 'group' of images (images in the same filter), we can
    /// assume they share some characteristics (e.g., center, shape).
    ///
    /// Defaults to iteratively calling measureOne. However, if the measurement
    /// cannot be obtained by merely averaging the outputs of a single
    /// measurement, e.g., some measured parameters are made across all
    /// exposures as part of the measurement (e.g., a common shape), then the
    /// Algorithm needs to define this method.
    virtual PTR(MeasurementT) measureGroup(GroupT const& group,
                                           afwDet::Source const& source) const {
        PTR(MeasurementT) meas(new MeasurementT());
        for (typename GroupT::const_iterator iter = group.begin(); iter != group.end(); ++iter) {
            PTR(MeasurementT) val;
            try {
                PTR(PatchT) patch = *iter;
                val = measureOne(*patch, source);
            } catch (lsst::pex::exceptions::Exception const& e) {
#if 0
                std::cerr << (boost::format("Measuring single %s at (%d,%d): %s") %
                              getName() % source.getXAstrom() % source.getYAstrom() %
                              e.what()).str() << std::endl;
#endif
                val = measureNull();
            }
            val->setAlgorithm(getName());
            meas->add(val);
        }
        return meas->average();
    }
    
    /// Measure multiple values from groups of multiple images.
    /// Returns composite of MeasurementT (one measurement per group).
    /// Defaults to iteratively calling measureGroup and then averaging the results for each group.
    /// However, if the measurement requires linking a measurement across the groups (e.g., a fitted center
    /// from treating all the data), then the Algorithm needs to define this method.
    virtual PTR(MeasurementT) measureGroups(GroupSetT const& groups, afwDet::Source const& source) const {
        PTR(MeasurementT) meas(new MeasurementT());
        for (typename GroupSetT::const_iterator iter = groups.begin(); iter != groups.end(); ++iter) {
            PTR(MeasurementT) val;
            try {
                PTR(GroupT) group = *iter;
                val = measureGroup(*group, source);
            } catch (lsst::pex::exceptions::Exception const& e) {
#if 0
                std::cerr << (boost::format("Measuring single %s at (%d,%d): %s") %
                              getName() % source.getXAstrom() % source.getYAstrom() % 
                              e.what()).str() << std::endl;
#endif
                val = measureNull();
            }
            val->setAlgorithm(getName());
            meas->add(val);
        }
        return meas;
    }

    /// Return a null measurement
    ///
    /// This is called when we hit an exception.
    virtual PTR(MeasurementT) measureNull(void) const {
        return MeasurementT::null();
    }
    
    /// Configure the algorithm
    virtual void configure(pexPolicy::Policy const&) {};

    /// Name of the algorithm
    virtual std::string getName() const = 0;

    /// Clone algorithm
    virtual PTR(Algorithm<MeasurementT, ExposureT>) clone() const = 0;
};


}}} // namespace

#endif
