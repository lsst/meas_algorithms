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

#include "lsst/base.h"
#include "lsst/pex/policy.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/image/Filter.h"

namespace lsst { namespace meas { namespace algorithms {

namespace pexPolicy = lsst::pex::policy;
namespace afwDet = lsst::afw::detection;
namespace afwImage = lsst::afw::image;

/// A convenience container for the exposure and footprint that will be measured.
///
/// This is clearer than a std::pair or similar.
template<typename ExposureT>
class ExposurePatch {
    /// Flag values, indicating which measurement is bad
    enum { NONE       = 0x00,           /// None bad
           ASTROMETRY = 0x01,           /// Bad astrometry
           SHAPE      = 0x02,           /// Bad shapes
           PHOTOMETRY = 0x04,           /// Bad photometry
           ALL        = 0xFF            /// All are bad
    };
public:
    typedef unsigned char FlagT;

    /// Constructor
    ExposurePatch(CONST_PTR(ExposureT) exp, 
                  CONST_PTR(afwDet::Footprint) foot=CONST_PTR(afwDet::Footprint)(),
                  FlagT flags=NONE) :
        _exp(exp), _foot(foot), _flags(flags) {}

    /// Destructor
    ~ExposurePatch();

    /// Accessors
    CONST_PTR(ExposureT) getExposure() const { return _exp; }
    CONST_PTR(afwDet::Footprint) getFootprint() const { return _foot; }
    bool getFlags() const { return _flags; }

    /// Modifiers
    void setExposure(CONST_PTR(ExposureT) exp) { _exp = exp; }
    void setFootprint(CONST_PTR(afwDet::Footprint) foot) { _foot = foot; }
    void setFlags(FlagT flags) { _flags = flags; }    

private:
    CONST_PTR(ExposureT) _exp;          // Exposure to be measured
    CONST_PTR(afwDet::Footprint) _foot; // Footprint to be measured, or NULL
    FlagT _flags;                       // Flags indicating which measurement is bad
};

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
class ExposureGroup {
public:
    typedef ExposurePatch<ExposureT> PatchT;
    typedef typename std::vector<PatchT>::iterator iterator;
    typedef typename std::vector<PatchT>::const_iterator const_iterator;

    /// Constructor
    ///
    /// Don't want a constructor like ExposureGroup(std::vector<PatchT>
    /// &patches) because what happens if the exposures don't all have the same
    /// filter?  We should avoid throwing an exception in the constructor.
    ExposureGroup() : _patches(), _filter() {}

    /// Destructor
    ~ExposureGroup() {}

    /// Iterators
    ///
    /// Access to the exposures is through the iterators.
    iterator begin() { return _patches.begin(); }
    const_iterator begin() const { return _patches.begin(); }
    iterator end() { return _patches.end(); }
    const_iterator end() const { return _patches.end(); }

    /// Accessors
    CONST_PTR(afwImage::Filter) getFilter() const { return _filter; }

    /// Append to the list of exposure patches.
    void append(PatchT& patch) {
        CONST_PTR(afwImage::Filter) patchFilter = patch.getExposure()->getFilter();
        if (_filter->getId() == afwImage::Filter::UNKNOWN) {
            *_filter = *patchFilter;
        } else if (_filter->getFilterProperty() != patchFilter->getFilterProperty()) {
            // XXX Throw exception
            abort();
        }
        _patches.push_back(patch);
    }

private:
    std::vector<PatchT> _patches;       // Exposure patches of interest
    CONST_PTR(afwImage::Filter) _filter; // Common filter for patches
};


/// Base class for algorithms for measuring MeasurementT (e.g., Photometry)
template<typename MeasurementT, typename ExposureT>
class Algorithm {
public:
    typedef ExposurePatch<ExposureT> PatchT;
    typedef ExposureGroup<ExposureT> GroupT;
    typedef std::vector<ExposureGroup<ExposureT> > GroupSetT;

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
    virtual PTR(MeasurementT) measureOne(PatchT const&, afwDet::Source const&) = 0;
    
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
                                           afwDet::Source const& source) {
        PTR(MeasurementT) meas(new MeasurementT());
        for (typename GroupT::iterator iter = group.begin(); iter != group.end(); ++iter) {
            meas.add(measureOne(*iter, source));
        }
        return meas;
    }
    
    /// Measure multiple values from groups of multiple images.
    /// Returns composite of MeasurementT (one measurement per group).
    /// Defaults to iteratively calling measureGroup and then averaging the results for each group.
    /// However, if the measurement requires linking a measurement across the groups (e.g., a fitted center
    /// from treating all the data), then the Algorithm needs to define this method.
    virtual PTR(MeasurementT) measureGroups(GroupSetT const& groups, afwDet::Source const& source) {
        PTR(MeasurementT) meas(new MeasurementT());
        for (typename GroupSetT::iterator iter = groups.begin(); iter != groups.end(); ++iter) {
            GroupT group = *iter;
            PTR(MeasurementT) measGroup = measureGroup(group, source);
            meas.add(measGroup->average());
        }
        return meas;
    }

    /// Return a null measurement
    ///
    /// This is called when we hit an exception.
    ///
    /// Pure-virtual, so subclasses MUST define.
    virtual PTR(MeasurementT) measureNull(void) = 0;
    
    /// Configure the algorithm
    virtual void configure(pexPolicy::Policy const&) {};
};


}}} // namespace

#endif
