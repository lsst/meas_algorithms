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
#include "lsst/afw/detection/Source.h"
#include "lsst/meas/algorithms/ExposurePatch.h"

namespace lsst { namespace meas { namespace algorithms {

/// Base class for algorithms for measuring MeasurementT (e.g., Photometry)
template<typename MeasurementT, typename ExposureT>
class Algorithm {
public:
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
    virtual PTR(MeasurementT) measureSingle(afw::detection::Source const& target,
                                            afw::detection::Source const& source,
                                            ExposurePatch<ExposureT> const& patch) const = 0;
    
    /// Measure a single value from multiple images.
    ///
    /// Returns leaf of MeasurementT (single measurement).
    ///
    /// Because it is a 'group' of images (images in the same filter), we can
    /// assume they share some characteristics (e.g., center, shape).
    ///
    /// Defaults to iteratively calling measureOne. However, if the measurement
    /// cannot be obtained by merely averaging the outputs of a single
    /// measurement, e.g., some measured parameters are made across all
    /// exposures as part of the measurement (e.g., a common shape), then the
    /// Algorithm needs to define this method.
    virtual PTR(MeasurementT) measureMultiple(afw::detection::Source const& target,
                                              afw::detection::Source const& source,
                                              std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches
        ) const {
        typedef std::vector<CONST_PTR(ExposurePatch<ExposureT>)> PatchVector;
        PTR(MeasurementT) meas(new MeasurementT());
        for (typename PatchVector::const_iterator iter = patches.begin(); iter != patches.end(); ++iter) {
            PTR(MeasurementT) val;
            try {
                CONST_PTR(ExposurePatch<ExposureT>) patch = *iter;
                val = measureSingle(target, source, *patch);
            } catch (pex::exceptions::Exception const& e) {
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

    /// Return a null measurement
    ///
    /// This is called when we hit an exception.
    virtual PTR(MeasurementT) measureNull(void) const {
        return MeasurementT::null();
    }
    
    /// Configure the algorithm
    virtual void configure(pex::policy::Policy const&) {};

    /// Name of the algorithm
    virtual std::string getName() const = 0;

    /// Clone algorithm
    virtual PTR(Algorithm<MeasurementT, ExposureT>) clone() const = 0;
};


}}} // namespace

#endif
