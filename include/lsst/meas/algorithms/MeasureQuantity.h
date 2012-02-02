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
 
#if !defined(LSST_MEAS_ALGORITHMS_MEASURE_QUANTITY_H)
#define LSST_MEAS_ALGORITHMS_MEASURE_QUANTITY_H

//!
// Measure properties of an image selected by a Footprint
//

#include <cmath>
#include "boost/multi_index_container.hpp"
#include "boost/multi_index/indexed_by.hpp"
#include "boost/multi_index/sequenced_index.hpp"
#include "boost/multi_index/hashed_index.hpp"
#include "boost/multi_index/mem_fun.hpp"
#include "boost/preprocessor/cat.hpp"
#include "boost/preprocessor/stringize.hpp"
#include "boost/serialization/force_include.hpp"
#include "lsst/base.h"
#include "lsst/utils/Demangle.h"
#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/daf/base/Citizen.h"
#include "lsst/afw/detection/Schema.h"
#include "lsst/afw/detection/Peak.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/meas/algorithms/Algorithm.h"

namespace lsst { namespace meas { namespace algorithms {


namespace detail {

#ifdef SWIG
template<typename PtrAlgorithmT, typename AlgorithmT>
class AlgorithmMap;
#else
/// An insertion-ordered map of algorithms
///
/// A convenience class, wallpapering over the ugliness of boost::multi_index
template<typename PtrAlgorithmT, typename AlgorithmT>
class AlgorithmMap : public lsst::daf::base::Citizen {
private:
    typedef typename boost::multi_index_container<PtrAlgorithmT,
                                                  boost::multi_index::indexed_by<
                                                      boost::multi_index::sequenced<>,
                                                      boost::multi_index::hashed_unique<
                                                          boost::multi_index::const_mem_fun<
                                                              AlgorithmT, std::string,
                                                              &AlgorithmT::getName> >
                                                      > > AlgorithmMapT;
public:
    typedef typename AlgorithmMapT::template nth_index<0>::type::iterator iterator;
    typedef typename AlgorithmMapT::template nth_index<0>::type::const_iterator const_iterator;

    /// Constructor
    AlgorithmMap() : lsst::daf::base::Citizen(typeid(this)), _map() {}
    ~AlgorithmMap() {}

    /// Iterators
    iterator begin() { return _map.template get<0>().begin(); }
    const_iterator begin() const { return _map.template get<0>().begin(); }
    iterator end() { return _map.template get<0>().end(); }
    const_iterator end() const { return _map.template get<0>().end(); }

    std::size_t size() const { return _map.template get<0>().size(); }

    /// Append to the list of algorithms
    void append(PtrAlgorithmT alg) {
        _map.template get<0>().insert(_map.end(), alg);
    }

    /// Check for existence of an algorithm
    bool exists(std::string const& name) const {
        typename AlgorithmMapT::template nth_index<1>::type::const_iterator iter = 
            _map.template get<1>().find(name);
        return !(iter == _map.template get<1>().end());
    }

    /// Find an algorithm by name
    PtrAlgorithmT find(std::string const& name) const {
        typename AlgorithmMapT::template nth_index<1>::type::const_iterator iter = 
            _map.template get<1>().find(name);
        if (iter == _map.template get<1>().end()) {
            throw LSST_EXCEPT(pex::exceptions::NotFoundException, 
                              (boost::format("Algorithm %s not registered.") % name).str());
        }
        return *iter;
    }

private:
    AlgorithmMapT _map;                 // Insertion-ordered map
};
#endif

/// Templated typedefs for insertion-ordered maps of pointers to algorithms
template<typename AlgorithmT>
struct AlgorithmMapTypes {
    typedef typename AlgorithmMap<PTR(AlgorithmT), AlgorithmT>::AlgorithmMap PtrAlgorithmMap;
    typedef typename AlgorithmMap<CONST_PTR(AlgorithmT), AlgorithmT>::AlgorithmMap ConstPtrAlgorithmMap;
};

/// Measuring sources on the same image on which they were detected.
template<typename ExposureT>
struct SingleMeasurer {
    typedef ExposurePatch<ExposureT> ExposureContainerT;
    typedef ExposurePatch<ExposureT> const ConstExposureContainerT;

    /// Execute the algorithm
    template<typename MeasurementT>
    static PTR(MeasurementT) algorithm(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                       afw::detection::Source const& target,
                                       afw::detection::Source const& source,
                                       ExposurePatch<ExposureT> const& patch) {
        return alg->measureSingle(target, source, patch);
    }                                       
};

/// Measuring a single source on multiple images
template<typename ExposureT>
struct MultipleMeasurer {
    typedef std::vector<PTR(ExposurePatch<ExposureT>)> ExposureContainerT;
    typedef std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const ConstExposureContainerT;
    
    /// Execute the algorithm
    template<typename MeasurementT>
    static PTR(MeasurementT) algorithm(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                       afw::detection::Source const& target,
                                       afw::detection::Source const& source,
                                       std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches) {
        return alg->measureMultiple(target, source, patches);
    }                                       
};

} // namespace detail



/************************************************************************************************************/
/*
 * Measure a quantity using a set of algorithms.  Each algorithm will fill one item in the returned
 * measurement
 */
template<typename MeasurementT, typename ExposureT>
class MeasureQuantity {
public:
    typedef PTR(MeasureQuantity) Ptr;
    typedef CONST_PTR(MeasureQuantity) ConstPtr;
    typedef Algorithm<MeasurementT, ExposureT> AlgorithmT;
    typedef ExposurePatch<ExposureT> ExposurePatchT;
private:
    typedef typename detail::AlgorithmMapTypes<AlgorithmT>::PtrAlgorithmMap PtrAlgorithmMapT;
    typedef typename detail::AlgorithmMapTypes<AlgorithmT>::ConstPtrAlgorithmMap ConstPtrAlgorithmMapT;

public:

    /// Constructor
    explicit MeasureQuantity() : _algorithms() {}

    virtual ~MeasureQuantity() {}

    /// Include the algorithm defined by the given control object in the list of measurement algorithms.
    void addAlgorithm(AlgorithmControl<MeasurementT> const & ctrl) {
        _algorithms.append(ctrl.template makeAlgorithm<ExposureT>());
    }

    /// Methods to make the measurements.

    /// Measure a single exposure
    PTR(MeasurementT) measure(
        afw::detection::Source const& target,   ///< Source being measured
        CONST_PTR(ExposureT) exp,       ///< Exposure to measure
        afw::geom::Point2D const& center, ///< Center of source
        pex::logging::Log &log=pex::logging::Log::getDefaultLog() ///< Log
        ) const {
        ExposurePatch<ExposureT> patch(exp, target.getFootprint(), center);
        return measureSingle(target, target, patch, log);
    }

    /// Measure multiple exposures
    PTR(MeasurementT) measure(
        afw::detection::Source const& target,   ///< Source being measured
        afw::detection::Source const& source,   ///< Original source
        afw::image::Wcs const& wcs,       ///< WCS for original source
        std::vector<typename ExposureT::ConstPtr> const& exposures, ///< Exposures being measured
        pex::logging::Log &log=pex::logging::Log::getDefaultLog() ///< Log
        ) const {
        size_t size = exposures.size();
        std::vector<CONST_PTR(ExposurePatch<ExposureT>)> patches(size);
        afw::geom::Point2D center(source.getXAstrom(), source.getYAstrom());
        for (size_t i = 0; i != size; ++i) {
            patches[i] = makeExposurePatch(exposures[i], *source.getFootprint(), center, wcs);
        }
        return measureMultiple(target, source, patches, log);
    }

    /// Measure a single exposure patch
    PTR(MeasurementT) measureSingle(
        afw::detection::Source const& target,   ///< Source being measured
        afw::detection::Source const& source,   ///< Original source
        ExposurePatchT const& patch,    ///< Patch being measured
        pex::logging::Log &log=pex::logging::Log::getDefaultLog() ///< Log
        ) const {
        return _measure<detail::SingleMeasurer<ExposureT> >(target, source, patch, log);
    }
    
    /// Measure multiple exposure patches
    PTR(MeasurementT) measureMultiple(
        afw::detection::Source const& target,   ///< Source being measured
        afw::detection::Source const& source,   ///< Original source
        std::vector<CONST_PTR(ExposurePatch<ExposureT>)> const& patches, ///< Patches being measured
        pex::logging::Log &log=pex::logging::Log::getDefaultLog()            ///< Log
        ) const {
        return _measure<detail::MultipleMeasurer<ExposureT> >(target, source, patches, log);
    }

    /// List activated algorithms
    std::vector<std::string> listActive() const {
        std::vector<std::string> names = std::vector<std::string>();
        for (typename PtrAlgorithmMapT::const_iterator iter = _algorithms.begin(); 
             iter != _algorithms.end(); ++iter) {
            CONST_PTR(AlgorithmT) alg = *iter;
            names.push_back(alg->getName());
        }
        return names;
    }

private:
    /// Measure the appropriate exposure container with each algorithm
    ///
    /// Measurer is a class that does the appropriate measurement for
    /// ExposureContainerT as a static method called algorithm().
    template<class Measurer>
    PTR(MeasurementT) _measure(
        afw::detection::Source const& target, ///< Source being measured
        afw::detection::Source const& source, ///< Original source
        typename Measurer::ConstExposureContainerT const& exp, ///< Exposure patch to measure
        pex::logging::Log &log                           ///< Log
        ) const {
        PTR(MeasurementT) values = boost::make_shared<MeasurementT>();        
        for (typename PtrAlgorithmMapT::const_iterator iter = _algorithms.begin();
             iter != _algorithms.end(); ++iter) {
            CONST_PTR(AlgorithmT) alg = *iter; // Algorithm to execute
            std::string const name = alg->getName(); // Name of algorithm
            PTR(MeasurementT) val;                    // Value measured by algorithm
            try {
#if 0
                std::cerr << (boost::format("Measuring %s at %f,%f") %
                              name % source.getXAstrom() % source.getYAstrom()).str() << std::endl;
#endif
                val = Measurer::algorithm(alg, target, source, exp);
            } catch (pex::exceptions::Exception const& e) {
                // Swallow all exceptions, because one bad measurement shouldn't affect all others
                log.log(pex::logging::Log::DEBUG, boost::format("Measuring %s at (%d,%d): %s") %
                        name % source.getXAstrom() % source.getYAstrom() % e.what());
                val = alg->measureNull();
            }
        
            val->setAlgorithm(name);
            values->add(val);
        }
    
        return values;
    }

    PtrAlgorithmMapT _algorithms;     /// The list of algorithms that we wish to use

};


/// Specialisation of MeasureQuantity for Astrometry measurements
template<typename ExposureT>
class MeasureAstrometry : public MeasureQuantity<afw::detection::Astrometry, ExposureT> {
public:
    typedef PTR(MeasureAstrometry) Ptr;

    MeasureAstrometry() : MeasureQuantity<afw::detection::Astrometry, ExposureT>() {}
};

/// Specialisation of MeasureQuantity for Shape measurements
template<typename ExposureT>
class MeasureShape : public MeasureQuantity<afw::detection::Shape, ExposureT> {
public:
    typedef PTR(MeasureShape) Ptr;

    MeasureShape() : MeasureQuantity<afw::detection::Shape, ExposureT>() {}
};

/// Specialisation of MeasureQuantity for Photometry measurements
template<typename ExposureT>
class MeasurePhotometry : public MeasureQuantity<afw::detection::Photometry, ExposureT> {
public:
    typedef PTR(MeasurePhotometry) Ptr;

    MeasurePhotometry() : MeasureQuantity<afw::detection::Photometry, ExposureT>() {}
};

/// Factory functions for MeasureQuantity specialisations
template<typename ExposureT>
PTR(MeasureAstrometry<ExposureT>) makeMeasureAstrometry(ExposureT const& exp) {
    return boost::make_shared<MeasureAstrometry<ExposureT> >();
}
template<typename ExposureT>
PTR(MeasurePhotometry<ExposureT>) makeMeasurePhotometry(ExposureT const& exp) {
    return boost::make_shared<MeasurePhotometry<ExposureT> >();
}
template<typename ExposureT>
PTR(MeasureShape<ExposureT>) makeMeasureShape(ExposureT const& exp) {
    return boost::make_shared<MeasureShape<ExposureT> >();
}

}}} // namespace

#endif
