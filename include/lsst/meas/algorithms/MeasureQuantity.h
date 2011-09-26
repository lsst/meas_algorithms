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
#include "lsst/base.h"
#include "lsst/utils/Demangle.h"
#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/afw/detection/Schema.h"
#include "lsst/afw/detection/Peak.h"
#include "lsst/afw/detection/Measurement.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/meas/algorithms/Algorithm.h"


namespace lsst { namespace meas { namespace algorithms {

namespace pexLogging = lsst::pex::logging;
namespace pexPolicy = lsst::pex::policy;
namespace afwDet = lsst::afw::detection;

namespace {

/// An insertion-ordered map of algorithms
///
/// A convenience class, wallpapering over the ugliness of boost::multi_index
template<typename PtrAlgorithmT, typename AlgorithmT>
class AlgorithmMap {
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
    AlgorithmMap() : _map() {}
    ~AlgorithmMap() {}

    /// Iterators
    iterator begin() { return _map.template get<0>().begin(); }
    const_iterator begin() const { return _map.template get<0>().begin(); }
    iterator end() { return _map.template get<0>().end(); }
    const_iterator end() const { return _map.template get<0>().end(); }

    /// Append to the list of algorithms
    void append(PtrAlgorithmT alg) {
        _map.template get<0>().insert(_map.end(), alg);
    }

    /// Find an algorithm by name
    PtrAlgorithmT find(std::string const& name) const {
        typename AlgorithmMapT::template nth_index<1>::type::const_iterator iter = 
            _map.template get<1>().find(name);
        if (iter == _map.template get<1>().end()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException, 
                              (boost::format("Algorithm %s not registered.") % name).str());
        }
        return *iter;
    }

private:
    AlgorithmMapT _map;                 // Insertion-ordered map
};


/// Templated typedefs for insertion-ordered maps of pointers to algorithms
template<typename AlgorithmT>
struct AlgorithmMapTypes {
    typedef typename AlgorithmMap<PTR(AlgorithmT), AlgorithmT>::AlgorithmMap PtrAlgorithmMap;
    typedef typename AlgorithmMap<CONST_PTR(AlgorithmT), AlgorithmT>::AlgorithmMap ConstPtrAlgorithmMap;
};


/// Singleton algorithm registry
template<typename AlgorithmT>
typename AlgorithmMapTypes<AlgorithmT>::ConstPtrAlgorithmMap& _getRegisteredAlgorithms() {
    static typename AlgorithmMapTypes<AlgorithmT>::ConstPtrAlgorithmMap *registeredAlgorithms = 
        new typename AlgorithmMapTypes<AlgorithmT>::ConstPtrAlgorithmMap();
    return *registeredAlgorithms;
}

/// Measurers for use with MeasureQuantity::_measure
///
/// These make up for the lack of support for passing method names: we simply
/// template on the Measurer class, and call a class (static) method of a
/// standard name.  The optimised compiled code should end up calling the
/// appropriate method name directly.
template<typename MeasurementT, typename ExposureT>
struct SingleAlgMeasurer {
    typedef ExposurePatch<ExposureT> ExposureContainerT;
    static PTR(MeasurementT) measure(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                     ExposureContainerT const& exp,
                                     afwDet::Source const& source) {
        return alg->measureOne(exp, source);
    }
};
template<typename MeasurementT, typename ExposureT>
struct GroupAlgMeasurer {
    typedef ExposureGroup<ExposureT> ExposureContainerT;
    static PTR(MeasurementT) measure(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                     ExposureContainerT const& group,
                                     afwDet::Source const& source) {
        return alg->measureGroup(group, source);
    }
};
template<typename MeasurementT, typename ExposureT>
struct GroupsAlgMeasurer {
    typedef typename std::vector<ExposureGroup<ExposureT> > ExposureContainerT;
    static PTR(MeasurementT) measure(CONST_PTR(Algorithm<MeasurementT, ExposureT>) alg,
                                     ExposureContainerT const& groups,
                                     afwDet::Source const& source) {
        return alg->measureGroups(groups, source);
    }
};


} // Anonymous namespace



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
    typedef ExposureGroup<ExposureT> ExposureGroupT;
private:
    typedef typename AlgorithmMapTypes<AlgorithmT>::PtrAlgorithmMap PtrAlgorithmMapT;
    typedef typename AlgorithmMapTypes<AlgorithmT>::ConstPtrAlgorithmMap ConstPtrAlgorithmMapT;

public:

    MeasureQuantity(CONST_PTR(pexPolicy::Policy) policy=CONST_PTR(pexPolicy::Policy)()) : 
        _algorithms() {
        configure(policy);
    }
    virtual ~MeasureQuantity() {}

    /// Class methods to instantiate based on type of image container
    ///
    /// These do nothing with the exposure (container) except use it for type determination
    MeasureQuantity<MeasurementT, ExposureT> create(
        ExposureT const& exp, 
        CONST_PTR(pexPolicy::Policy) policy=CONST_PTR(pexPolicy::Policy)()) {
        return MeasureQuantity<MeasurementT, ExposureT>(policy);
    }
    MeasureQuantity<MeasurementT, ExposureT> create(
        ExposureGroupT const& group, 
        CONST_PTR(pexPolicy::Policy) policy=CONST_PTR(pexPolicy::Policy)()) {
        return MeasureQuantity<MeasurementT, ExposureT>(policy);
    }
    MeasureQuantity<MeasurementT, ExposureT> create(
        std::vector<ExposureGroupT> const& groups, 
        CONST_PTR(pexPolicy::Policy) policy=CONST_PTR(pexPolicy::Policy)()) {
        return MeasureQuantity<MeasurementT, ExposureT>(policy);
    }

    /// Include the algorithm called name in the list of measurement algorithms to use
    ///
    /// This name is looked up in the registry (\sa declare), and used as the name of the
    /// measurement if you wish to retrieve it using the schema
    ///
    PTR(AlgorithmT) addAlgorithm(std::string const& name ///< The name of the algorithm
                     ) {
        CONST_PTR(AlgorithmT) regAlg = _lookupRegisteredAlgorithm(name); // Registered algorithm
        PTR(AlgorithmT) alg = regAlg->clone();                           // Algorithm to use
        _algorithms.append(alg);
        return alg;
    }

    PTR(MeasurementT) measure(ExposurePatchT const& patch,
                              afwDet::Source const& source,
                              pexLogging::Log &log=pexLogging::Log::getDefaultLog()) const {
        return _measure<SingleAlgMeasurer<MeasurementT, ExposureT> >(patch, source, log);
    }
    
    PTR(MeasurementT) measureGroup(ExposureGroupT const& group,
                                   afwDet::Source const& source,
                                   pexLogging::Log &log=pexLogging::Log::getDefaultLog()) const {
        return _measure<GroupAlgMeasurer<MeasurementT, ExposureT> >(group, source, log);
    }
    
    PTR(MeasurementT) measureGroups(std::vector<ExposureGroupT> const& groups,
                                    afwDet::Source const& source,
                                    pexLogging::Log &log=pexLogging::Log::getDefaultLog()) const {
        return _measure<GroupsAlgMeasurer<MeasurementT, ExposureT> >(groups, source, log);
    }
    
    /// Configure active algorithms and their parameters
    void configure(CONST_PTR(pexPolicy::Policy) policy ///< The Policy to configure algorithms
        ) {
        pexPolicy::Policy::StringArray names = policy->policyNames(false);

        for (pexPolicy::Policy::StringArray::iterator iter = names.begin();
             iter != names.end(); ++iter) {
            std::string const name = *iter;
            pexPolicy::Policy::ConstPtr subPol = policy->getPolicy(name);
            if (!subPol->exists("enabled") || subPol->getBool("enabled")) {
                PTR(AlgorithmT) alg = addAlgorithm(name);
                alg->configure(*subPol);
            }
        }
    }

    /// Declare an algorithm's existence
    static bool declare(CONST_PTR(AlgorithmT) alg) {
        ConstPtrAlgorithmMapT &registered = _getRegisteredAlgorithms<AlgorithmT>();
        registered.append(alg);
        return true;
    }

private:
    typename ExposureT::ConstPtr _im;   // The data that we wish to measure
    PtrAlgorithmMapT _algorithms;     // The list of algorithms that we wish to use

    /// Lookup a registered algorithm
    static CONST_PTR(AlgorithmT) _lookupRegisteredAlgorithm(std::string const& name) {
        ConstPtrAlgorithmMapT const& registered = _getRegisteredAlgorithms<AlgorithmT>();
        return registered.find(name);
    }

public:
    /// Measure the appropriate exposure container with each algorithm
    ///
    /// Measurer is a class that does the appropriate measurement for
    /// ExposureContainerT as a static method called measure().
    template<class Measurer>
    PTR(MeasurementT) _measure(typename Measurer::ExposureContainerT const& exp, afwDet::Source const& source,
                               pexLogging::Log &log) const {
        PTR(MeasurementT) values = boost::make_shared<MeasurementT>();
        
        for (typename PtrAlgorithmMapT::const_iterator iter = _algorithms.begin(); 
             iter != _algorithms.end(); ++iter) {
            CONST_PTR(AlgorithmT) alg = *iter; // Algorithm to execute
            std::string const name = alg->getName(); // Name of algorithm
            PTR(MeasurementT) val;                    // Value measured by algorithm
            try {
                val = Measurer::measure(alg, exp, source);
            } catch  (lsst::pex::exceptions::Exception const& e) {
                // Swallow all exceptions, because one bad measurement shouldn't affect all others
                log.log(pexLogging::Log::DEBUG, boost::format("Measuring %s at (%d,%d): %s") %
                        name % source.getXAstrom() % source.getYAstrom() % e.what());
                val = alg->measureNull();
            }
            val->getSchema()->setComponent(name); // name this type of measurement (e.g. psf)
            values->add(val);
        }

        return values;
    }

};

/// Declare the availability of an algorithm
///
/// Instantiates and registers the algorithm
#define DECLARE_ALGORITHM(MEASUREMENT, ALGORITHM, PIXEL) \
namespace { \
    static volatile PTR(ALGORITHM) instance(new ALGORITHM); \
    MeasureQuantity<MEASUREMENT, lsst::afw::image::Exposure<PIXEL> >::declare(instance); \
}


/// Specialisation of MeasureQuantity for Astrometry measurements
template<typename ExposureT>
class MeasureAstrometry : public MeasureQuantity<lsst::afw::detection::Astrometry, ExposureT> {
public:
    typedef PTR(MeasureAstrometry) Ptr;

    MeasureAstrometry(ExposureT const& exp,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Astrometry, ExposureT>(policy) {}
    MeasureAstrometry(ExposureGroup<ExposureT> const& group,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Astrometry, ExposureT>(policy) {}
    MeasureAstrometry(std::vector<ExposureGroup<ExposureT> > const& groups,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Astrometry, ExposureT>(policy) {}
};

/// Specialisation of MeasureQuantity for Shape measurements
template<typename ExposureT>
class MeasureShape : public MeasureQuantity<lsst::afw::detection::Shape, ExposureT> {
public:
    typedef PTR(MeasureShape) Ptr;

    MeasureShape(ExposureT const& exp,
                 CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Shape, ExposureT>(policy) {}
    MeasureShape(ExposureGroup<ExposureT> const& group,
                 CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Shape, ExposureT>(policy) {}
    MeasureShape(std::vector<ExposureGroup<ExposureT> > const& groups,
                 CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Shape, ExposureT>(policy) {}
};

/// Specialisation of MeasureQuantity for Photometry measurements
template<typename ExposureT>
class MeasurePhotometry : public MeasureQuantity<lsst::afw::detection::Photometry, ExposureT> {
public:
    typedef PTR(MeasurePhotometry) Ptr;
    MeasurePhotometry(ExposureT const& exp,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Photometry, ExposureT>(exp, policy) {}
    MeasurePhotometry(ExposureGroup<ExposureT> const& group,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Photometry, ExposureT>(group, policy) {}
    MeasurePhotometry(std::vector<ExposureGroup<ExposureT> > const& groups,
                      CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)()) : 
        MeasureQuantity<lsst::afw::detection::Photometry, ExposureT>(groups, policy) {}
};

}}} // namespace

#endif
