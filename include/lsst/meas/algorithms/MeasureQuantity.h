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
                                                              AlgorithmT, std::string&,
                                                              &AlgorithmT::getName> >
                                                      > > AlgorithmMapT;
public:
    typedef typename AlgorithmMapT::template nth_index<0>::type::iterator iterator;
    typedef typename AlgorithmMapT::template nth_index<0>::type::const_iterator const_iterator;

    /// Constructor
    AlgorithmMap() : _map() {}
    ~AlgorithmMap() {}

    /// Iterators
    iterator begin() { return _map.template get<0>.begin(); }
    const_iterator begin() const { return _map.template get<0>.begin(); }
    iterator end() { return _map.template get<0>.end(); }
    const_iterator end() const { return _map.template get<0>.end(); }

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

} // Anonymous namespace



/************************************************************************************************************/
/*
 * Measure a quantity using a set of algorithms.  Each algorithm will fill one item in the returned
 * measurement
 */
template<typename MeasurementT, typename ExposureT>
class MeasureQuantity {
public:
    typedef Algorithm<MeasurementT, ExposureT> AlgorithmT;
    typedef ExposurePatch<ExposureT> ExposurePatchT;
private:
    typedef typename AlgorithmMapTypes<AlgorithmT>::PtrAlgorithmMap PtrAlgorithmMapT;
    typedef typename AlgorithmMapTypes<AlgorithmT>::ConstPtrAlgorithmMap ConstPtrAlgorithmMapT;

public:

    MeasureQuantity(typename ExposureT::ConstPtr im,
                    CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)())
        : _im(im), _algorithms()
    {
        if (policy) {
            configure(*policy);
        }
    }
    virtual ~MeasureQuantity() {}

    /**
     * Return the image data that we are measuring
     */
    typename ExposureT::ConstPtr getImage() const {
        return _im;
    }
    /**
     * (Re)set the data that we are measuring
     */
    void setImage(typename ExposureT::ConstPtr im) {
        _im = im;
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

    /// Actually measure im using all requested algorithms, returning the result
    PTR(MeasurementT) measureOne(ExposurePatchT const& patch, afwDet::Source const& source,
                                 pexLogging::Log &log=pexLogging::Log::getDefaultLog()
        ) {
        PTR(MeasurementT) values = boost::make_shared<MeasurementT>();
        
        if (!_im) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              "I cannot measure a NULL image");
        }

        for (typename PtrAlgorithmMapT::const_iterator iter = _algorithms.begin(); 
             iter != _algorithms.end(); ++iter) {
            CONST_PTR(AlgorithmT) alg = *iter; // Algorithm to execute
            std::string const& name = alg->getName(); // Name of algorithm
            PTR(MeasurementT) val;                    // Value measured by algorithm
            try {
                val = alg->measureOne(patch, source);
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

    /// Configure active algorithms and their parameters
    void configure(lsst::pex::policy::Policy const& policy ///< The Policy to configure algorithms
        ) {
        lsst::pex::policy::Policy::StringArray names = policy.policyNames(false);

        for (lsst::pex::policy::Policy::StringArray::iterator iter = names.begin();
             iter != names.end(); ++iter) {
            std::string const name = *iter;
            lsst::pex::policy::Policy::ConstPtr subPol = policy.getPolicy(name);
            if (!subPol->exists("enabled") || subPol->getBool("enabled")) {
                PTR(AlgorithmT) alg = addAlgorithm(name);
                alg->configure(*subPol);
            }
        }
    }

    /// Declare an algorithm's existence
    static void declare(CONST_PTR(AlgorithmT) alg) {
        ConstPtrAlgorithmMapT &registered = _getRegisteredAlgorithms<AlgorithmT>();
        registered.append(alg);
    }

private:
    typename ExposureT::ConstPtr _im;   // The data that we wish to measure
    PtrAlgorithmMapT _algorithms;     // The list of algorithms that we wish to use

    /// Lookup a registered algorithm
    static CONST_PTR(AlgorithmT) _lookupRegisteredAlgorithm(std::string const& name) {
        ConstPtrAlgorithmMapT const& registered = _getRegisteredAlgorithms<AlgorithmT>();
        return registered.find(name);
    }
};


#define DECLARE_ALGORITHM(MEASUREMENT, ALGORITHM, PIXEL) \
namespace { \
    typedef lsst::afw::image::Exposure<PIXEL> ExposureT; \
    static volatile PTR(ALGORITHM<MEASUREMENT, ExposureT>) instance(new ALGORITHM<MEASUREMENT, ExposureT>); \
    MeasureQuantity<MEASUREMENT, ExposureT>::declare(instance); \
}


}}} // namespace

#endif
