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

namespace lsst { namespace meas { namespace algorithms {

namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;

/************************************************************************************************************/
/*
 * Measure a quantity using a set of algorithms.  Each algorithm will fill one item in the returned
 * Values (a Measurement)
 */
template<typename T, typename ImageT>
class MeasureQuantity {
public:
    typedef afwDet::Measurement<T> Values;
    typedef boost::shared_ptr<T> (*makeMeasureQuantityFunc)(typename ImageT::ConstPtr,
                                                            CONST_PTR(afwDet::Peak), 
                                                            CONST_PTR(afwDet::Source));
    typedef bool (*configureMeasureQuantityFunc)(lsst::pex::policy::Policy const&);
    typedef std::pair<makeMeasureQuantityFunc, configureMeasureQuantityFunc> measureQuantityFuncs;
private:
    typedef std::map<std::string, measureQuantityFuncs> AlgorithmList;
public:

    MeasureQuantity(typename ImageT::ConstPtr im,
                    CONST_PTR(lsst::pex::policy::Policy) policy=CONST_PTR(lsst::pex::policy::Policy)())
        : _im(im), _algorithms()
    {
        if (policy) {
            lsst::pex::policy::Policy::StringArray names = policy->policyNames(false);

            for (lsst::pex::policy::Policy::StringArray::iterator ptr = names.begin();
                 ptr != names.end(); ++ptr) {
                lsst::pex::policy::Policy::ConstPtr subPol = policy->getPolicy(*ptr);

                if (!subPol->exists("enabled") || subPol->getBool("enabled")) {
                    addAlgorithm(*ptr);
                }
            }

            configure(*policy);
        }
    }
    virtual ~MeasureQuantity() {}

    /**
     * Return the image data that we are measuring
     */
    typename ImageT::ConstPtr getImage() const {
        return _im;
    }
    /**
     * (Re)set the data that we are measuring
     */
    void setImage(typename ImageT::ConstPtr im) {
        _im = im;
    }

    /// Include the algorithm called name in the list of measurement algorithms to use
    ///
    /// This name is looked up in the registry (\sa declare), and used as the name of the
    /// measurement if you wish to retrieve it using the schema
    ///
    void addAlgorithm(std::string const& name ///< The name of the algorithm
                     ) {
        _algorithms[name] = _lookupAlgorithm(name);
    }
    /// Actually measure im using all requested algorithms, returning the result
    PTR(Values) measure(CONST_PTR(afwDet::Peak) peak=PTR(afwDet::Peak)(),
                        CONST_PTR(afwDet::Source) src=PTR(afwDet::Source)(),
                        pexLogging::Log &log=pexLogging::Log::getDefaultLog()
                       ) {
        PTR(Values) values = boost::make_shared<Values>();

        if (!_im) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              "I cannot measure a NULL image");
        }

        for (typename AlgorithmList::iterator ptr = _algorithms.begin(); ptr != _algorithms.end(); ++ptr) {
            boost::shared_ptr<T> val;
            try {
                val = ptr->second.first(_im, peak, src);
            } catch (lsst::pex::exceptions::Exception const& e) {
                // Swallow all exceptions, because one bad measurement shouldn't affect all others
                log.log(pexLogging::Log::DEBUG, boost::format("Measuring %s at (%d,%d): %s") %
                        ptr->first % peak->getIx() % peak->getIy() % e.what());
                // Blank measure should set blank values
                val = ptr->second.first(_im, boost::shared_ptr<afwDet::Peak>(), 
                                        boost::shared_ptr<afwDet::Source>());
            }
            val->getSchema()->setComponent(ptr->first); // name this type of measurement (e.g. psf)
            values->add(val);
        }

        return values;
    }
    PTR(Values) measure(pexLogging::Log &log) {
        return measure(PTR(afwDet::Peak)(), PTR(afwDet::Source)(), log);
    }

    /// Configure the behaviour of the algorithm
    bool configure(lsst::pex::policy::Policy const& policy ///< The Policy to configure algorithms
                  ) {
        bool value = true;

        for (typename AlgorithmList::iterator ptr = _algorithms.begin(); ptr != _algorithms.end(); ++ptr) {
            if (policy.exists(ptr->first)) {
                lsst::pex::policy::Policy::ConstPtr subPol = policy.getPolicy(ptr->first);
                if (!subPol->exists("enabled") || subPol->getBool("enabled")) {
                    value = ptr->second.second(*subPol) && value; // don't short-circuit the call
                }
            }
        }

        return value;
    }

    static bool declare(std::string const& name,
        typename MeasureQuantity<T, ImageT>::makeMeasureQuantityFunc makeFunc,
        typename MeasureQuantity<T, ImageT>::configureMeasureQuantityFunc configFunc=_iefbr15);
private:
    //
    // The data that we wish to measure
    //
    typename ImageT::ConstPtr _im;
    //
    // The list of algorithms that we wish to use
    //
    AlgorithmList _algorithms;
    //
    // A mapping from names to algorithms
    //
    // _registryWorker must be inline as it contains a critical static variable, _registry
    //    
    typedef std::map<std::string, measureQuantityFuncs> AlgorithmRegistry;

    static inline measureQuantityFuncs _registryWorker(std::string const& name,
        typename MeasureQuantity<T, ImageT>::makeMeasureQuantityFunc makeFunc,
        typename MeasureQuantity<T, ImageT>::configureMeasureQuantityFunc configFunc
                                                      );
    static measureQuantityFuncs _lookupAlgorithm(std::string const& name);
    /// The unknown algorithm; used to allow _lookupAlgorithm use _registryWorker
    static boost::shared_ptr<T> _iefbr14(typename ImageT::ConstPtr, CONST_PTR(afwDet::Peak), 
                                         CONST_PTR(afwDet::Source)) {
        return boost::shared_ptr<T>();
    }
public:                                 // needed for swig to support keyword arguments
    static bool _iefbr15(lsst::pex::policy::Policy const &) {
        return true;
    }
private:
    //
    // Do the real work of measuring things
    //
    // Can't be pure virtual as we create a do-nothing MeasureQuantity which we then add to
    //
    virtual boost::shared_ptr<T> doMeasure(CONST_PTR(ImageT),
                                           CONST_PTR(afwDet::Peak),
                                           CONST_PTR(afwDet::Source) src=PTR(afwDet::Source)()
                                          ) {
        return boost::shared_ptr<T>();
    }
};

/**
 * Support the algorithm registry
 */
template<typename T, typename ImageT>
typename MeasureQuantity<T, ImageT>::measureQuantityFuncs
MeasureQuantity<T, ImageT>::_registryWorker(
        std::string const& name,
        typename MeasureQuantity<T, ImageT>::makeMeasureQuantityFunc makeFunc,
        typename MeasureQuantity<T, ImageT>::configureMeasureQuantityFunc configureFunc
                                                  )
{
    // N.b. This is a pointer rather than an object as this helps the intel compiler generate a
    // single copy of the _registry across multiple dynamically loaded libraries.  The intel
    // bug ID for RHL's report is 580524
    static typename MeasureQuantity<T, ImageT>::AlgorithmRegistry *_registry = NULL;

    if (!_registry) {
        _registry = new typename MeasureQuantity<T, ImageT>::AlgorithmRegistry;
    }

    if (makeFunc == _iefbr14) {     // lookup functions
        typename MeasureQuantity<T, ImageT>::AlgorithmRegistry::const_iterator ptr =
            _registry->find(name);
        
        if (ptr == _registry->end()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                              (boost::format("Unknown algorithm %s for image of type %s")
                               % name % lsst::utils::demangleType(typeid(ImageT).name())).str());
        }
        
        return ptr->second;
    } else {                            // register functions
        typename MeasureQuantity<T, ImageT>::measureQuantityFuncs funcs = 
            std::make_pair(makeFunc, configureFunc);            

        (*_registry)[name] = funcs;

        return funcs;
    }
}

/**
 * Register the factory function for a named algorithm
 */
template<typename T, typename ImageT>
bool MeasureQuantity<T, ImageT>::declare(
        std::string const& name,
        typename MeasureQuantity<T, ImageT>::makeMeasureQuantityFunc makeFunc,
        typename MeasureQuantity<T, ImageT>::configureMeasureQuantityFunc configFunc
                                               )
{
    _registryWorker(name, makeFunc, configFunc);

    return true;
}

/**
 * Return the factory function for a named algorithm
 */
template<typename T, typename ImageT>
typename MeasureQuantity<T, ImageT>::measureQuantityFuncs
MeasureQuantity<T, ImageT>::_lookupAlgorithm(std::string const& name)
{
    return _registryWorker(name, _iefbr14, _iefbr15);
}


}}} // namespace

#endif
