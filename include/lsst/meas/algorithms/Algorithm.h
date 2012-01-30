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
#include "boost/make_shared.hpp"

#include "lsst/base.h"
#include "lsst/pex/logging/Log.h"
#include "lsst/pex/config.h"
#include "lsst/pex/policy.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/meas/algorithms/ExposurePatch.h"
#include "lsst/afw/table/Source.h"

namespace lsst { namespace meas { namespace algorithms {

class AlgorithmControl;

/// Base class for algorithms for measuring sources
template<typename ExposureT>
class Algorithm {
public:

    explicit Algorithm(AlgorithmControl const & ctrl);

    /// Destructor
    virtual ~Algorithm() {}

    std::string const & getName() const { return _name; }

    virtual void apply(
        afw::table::SourceRecord & source,
        ExposurePatch<ExposureT> const& exposure
    ) const = 0;

private:
    std::string _name;
};

#define LSST_ALGORITHM_CONTROL_PRIVATE_DECL_PIXEL(PIXEL)    \
    virtual PTR(Algorithm< afw::image::Exposure< PIXEL > >) \
    _makeAlgorithm(afw::image::Exposure< PIXEL > *, afw::table::Schema & schema) const
#define LSST_ALGORITHM_CONTROL_PRIVATE_DECL()           \
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL_PIXEL(float);   \
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL_PIXEL(double);
#define LSST_ALGORITHM_CONTROL_PRIVATE_IMPL_PIXEL(CTRL_CLS, ALG_CLS, PIXEL) \
    PTR(Algorithm< afw::image::Exposure< PIXEL > >) \
    CTRL_CLS::_makeAlgorithm(afw::image::Exposure< PIXEL > *, afw::table::Schema & schema) const { \
        return boost::make_shared< ALG_CLS< afw::image::Exposure< PIXEL > > >(*this, boost::ref(schema)); \
    }
#define LSST_ALGORITHM_CONTROL_PRIVATE_IMPL(CTRL_CLS, ALG_CLS)   \
    LSST_ALGORITHM_CONTROL_PRIVATE_IMPL_PIXEL(CTRL_CLS, ALG_CLS, float) \
    LSST_ALGORITHM_CONTROL_PRIVATE_IMPL_PIXEL(CTRL_CLS, ALG_CLS, double)

/**
 *  @brief Base class for measurement algorithm control objects.
 *
 *  Control objects are the public face of the measurement algorithms, and the only
 *  part available in Python.  Derived classes must be default-constructable,
 *  and will often be simple structs with public data members.  Data members in a control
 *  object will generally have a one-to-one correspondence with the fields in a Python
 *  Config object; it's unfortunate that these definitions are duplicated, and must
 *  be kept in sync, but it seems better than finding a way to generate both from a
 *  third definition.
 *
 *  This is a polymorphic class hierarchy because control objects are also factories
 *  for algorithms - but this is considered an implementation detail, and only matters
 *  to algorithm writers, who must implement the protected algorithm factory functions.
 *  The advantage of this approach is that we don't have to SWIG the algorithm classes
 *  at all, and we since the control classes aren't templated on the exposure type,
 *  the Python world doesn't need to know about that template parameter at all.  Instead,
 *  we can just control objects (by base class pointer) directly to MeasureSources and
 *  MeasureQuantity.
 */
class AlgorithmControl {
public:

    LSST_CONTROL_FIELD(name, std::string, "name used for the algorithm's measurements");

    virtual ~AlgorithmControl() {}
    
protected:

    explicit AlgorithmControl(std::string const & name_) : name(name_) {}

    LSST_ALGORITHM_CONTROL_PRIVATE_DECL_PIXEL(float) = 0;
    LSST_ALGORITHM_CONTROL_PRIVATE_DECL_PIXEL(double) = 0;
    
private:

    template <typename ExposureT> friend class MeasureSources;

    template <typename ExposureT>
    PTR(Algorithm<ExposureT>) makeAlgorithm(afw::table::Schema & schema) const {
        return _makeAlgorithm((ExposureT*)0, schema);
    }

};

template <typename ExposureT>
inline Algorithm<ExposureT>::Algorithm(AlgorithmControl const & ctrl) : _name(ctrl.name) {}

}}} // namespace lsst::meas::algorithms

#endif
