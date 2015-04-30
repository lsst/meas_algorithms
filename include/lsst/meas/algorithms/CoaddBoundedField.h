// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#ifndef LSST_MEAS_ALGORITHMS_CoaddBoundedField_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_CoaddBoundedField_h_INCLUDED

#include "lsst/afw/math/BoundedField.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/polygon/Polygon.h"

namespace lsst { namespace meas { namespace algorithms {

/// Struct used to hold one Exposure's data in a CoaddBoundedField
struct CoaddBoundedFieldElement {

    CoaddBoundedFieldElement(
        PTR(afw::math::BoundedField) field_,
        PTR(afw::image::Wcs const) wcs_,
        PTR(afw::geom::polygon::Polygon const) validPolygon_,
        double weight_=1.0
        ) : field(field_), wcs(wcs_), validPolygon(validPolygon_),
            weight(weight_) {}

    PTR(afw::math::BoundedField) field;
    PTR(afw::image::Wcs const) wcs;
    PTR(afw::geom::polygon::Polygon const) validPolygon;
    double weight;
};

class CoaddBoundedField :
    public afw::table::io::PersistableFacade<CoaddBoundedField>,
    public afw::math::BoundedField
{
public:

    typedef CoaddBoundedFieldElement Element;
    typedef std::vector<Element> ElementVector;

    explicit CoaddBoundedField(
        afw::geom::Box2I const & bbox,
        PTR(afw::image::Wcs const) coaddWcs,
        ElementVector const & elements
    );

    explicit CoaddBoundedField(
        afw::geom::Box2I const & bbox,
        PTR(afw::image::Wcs const) coaddWcs,
        ElementVector const & elements,
        double default_
    );

    /// @copydoc afw::math::BoundedField::evaluate
    virtual double evaluate(afw::geom::Point2D const & position) const;

    /**
     *  @brief Return true if the CoaddBoundedField persistable (always true).
     */
    virtual bool isPersistable() const { return true; }

    // Factory used to read CoaddPsf from an InputArchive; defined only in the source file.
    class Factory;

    virtual PTR(afw::math::BoundedField) operator*(double const scale) const;

protected:

    // See afw::table::io::Persistable::getPersistenceName
    virtual std::string getPersistenceName() const;

    // See afw::table::io::Persistable::getPythonModule
    virtual std::string getPythonModule() const;

    // See afw::table::io::Persistable::write
    virtual void write(OutputArchiveHandle & handle) const;

private:

    bool _throwOnMissing;     // instead of using _default, raise an exception
    double _default;          // when none of the elements contribute at a point, return this value
    PTR(afw::image::Wcs const) _coaddWcs;  // coordinate system this field is defined in
    ElementVector _elements;  // vector of constituent fields being coadded
};

}}} // namespace lsst::meas::algorithms

#endif // !LSST_MEAS_ALGORITHMS_CoaddBoundedField_h_INCLUDED
