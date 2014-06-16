// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014, 2010 LSST Corporation.
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

#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/meas/algorithms/CoaddBoundedField.h"

namespace lsst { namespace meas { namespace algorithms {


CoaddBoundedField::CoaddBoundedField(
    afw::geom::Box2I const & bbox,
    PTR(afw::image::Wcs const) coaddWcs,
    ElementVector const & elements
) :
    afw::math::BoundedField(bbox),
    _throwOnMissing(true),
    _default(0.0), // unused
    _coaddWcs(coaddWcs),
    _elements(elements)
{}

CoaddBoundedField::CoaddBoundedField(
    afw::geom::Box2I const & bbox,
    PTR(afw::image::Wcs const) coaddWcs,
    ElementVector const & elements,
    double default_
) :
    afw::math::BoundedField(bbox),
    _throwOnMissing(false),
    _default(default_),
    _coaddWcs(coaddWcs),
    _elements(elements)
{}

double CoaddBoundedField::evaluate(afw::geom::Point2D const & position) const {
    PTR(afw::coord::Coord) coord = _coaddWcs->pixelToSky(position);
    double sum = 0.0;
    double wSum = 0.0;
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        afw::geom::Point2D transformedPosition = i->wcs->skyToPixel(*coord);
        if (afw::geom::Box2D(i->field->getBBox()).contains(transformedPosition)) {
            sum += i->weight * i->field->evaluate(transformedPosition);
            wSum += i->weight;
        }
    }
    if (wSum == 0.0) {
        if (_throwOnMissing) {
            throw LSST_EXCEPT(
                pex::exceptions::DomainErrorException,
                (boost::format("No constituent fields to evaluate at point %f, %f")
                 % position.getX() % position.getY()).str()
            );
        } else {
            return _default;
        }
    }
    return sum / wSum;
}

// ---------- Persistence -----------------------------------------------------------------------------------

// For persistence of CoaddBoundedField, we have two catalogs: the first has just one record, and contains
// the archive ID of the coadd WCS and the parameters that control missing data.  The other catalog
// has one record for each element, with fields corresponding to the data members of the Element struct.

namespace {

namespace tbl = afw::table;

// Singleton class that manages the first persistence catalog's schema and keys
class CoaddBoundedFieldPersistenceKeys1 : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<tbl::Point<int> > bboxMin;
    tbl::Key<tbl::Point<int> > bboxMax;
    tbl::Key<int> coaddWcs;
    tbl::Key<tbl::Flag> throwOnMissing;
    tbl::Key<double > default_;

    static CoaddBoundedFieldPersistenceKeys1 const & get() {
        static CoaddBoundedFieldPersistenceKeys1 const instance;
        return instance;
    }

private:
    CoaddBoundedFieldPersistenceKeys1() :
        schema(),
        bboxMin(schema.addField< tbl::Point<int> >("bbox.min", "minimum point of bounding box")),
        bboxMax(schema.addField< tbl::Point<int> >("bbox.max", "maximum point of bounding box")),
        coaddWcs(schema.addField<int>("coaddWcs", "archive ID of the coadd's WCS")),
        throwOnMissing(schema.addField<tbl::Flag>("throwOnMissing",
                                                  "whether to throw an exception on missing data")),
        default_(schema.addField<double>("default", "default value to use when throwOnMissing is False"))
    {
        schema.getCitizen().markPersistent();
    }
};

// Singleton class that manages the second persistence catalog's schema and keys
class CoaddBoundedFieldPersistenceKeys2 : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<int> field;
    tbl::Key<int> wcs;
    tbl::Key<double> weight;

    static CoaddBoundedFieldPersistenceKeys2 const & get() {
        static CoaddBoundedFieldPersistenceKeys2 const instance;
        return instance;
    }

private:
    CoaddBoundedFieldPersistenceKeys2() :
        schema(),
        field(schema.addField<int>("field", "archive ID of the BoundedField to be coadded")),
        wcs(schema.addField<int>("wcs", "archive ID of the Wcs associated with this element")),
        weight(schema.addField<double>("weight", "weight value for this element"))
    {
        schema.getCitizen().markPersistent();
    }
};

} // anonymous

class CoaddBoundedField::Factory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        CoaddBoundedFieldPersistenceKeys1 const & keys1 = CoaddBoundedFieldPersistenceKeys1::get();
        CoaddBoundedFieldPersistenceKeys2 const & keys2 = CoaddBoundedFieldPersistenceKeys2::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 2u);
        LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys1.schema);
        LSST_ARCHIVE_ASSERT(catalogs.back().getSchema() == keys2.schema);
        tbl::BaseRecord const & record1 = catalogs.front().front();
        ElementVector elements;
        elements.reserve(catalogs.back().size());
        for (tbl::BaseCatalog::const_iterator i = catalogs.back().begin(); i != catalogs.back().end(); ++i) {
            elements.push_back(
                Element(
                    archive.get<afw::math::BoundedField>(i->get(keys2.field)),
                    archive.get<afw::image::Wcs>(i->get(keys2.wcs)),
                    i->get(keys2.weight)
                )
            );
        }
        return boost::make_shared<CoaddBoundedField>(
            afw::geom::Box2I(record1.get(keys1.bboxMin), record1.get(keys1.bboxMax)),
            archive.get<afw::image::Wcs>(record1.get(keys1.coaddWcs)),
            elements,
            record1.get(keys1.default_)
        );
    }

    Factory(std::string const & name) : tbl::io::PersistableFactory(name) {}

};

namespace {

std::string getCoaddBoundedFieldPersistenceName() { return "CoaddBoundedField"; }

CoaddBoundedField::Factory registration(getCoaddBoundedFieldPersistenceName());

} // anonymous

std::string CoaddBoundedField::getPersistenceName() const { return getCoaddBoundedFieldPersistenceName(); }

std::string CoaddBoundedField::getPythonModule() const { return "lsst.meas.algorithms"; }

void CoaddBoundedField::write(OutputArchiveHandle & handle) const {
    CoaddBoundedFieldPersistenceKeys1 const & keys1 = CoaddBoundedFieldPersistenceKeys1::get();
    CoaddBoundedFieldPersistenceKeys2 const & keys2 = CoaddBoundedFieldPersistenceKeys2::get();
    tbl::BaseCatalog cat1 = handle.makeCatalog(keys1.schema);
    PTR(tbl::BaseRecord) record1 = cat1.addNew();
    record1->set(keys1.bboxMin, getBBox().getMin());
    record1->set(keys1.bboxMax, getBBox().getMax());
    record1->set(keys1.coaddWcs, handle.put(_coaddWcs));
    record1->set(keys1.default_, _default);
    handle.saveCatalog(cat1);
    tbl::BaseCatalog cat2 = handle.makeCatalog(keys2.schema);
    for (ElementVector::const_iterator i = _elements.begin(); i != _elements.end(); ++i) {
        PTR(tbl::BaseRecord) record2 = cat2.addNew();
        record2->set(keys2.field, handle.put(i->field));
        record2->set(keys2.wcs, handle.put(i->wcs));
        record2->set(keys2.weight, i->weight);
    }
    handle.saveCatalog(cat2);
}

}}} // namespace lsst::meas::algorithms


