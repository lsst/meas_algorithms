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

/*
 * Represent a PSF as for a Coadd based on the James Jee stacking
 * algorithm which was extracted from Stackfit.
 */
#include <cmath>
#include <sstream>
#include <iostream>
#include <numeric>
#include "boost/iterator/iterator_adaptor.hpp"
#include "boost/iterator/transform_iterator.hpp"
#include "ndarray/eigen.h"
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/meas/algorithms/WarpedPsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

// Struct used to simplify calculations in computeAveragePosition; lets us use
// std::accumulate instead of explicit for loop.
struct AvgPosItem {
    double wx;   // weighted x position
    double wy;   // weighted y position
    double w;    // weight value

    explicit AvgPosItem(double wx_=0.0, double wy_=0.0, double w_=0.0) : wx(wx_), wy(wy_), w(w_) {}

    // return point, assuming this is a sum of many AvgPosItems
    afw::geom::Point2D getPoint() const { return afw::geom::Point2D(wx/w, wy/w); }

    // comparison so we can sort by weights
    bool operator<(AvgPosItem const & other) const {
        return w < other.w;
    }

    AvgPosItem & operator+=(AvgPosItem const & other) {
        wx += other.wx;
        wy += other.wy;
        w += other.w;
        return *this;
    }

    AvgPosItem & operator-=(AvgPosItem const & other) {
        wx -= other.wx;
        wy -= other.wy;
        w -= other.w;
        return *this;
    }

    friend AvgPosItem operator+(AvgPosItem a, AvgPosItem const & b) { return a += b; }

    friend AvgPosItem operator-(AvgPosItem a, AvgPosItem const & b) { return a -= b; }
};

afw::geom::Point2D computeAveragePosition(
    afw::table::ExposureCatalog const & catalog,
    afw::image::Wcs const & coaddWcs,
    afw::table::Key<double> weightKey
) {
    afw::table::Key<int> goodPixKey;
    try {
        goodPixKey = catalog.getSchema()["goodpix"];
    } catch (pex::exceptions::NotFoundError &) {}
    std::vector<AvgPosItem> items;
    items.reserve(catalog.size());
    for (afw::table::ExposureCatalog::const_iterator i = catalog.begin(); i != catalog.end(); ++i) {
        afw::geom::Point2D p = coaddWcs.skyToPixel(
            *i->getWcs()->pixelToSky(
                i->getPsf()->getAveragePosition()
            )
        );
        AvgPosItem item(p.getX(), p.getY(), i->get(weightKey));
        if (goodPixKey.isValid()) {
            item.w *= i->get(goodPixKey);
        }
        item.wx *= item.w;
        item.wy *= item.w;
        items.push_back(item);
    }
    // This is a bit pessimistic - we save and sort all the weights all the time,
    // even though we'll only need them if the average position from all of them
    // is invalid.  But it makes for simpler code, and it's not that expensive
    // computationally anyhow.
    std::sort(items.begin(), items.end());
    AvgPosItem result = std::accumulate(items.begin(), items.end(), AvgPosItem());
    // If the position isn't valid (no input frames contain it), we remove frames
    // from the average until it does.
    for (
        std::vector<AvgPosItem>::iterator iter = items.begin();
        catalog.subsetContaining(result.getPoint(), coaddWcs, true).empty();
        ++iter
    ) {
        if (iter == items.end()) {
            // This should only happen if there are no inputs at all,
            // or if constituent Psfs have a badly-behaved implementation
            // of getAveragePosition().
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeError,
                "Could not find a valid average position for CoaddPsf"
            );
        }
        result -= *iter;
    }
    return result.getPoint();
}

} // anonymous

CoaddPsf::CoaddPsf(
    afw::table::ExposureCatalog const & catalog,
    afw::image::Wcs const & coaddWcs,
    std::string const & weightFieldName,
    std::string const & warpingKernelName,
    int cacheSize
) :
    _coaddWcs(coaddWcs.clone()),
    _warpingKernelName(warpingKernelName),
    _warpingControl(std::make_shared<afw::math::WarpingControl>(warpingKernelName, "", cacheSize))
{
    afw::table::SchemaMapper mapper(catalog.getSchema());
    mapper.addMinimalSchema(afw::table::ExposureTable::makeMinimalSchema(), true);

    // copy the field "goodpix", if available, for computeAveragePosition to use
    try {
        afw::table::Key<int> goodPixKey = catalog.getSchema()["goodpix"]; // auto does not work
        mapper.addMapping(goodPixKey, true);
    } catch (pex::exceptions::NotFoundError &) {}

    // copy the field specified by weightFieldName to field "weight"
    afw::table::Field<double> weightField = afw::table::Field<double>("weight", "Coadd weight");
    afw::table::Key<double> weightKey = catalog.getSchema()[weightFieldName];
    _weightKey = mapper.addMapping(weightKey, weightField);

    _catalog = afw::table::ExposureCatalog(mapper.getOutputSchema());
    for (afw::table::ExposureCatalog::const_iterator i = catalog.begin(); i != catalog.end(); ++i) {
         PTR(afw::table::ExposureRecord) record = _catalog.getTable()->makeRecord();
         record->assign(*i, mapper);
         _catalog.push_back(record);
    }
    _averagePosition = computeAveragePosition(_catalog, *_coaddWcs, _weightKey);
}

PTR(afw::detection::Psf) CoaddPsf::clone() const {
    return std::make_shared<CoaddPsf>(*this);
}

PTR(afw::detection::Psf) CoaddPsf::resized(int width, int height) const {
    // Not implemented for WarpedPsf
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
    }

// Read all the images from the Image Vector and return the BBox in xy0 offset coordinates

afw::geom::Box2I getOverallBBox(std::vector<PTR(afw::image::Image<double>)> const & imgVector) {

    afw::geom::Box2I bbox;
    // Calculate the box which will contain them all
    for (unsigned int i = 0; i < imgVector.size(); i ++) {
        PTR(afw::image::Image<double>) componentImg = imgVector[i];
        afw::geom::Box2I cBBox = componentImg->getBBox();
        bbox.include(cBBox); // JFB: this works even on empty bboxes
    }
    return bbox;
}


// Read all the images from the Image Vector and add them to image

void addToImage(
    PTR(afw::image::Image<double>) image,
    std::vector<PTR(afw::image::Image<double>)> const & imgVector,
    std::vector<double> const & weightVector
) {
    assert(imgVector.size() == weightVector.size());
    for (unsigned int i = 0; i < imgVector.size(); i ++) {
        PTR(afw::image::Image<double>) componentImg = imgVector[i];
        double weight = weightVector[i];
        double sum = componentImg->getArray().asEigen().sum();

        // Now get the portion of the component image which is appropriate to add
        // If the default image size is used, the component is guaranteed to fit,
        // but not if a size has been specified.
        afw::geom::Box2I cBBox = componentImg->getBBox();
        afw::geom::Box2I overlap(cBBox);
        overlap.clip(image->getBBox());
        // JFB: A subimage view of the image we want to add to, containing only the overlap region.
        afw::image::Image<double> targetSubImage(*image, overlap);
        // JFB: A subimage view of the image we want to add from, containing only the overlap region.
        afw::image::Image<double> cSubImage(*componentImg, overlap);
        targetSubImage.scaledPlus(weight/sum, cSubImage);
    }
}


afw::geom::Box2I CoaddPsf::doComputeBBox(
    afw::geom::Point2D const & ccdXY,
    afw::image::Color const & color
) const {
    afw::table::ExposureCatalog subcat = _catalog.subsetContaining(ccdXY, *_coaddWcs, true);
    if (subcat.empty()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("Cannot compute BBox at point %s; no input images at that point.")
             % ccdXY).str());
    }

    afw::geom::Box2I ret;
    for (auto const & exposureRecord : subcat) {
        PTR(afw::geom::XYTransform) xytransform(
            new afw::image::XYTransformFromWcsPair(_coaddWcs, exposureRecord.getWcs()));
        WarpedPsf warpedPsf = WarpedPsf(exposureRecord.getPsf(), xytransform, _warpingControl);
        afw::geom::Box2I componentBBox = warpedPsf.computeBBox(ccdXY, color);
        ret.include(componentBBox);
    }

    return ret;
}

PTR(afw::detection::Psf::Image) CoaddPsf::doComputeKernelImage(
    afw::geom::Point2D const & ccdXY,
    afw::image::Color const & color
) const {
    // Get the subset of expoures which contain our coordinate within their validPolygons.
    afw::table::ExposureCatalog subcat = _catalog.subsetContaining(ccdXY, *_coaddWcs, true);
    if (subcat.empty()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("Cannot compute CoaddPsf at point %s; no input images at that point.")
             % ccdXY).str()
        );
    }
    double weightSum = 0.0;

    // Read all the Psf images into a vector.  The code is set up so that this can be done in chunks,
    // with the image modified to accomodate
    // However, we currently read all of the images.
    std::vector<PTR(afw::image::Image<double>)> imgVector;
    std::vector<double> weightVector;

    for (auto const & exposureRecord : subcat) {
        PTR(afw::geom::XYTransform) xytransform(
            new afw::image::XYTransformFromWcsPair(_coaddWcs, exposureRecord.getWcs())
        );
        PTR(afw::image::Image<double>) componentImg;
        try {
            WarpedPsf warpedPsf = WarpedPsf(exposureRecord.getPsf(), xytransform, _warpingControl);
            componentImg = warpedPsf.computeKernelImage(ccdXY, color);
        } catch (pex::exceptions::RangeError & exc) {
            LSST_EXCEPT_ADD(exc, (boost::format("Computing WarpedPsf kernel image for id=%d") %
                                  exposureRecord.getId()).str());
            throw exc;
        }
        imgVector.push_back(componentImg);
        weightSum += exposureRecord.get(_weightKey);
        weightVector.push_back(exposureRecord.get(_weightKey));
    }

    afw::geom::Box2I bbox = getOverallBBox(imgVector);

    // create a zero image of the right size to sum into
    PTR(afw::detection::Psf::Image) image = std::make_shared<afw::detection::Psf::Image>(bbox);
    *image = 0.0;
    addToImage(image, imgVector, weightVector);
    *image /= weightSum;
    return image;
}

int CoaddPsf::getComponentCount() const {
    return _catalog.size();
}

CONST_PTR(afw::detection::Psf) CoaddPsf::getPsf(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getPsf();
}

CONST_PTR(afw::image::Wcs) CoaddPsf::getWcs(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getWcs();
}

CONST_PTR(afw::geom::polygon::Polygon) CoaddPsf::getValidPolygon(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getValidPolygon();
}

double CoaddPsf::getWeight(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "index of CoaddPsf component out of range");
    }
    return _catalog[index].get(_weightKey);
}

afw::table::RecordId CoaddPsf::getId(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getId();
}

afw::geom::Box2I CoaddPsf::getBBox(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeError, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getBBox();
}

// ---------- Persistence -----------------------------------------------------------------------------------

// For persistence of CoaddPsf, we have two catalogs: the first has just one record, and contains
// the archive ID of the coadd WCS, the size of the warping cache, the name of the warping kernel,
// and the average position.  The latter is simply the ExposureCatalog.

namespace {

namespace tbl = afw::table;

// Singleton class that manages the first persistence catalog's schema and keys
class CoaddPsfPersistenceHelper {
public:
    tbl::Schema schema;
    tbl::Key<int> coaddWcs;
    tbl::Key<int> cacheSize;
    tbl::PointKey<double> averagePosition;
    tbl::Key<std::string> warpingKernelName;

    static CoaddPsfPersistenceHelper const & get() {
        static CoaddPsfPersistenceHelper const instance;
        return instance;
    }

private:
    CoaddPsfPersistenceHelper() :
        schema(),
        coaddWcs(schema.addField<int>("coaddwcs", "archive ID of the coadd's WCS")),
        cacheSize(schema.addField<int>("cachesize", "size of the warping cache")),
        averagePosition(tbl::PointKey<double>::addFields(
            schema, "avgpos", "PSF accessors default position", "pixel"
        )),
        warpingKernelName(schema.addField<std::string>("warpingkernelname", "warping kernel name", 32))
    {
        schema.getCitizen().markPersistent();
    }
};

} // anonymous

class CoaddPsf::Factory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        if (catalogs.size() == 1u) {
            // Old CoaddPsfs were saved in only one catalog, because we didn't
            // save the warping parameters and average position, and we could
            // save the coadd Wcs in a special final record.
            return readV0(archive, catalogs);
        }
        LSST_ARCHIVE_ASSERT(catalogs.size() == 2u);
        CoaddPsfPersistenceHelper const & keys1 = CoaddPsfPersistenceHelper::get();
        LSST_ARCHIVE_ASSERT(catalogs.front().getSchema() == keys1.schema);
        tbl::BaseRecord const & record1 = catalogs.front().front();
        return PTR(CoaddPsf)(
            new CoaddPsf(
                tbl::ExposureCatalog::readFromArchive(archive, catalogs.back()),
                archive.get<afw::image::Wcs>(record1.get(keys1.coaddWcs)),
                record1.get(keys1.averagePosition),
                record1.get(keys1.warpingKernelName),
                record1.get(keys1.cacheSize)
            )
        );
    }


    // Backwards compatibility for files saved before meas_algorithms commit
    // 53e61fae (7/10/2013).  Prior to that change, the warping configuration
    // and the average position were not saved at all, making it impossible to
    // reconstruct the average position exactly, but it's better to
    // approximate than to fail completely.
    std::shared_ptr<tbl::io::Persistable>
    readV0(InputArchive const & archive, CatalogVector const & catalogs) const {
        auto internalCat = tbl::ExposureCatalog::readFromArchive(archive, catalogs.front());
        // Coadd WCS is stored in a special last record.
        auto coaddWcs = internalCat.back().getWcs();
        internalCat.pop_back();
        // Attempt to reconstruct the average position.  We can't do this
        // exactly, since the catalog we saved isn't the same one that was
        // used to compute the original average position.
        tbl::Key<double> weightKey;
        try {
             weightKey = internalCat.getSchema()["weight"];
        } catch (pex::exceptions::NotFoundError &) {}
        auto averagePos = computeAveragePosition(internalCat, *coaddWcs, weightKey);
        return std::shared_ptr<CoaddPsf>(new CoaddPsf(internalCat, coaddWcs, averagePos));
    }

    Factory(std::string const & name) : tbl::io::PersistableFactory(name) {}

};

namespace {

std::string getCoaddPsfPersistenceName() { return "CoaddPsf"; }

CoaddPsf::Factory registration(getCoaddPsfPersistenceName());

} // anonymous

std::string CoaddPsf::getPersistenceName() const { return getCoaddPsfPersistenceName(); }

std::string CoaddPsf::getPythonModule() const { return "lsst.meas.algorithms"; }

void CoaddPsf::write(OutputArchiveHandle & handle) const {
    CoaddPsfPersistenceHelper const & keys1 = CoaddPsfPersistenceHelper::get();
    tbl::BaseCatalog cat1 = handle.makeCatalog(keys1.schema);
    PTR(tbl::BaseRecord) record1 = cat1.addNew();
    record1->set(keys1.coaddWcs, handle.put(_coaddWcs));
    record1->set(keys1.cacheSize, _warpingControl->getCacheSize());
    record1->set(keys1.averagePosition, _averagePosition);
    record1->set(keys1.warpingKernelName, _warpingKernelName);
    handle.saveCatalog(cat1);
    _catalog.writeToArchive(handle, false);
}

CoaddPsf::CoaddPsf(
    afw::table::ExposureCatalog const & catalog,
    PTR(afw::image::Wcs const) coaddWcs,
    afw::geom::Point2D const & averagePosition,
    std::string const & warpingKernelName,
    int cacheSize
) :
    _catalog(catalog), _coaddWcs(coaddWcs), _weightKey(_catalog.getSchema()["weight"]),
    _averagePosition(averagePosition), _warpingKernelName(warpingKernelName),
    _warpingControl(new afw::math::WarpingControl(warpingKernelName, "", cacheSize))
{}

}}} // namespace lsst::meas::algorithms


