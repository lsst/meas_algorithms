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

afw::geom::Point2D computeAveragePosition(
    afw::table::ExposureCatalog const & catalog,
    afw::image::Wcs const & coaddWcs,
    afw::table::Key<double> weightKey
) {
    double avgX = 0.0;
    double avgY = 0.0;
    double wSum = 0.0;
    for (afw::table::ExposureCatalog::const_iterator i = catalog.begin(); i != catalog.end(); ++i) {
        afw::geom::Point2D p = coaddWcs.skyToPixel(
            *i->getWcs()->pixelToSky(
                i->getPsf()->getAveragePosition()
            )
        );
        double w = i->get(weightKey);
        avgX += p.getX() * w;
        avgY += p.getY() * w;
        wSum += w;
    }
    avgX /= wSum;
    avgY /= wSum;
    return afw::geom::Point2D(avgX, avgY);
}

} // anonymous

CoaddPsf::CoaddPsf(
    afw::table::ExposureCatalog const & catalog,
    afw::image::Wcs const & coaddWcs,
    std::string const & weightFieldName
) {
    _coaddWcs = coaddWcs.clone();
    afw::table::SchemaMapper mapper(catalog.getSchema());
    mapper.addMinimalSchema(afw::table::ExposureTable::makeMinimalSchema(), true);
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
    return boost::make_shared<CoaddPsf>(*this);
}


// Read all the images from the Image Vector and return the BBox in xy0 offset coordinates

afw::geom::Box2I getOverallBBox(std::vector<PTR(afw::image::Image<double>)> const & imgVector) {

    afw::geom::Box2I bbox;
    // Calculate the box which will contain them all
    for (unsigned int i = 0; i < imgVector.size(); i ++) {
        PTR(afw::image::Image<double>) componentImg = imgVector[i];
        afw::geom::Box2I cBBox = componentImg->getBBox(afw::image::PARENT);
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
        afw::geom::Box2I cBBox = componentImg->getBBox(afw::image::PARENT);
        afw::geom::Box2I overlap(cBBox);
        overlap.clip(image->getBBox(afw::image::PARENT));
        // JFB: A subimage view of the image we want to add to, containing only the overlap region.
        afw::image::Image<double> targetSubImage(*image, overlap, afw::image::PARENT);
        // JFB: A subimage view of the image we want to add from, containing only the overlap region.
        afw::image::Image<double> cSubImage(*componentImg, overlap, afw::image::PARENT);
        targetSubImage.scaledPlus(weight/sum, cSubImage);
    }
}


PTR(afw::detection::Psf::Image) CoaddPsf::doComputeKernelImage(
    afw::geom::Point2D const& ccdXY,
    afw::image::Color const& color
) const {
    // get the subset of exposures which contain our coordinate
    afw::table::ExposureCatalog subcat = _catalog.subsetContaining(ccdXY, *_coaddWcs);
    if (subcat.empty()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
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

    for (afw::table::ExposureCatalog::const_iterator i = subcat.begin(); i != subcat.end(); ++i) {
        PTR(afw::geom::XYTransform) xytransform(
            new afw::image::XYTransformFromWcsPair(_coaddWcs, i->getWcs())
        );
        WarpedPsf warpedPsf = WarpedPsf(i->getPsf(), xytransform);
        PTR(afw::image::Image<double>) componentImg = warpedPsf.computeKernelImage(ccdXY, color);
        imgVector.push_back(componentImg);
        weightSum += i->get(_weightKey);
        weightVector.push_back(i->get(_weightKey));
    }

    afw::geom::Box2I bbox = getOverallBBox(imgVector);

    // create a zero image of the right size to sum into
    PTR(afw::detection::Psf::Image) image = boost::make_shared<afw::detection::Psf::Image>(bbox);
    *image = 0.0;
    addToImage(image, imgVector, weightVector);
    *image /= weightSum;
    return image;
}

/**
 * getComponentCount() - get the number of component Psf's in this CoaddPsf
 */
int CoaddPsf::getComponentCount() const {
    return _catalog.size();
}

/**
 * getPsf - get the Psf of the component at position index
 */
CONST_PTR(afw::detection::Psf) CoaddPsf::getPsf(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getPsf();
}

/**
 * getWcs - get the Wcs of the component at position index
 */
CONST_PTR(afw::image::Wcs) CoaddPsf::getWcs(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getWcs();
}

/**
 * getWeight - get the coadd weight of the component at position index
 */
double CoaddPsf::getWeight(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
    }
    return _catalog[index].get(_weightKey);
}

/**
 * getId - get the long id of the component at position index
 */
afw::table::RecordId CoaddPsf::getId(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getId();
}

afw::geom::Box2I CoaddPsf::getBBox(int index) {
    if (index < 0 || index > getComponentCount()) {
        throw LSST_EXCEPT(pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
    }
    return _catalog[index].getBBox();
}

// ---------- Persistence -----------------------------------------------------------------------------------

// For persistence of CoaddPsf, we append one record to the ExposureCatalog to hold the coadd Wcs;
// this saves us from having to add another catalog in persistence just to hold the Wcs archive ID
// (a single integer).  That in turn saves us from having to write a separate FITS HDU just for that
// integer.

class CoaddPsf::Factory : public afw::table::io::PersistableFactory {
public:

    virtual PTR(afw::table::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        return PTR(CoaddPsf)(
            new CoaddPsf(afw::table::ExposureCatalog::readFromArchive(archive, catalogs.front()))
        );
    }

    Factory(std::string const & name) : afw::table::io::PersistableFactory(name) {}

};

namespace {

std::string getCoaddPsfPersistenceName() { return "CoaddPsf"; }

CoaddPsf::Factory registration(getCoaddPsfPersistenceName());

} // anonymous

std::string CoaddPsf::getPersistenceName() const { return getCoaddPsfPersistenceName(); }

std::string CoaddPsf::getPythonModule() const { return "lsst.meas.algorithms"; }

void CoaddPsf::write(OutputArchiveHandle & handle) const {
    afw::table::ExposureCatalog catCopy(_catalog); // copy shares record data, but not container
    PTR(afw::table::ExposureRecord) coaddWcsRecord = catCopy.addNew();
    coaddWcsRecord->setWcs(_coaddWcs);
    catCopy.writeToArchive(handle, false);
}

CoaddPsf::CoaddPsf(afw::table::ExposureCatalog const & catalog) :
    _catalog(catalog), _coaddWcs(catalog.back().getWcs()), _weightKey(_catalog.getSchema()["weight"])
{
    _catalog.pop_back();
    _averagePosition = computeAveragePosition(_catalog, *_coaddWcs, _weightKey);
}

}}} // namespace lsst::meas::algorithms


