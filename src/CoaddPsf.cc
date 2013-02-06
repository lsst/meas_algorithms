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
 
/*!
 * @brief Represent a PSF as for a Coadd based on the James Jee stacking
 * algorithm which was extracted from Stackfit.
 *
 * Note that this Psf subclass only support computeImage, not the 
 * parameterization methodes defined on its super class.  In that sense,
 * it is not a true subclass.
 *
 * @file
 *
 * @ingroup algorithms
 */
#include <cmath>
#include <sstream>
#include <iostream>
#include <numeric>
#include "boost/iterator/iterator_adaptor.hpp"
#include "boost/iterator/transform_iterator.hpp"
#include "lsst/base.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/ImageUtils.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/PcaPsf.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/image/XYTransform.h"
#include "lsst/afw/detection/WarpedPsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
  * @brief CoaddPsf class
  *
  */


CoaddPsf::CoaddPsf(afw::table::ExposureCatalog const & catalog, afw::image::Wcs const & coaddWcs, std::string const & weightFieldName) {

    _coaddWcs = coaddWcs.clone();
    afw::table::SchemaMapper mapper(catalog.getSchema());
    mapper.addMinimalSchema(afw::table::ExposureTable::makeMinimalSchema(), true);
    afw::table::Field<double> weightField = afw::table::Field<double>("weight", "Coadd weight");
//    mapper.addOutputField(weightField);
    afw::table::Key<double> weightKey = catalog.getSchema()[weightFieldName];
    mapper.addMapping(weightKey, weightField);

    _catalog = afw::table::ExposureCatalog(mapper.getOutputSchema());
    for (lsst::afw::table::ExposureCatalog::const_iterator i = catalog.begin(); i != catalog.end(); ++i) {
         PTR(lsst::afw::table::ExposureRecord) record = _catalog.getTable()->makeRecord();
         record->assign(*i, mapper);
         _catalog.push_back(record);
    }
}

// Read all the images from the Image Vector and return the BBox in xy0 offset coordinates

afw::geom::Box2I getOverallBBox( std::vector<PTR(afw::image::Image<double>)> imgVector) {

    afw::geom::Box2I bbox;
    // Calculate the box which will contain them all
    for (unsigned int i = 0; i < imgVector.size(); i ++) {
        PTR(afw::image::Image<double>) componentImg = imgVector[i];
        afw::geom::Box2I cBBox = componentImg->getBBox();
        cBBox.shift(afw::geom::Extent2I(componentImg->getXY0()));
        if (i == 0) {
            bbox = cBBox;
        }
        else {
            bbox.include(cBBox);
        }
        //std::cout << "After image " << i << ", cbbox = " << cBBox << ", bbox = " << bbox << std::endl;
    }
    return bbox;
}


// Read all the images from the Image Vector and add them to image

void addToImage(PTR(afw::image::Image<double>) image, std::vector<PTR(afw::image::Image<double>)> imgVector, std::vector<double> weightVector) {
    assert(imgVector.size() == weightVector.size());
    for (unsigned int i = 0; i < imgVector.size(); i ++) {
        PTR(afw::image::Image<double>) componentImg = imgVector[i];
        double weight = weightVector[i];
        double sum = componentImg->getArray().asEigen().sum();
        componentImg->getArray()[ndarray::view(0,componentImg->getHeight())(0,componentImg->getWidth())] *= (weight/sum);
        std::stringstream ss;
        ss << i;
        componentImg->writeFits("comp" + ss.str() + ".fits");

        // Now get the portion of the component image which is appropriate to add
        // If the default image size is used, the component is guaranteed to fit,
        // but not if a size has been specified.
        afw::geom::Box2I cBBox = componentImg->getBBox();
        // std::cout << "Original: " << cBBox  << ", XY0 = " << componentImg->getXY0() << "\n";
        afw::geom::Extent2I relativeOffset = afw::geom::Extent2I(image->getXY0() - componentImg->getXY0());  // this has to be added to get us into the coadd image frame
        // std::cout << relativeOffset << "\n";
        cBBox.shift(relativeOffset);
        cBBox.clip(image->getBBox());    // This is now the destination bbox (in the image frame)
        int minx = cBBox.getBeginX();
        int miny = cBBox.getBeginY();
        int maxx = cBBox.getEndX();
        int maxy = cBBox.getEndY();
        int dx = relativeOffset.getX();
        int dy = relativeOffset.getY();
        // std::cout << "clipped component in image frame: " << cBBox  << ", dx = " << dx << ", dy=" << dy << "\n";
        // std::cout << "From component image: " << minx << ","<<maxx<<":" << miny << "," << maxy << " --> image: " << minx-dx << ","<<maxx-dx<<":" << miny-dy << "," << maxy-dy << "\n";
        image->writeFits("sum_pre" + ss.str() + ".fits");
        image->getArray()[ndarray::view(miny-dy,maxy-dy)(minx-dx,maxx-dx)] += componentImg->getArray()[ndarray::view(miny,maxy)(minx,maxx)];
    }
}

    /**
     *  @brief doComputeImage produces an estimate of the Psf at the given location, relative to the psf spatial modelj
     *   Still need to implement nomaliziePeak and distort 
     */

lsst::afw::detection::Psf::Image::Ptr CoaddPsf::doComputeImage(lsst::afw::image::Color const& color,
                                  lsst::afw::geom::Point2D const& ccdXY,
                                  lsst::afw::geom::Extent2I const& size,
                                  bool normalizePeak,
                                  bool distort
                                 ) const {
    // get the subset of exposures which contain our coordinate
    PTR(afw::coord::Coord) x = _coaddWcs->pixelToSky(ccdXY);
    afw::coord::Coord const & coord = *x;
    afw::table::ExposureCatalog subcat = _catalog.findContains(coord);
    afw::table::Key<double> weightKey = subcat.getSchema()["weight"];

    // this is the default xy0 and box for the coadd image.  If size is zero, it may be recalculated
    // bbox is the bounding box of our coadd image
    afw::geom::Box2I bbox = afw::geom::Box2I(afw::geom::Point2I((size.getX()-1)/2, (size.getY()-1)/2), size);

    // Read all the Psf images into a vector.  The code is set up so that this can be done in chunks, with the image modified to accomodate
    // However, we currently read all of the images.
    std::vector<PTR(afw::image::Image<double>)> imgVector;
    std::vector<double> weightVector;

    for (lsst::afw::table::ExposureCatalog::const_iterator i = subcat.begin(); i != subcat.end(); ++i) {
        boost::shared_ptr<afw::detection::Psf> psf = i->getPsf()->clone();
        boost::shared_ptr<afw::image::XYTransform> xytransform = boost::shared_ptr<afw::image::XYTransform> ( new afw::image::XYTransformFromWcsPair(_coaddWcs, i->getWcs()));
        afw::detection::WarpedPsf warpedPsf = afw::detection::WarpedPsf(psf, xytransform);
        PTR(afw::image::Image<double>) componentImg = warpedPsf.computeImage(ccdXY, bbox.getDimensions(), true, false);
        imgVector.push_back(componentImg);
        weightVector.push_back(i->get(weightKey));
    }

    // If size is not specified, calculate the box which will contain them all
    if(bbox.getWidth() == 0 || bbox.getHeight() ==0) {
        bbox = getOverallBBox(imgVector);
    }

    // std::cout << "Final bbox: "  << bbox << "\n";

    // create a zero image of the right size to sum into
    lsst::afw::detection::Psf::Image::Ptr image = boost::make_shared<lsst::afw::detection::Psf::Image>(bbox.getDimensions());
    *image *= 0.0;
    image->setXY0(bbox.getBegin());
    // std::cout << "Final bbox: "  << bbox << ", XY0 = " << image->getXY0() << "\n";
    addToImage(image, imgVector, weightVector);  

    // Not really sure what normalizePeak should do.  For now, set the max value to 1.0
    if (normalizePeak) {
        double max = image->getArray().asEigen().maxCoeff();
        *image *= 1.0/max;
    }
    image->setXY0(bbox.getBegin());
    //std::cout << "Image size = " << image->getDimensions() << ", image xy0 = " <<image->getXY0() << "\n";
    //image->writeFits("sum.fits");
    return image;
}



/**
 * @brief getComponentCount() - get the number of component Psf's in this CoaddPsf
 */
int CoaddPsf::getComponentCount() const {
    return _catalog.size();
}

/**
 * @brief getPsf - get the Psf of the component at position index
 */
afw::detection::Psf::ConstPtr CoaddPsf::getPsf(int index) {
    int count = 0;
    for (lsst::afw::table::ExposureCatalog::const_iterator i = _catalog.begin(); i != _catalog.end(); ++i) {
        if (count++ == index) return i->getPsf();
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
}

/**
 * @brief getWcs - get the Wcs of the component at position index
 */
afw::image::Wcs::ConstPtr CoaddPsf::getWcs(int index) {
    int count = 0;
    for (lsst::afw::table::ExposureCatalog::const_iterator i = _catalog.begin(); i != _catalog.end(); ++i) {
        if (count++ == index) return i->getWcs();
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
}

/**
 * @brief getWeight - get the coadd weight of the component at position index
 */
int CoaddPsf::getWeight(int index) {
    int count = 0;
    for (lsst::afw::table::ExposureCatalog::const_iterator i = _catalog.begin(); i != _catalog.end(); ++i) {
        if (count++ == index) {
            afw::table::Key<double> weightKey = _catalog.getSchema()["weight"];
            return i->get(weightKey);
        }
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
}

/**
 * @brief getId - get the long id of the component at position index
 */
long CoaddPsf::getId(int index) {
    int count = 0;
    for (lsst::afw::table::ExposureCatalog::const_iterator i = _catalog.begin(); i != _catalog.end(); ++i) {
        if (count++ == index) return i->getId();
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
}

afw::geom::Box2I CoaddPsf::getBBox(int index) {
    int count = 0;
    for (lsst::afw::table::ExposureCatalog::const_iterator i = _catalog.begin(); i != _catalog.end(); ++i) {
        if (count++ == index) return i->getBBox();
    }
    throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "index of CoaddPsf component out of range");
}

/*
//
// We need to make an instance here so as to register it with createPSF
//
// \cond
namespace {
    volatile bool isInstance =
        lsst::afw::detection::Psf::registerMe<CoaddPsf, PTR(lsst::afw::math::Kernel)>("COADD");
}
*/
// ---------- Persistence -----------------------------------------------------------------------------------

class CoaddPsf::Factory : public afw::table::io::PersistableFactory {
public:

    virtual PTR(afw::table::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        return PTR(CoaddPsf)(
            new CoaddPsf(
                afw::table::ExposureCatalog::readFromArchive(archive, catalogs.front()),
                true // invoke private persistence constructor
            )
        );
    }

    Factory(std::string const & name) : afw::table::io::PersistableFactory(name) {}

};

namespace {

std::string getCoaddPsfPersistenceName() { return "CoaddPsf"; }

CoaddPsf::Factory registration(getCoaddPsfPersistenceName());

} // anonymous

std::string CoaddPsf::getPersistenceName() const { return getCoaddPsfPersistenceName(); }

void CoaddPsf::write(OutputArchiveHandle & handle) const {
    _catalog.writeToArchive(handle, false);
}

CoaddPsf::CoaddPsf(afw::table::ExposureCatalog const & catalog, bool) : _catalog(catalog) {}

}}} // namespace lsst::meas::algorithms



// \endcond
