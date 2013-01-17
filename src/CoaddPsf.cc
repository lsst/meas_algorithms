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
    // get the WCS coord of the requested point <pgee> do we need to add the image xy?
    // find a subcat of images which contain this coord
    PTR(afw::coord::Coord) x = _coaddWcs->pixelToSky(ccdXY);
    afw::coord::Coord const & coord = *x;
    afw::table::ExposureCatalog subcat = _catalog.findContains(coord);
    afw::table::Key<double> weightKey = subcat.getSchema()["weight"];

    // create a zero image of the right size to sum into
    lsst::afw::detection::Psf::Image::Ptr image = boost::make_shared<lsst::afw::detection::Psf::Image>(size);
    *image *= 0.0;

    for (lsst::afw::table::ExposureCatalog::const_iterator i = subcat.begin(); i != subcat.end(); ++i) {

        // This is the code we used to use:
/*
        lsst::afw::geom::Box2I bbox =  r.getBBox();
        double xrel = ccdXY.getX() - bbox.getBeginX();
        double yrel = ccdXY.getY() - bbox.getBeginY();
*/
        // this code is a guess at what is needed:  my best guess is that we need to move the ccdXY point to the
        // corresponding point on the ccd of the original Psf.  Note sure what Kendrick's code will do here.
        // There is also the issue of the xy of the image.  Not sure if that affects the coordinate system of the Psf
        double weight = i->get(weightKey);
        CONST_PTR(lsst::afw::detection::Psf) psf = (i->getPsf());
        afw::geom::Point2D point = i->getWcs()->skyToPixel(coord);
        PTR(afw::image::Image<double>) ii = psf->computeImage(point, size, true, true);
        double sum = ii->getArray().asEigen().sum();
        std::cout << "sum = " << sum << "\n";  //*image *= 1/sum;
        image->scaledPlus(weight/sum, *ii);
    }

    // Not really sure what normalizePeak should do.  For now, set the max value to 1.0
    if (normalizePeak) {
        double max = image->getArray().asEigen().maxCoeff();
        *image *= 1.0/max;
    }
    return image;
}

int CoaddPsf::getComponentCount() const {
    return _catalog.size();
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
