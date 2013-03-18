// -*- lsst-c++ -*-
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

#if !defined(LSST_MEAS_ALGORITHMS_COADDPSF_H)
#define LSST_MEAS_ALGORITHMS_COADDPSF_H

#include <boost/make_shared.hpp>
#include "lsst/base.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/types.h"
#include "lsst/afw/geom/Box.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief CoaddPsf is the Psf descendent to be used for Coadd images.
 *  It incorporates the logic of James Jee's Stackfit algorithm
 *  for estimating the Psf of the stack of images (calexps)
 *  weighted by a given weighting vector
 *
 *  The user is expected to supply an Exposure Catalog
 *  which describes the images and Psf's which are to be stacked
 */

class CoaddPsf : public afw::table::io::PersistableFacade<CoaddPsf>, public afw::detection::Psf {
public:
    typedef PTR(CoaddPsf) Ptr;
    typedef CONST_PTR(CoaddPsf) ConstPtr;

    /**
     * @brief constructors for a CoaddPsf - The ExposureCatalog contains info about each visit/ccd in Coadd
     *                                      Must be provided on the constructor, and cannot be changed.
     *
     * Parameters:  ExposureCatalog containing the refid, bbox, wcs, psf and weight for each ccd/visit
     *              weightFieldName is optional.  Field is assumed to be a double of name "weight".
     */
    explicit CoaddPsf(afw::table::ExposureCatalog const & catalog, afw::image::Wcs const & coaddWcs,
                        std::string const & weightFieldName = "weight");

    virtual PTR(afw::detection::Psf) clone() const {
        return boost::make_shared<CoaddPsf>(*this);
    }

    /**
     * @brief getCoaddWcs - Wcs of the coadd - this is specified on the constructor
     */
    CONST_PTR(afw::image::Wcs) getCoaddWcs() {
        return _coaddWcs;
    }

    /**
     * @brief getComponentCount() - get the number of component Psf's in this CoaddPsf
     */
    int getComponentCount() const;

    /**
     * @brief getPsf - get the Psf of the component at position index
     */
    CONST_PTR(afw::detection::Psf) getPsf(int index);

    /**
     * @brief getWcs - get the Wcs of the component at position index
     */
    CONST_PTR(afw::image::Wcs) getWcs(int index);

    /**
     * @brief getWeight - get the coadd weight of the component at position index
     */
    double getWeight(int index);

    /**
     * @brief getId - get the long id of the component at position index
     */
    afw::table::RecordId getId(int index);

    /**
     * @brief getBBox - the bounding box for this component in its own Wcs
     */
    afw::geom::Box2I getBBox(int index);
    /**
     *  @brief Return true if the CoaddPsf persistable (always true).
     *
     *  While it's actually possible to construct a CoaddPsf that isn't persistable (because its nested
     *  Psfs and Wcss are not persistable) in artificial situations, in realistic situations it's
     *  pretty much impossible, because persistence is a necessary part of how CoaddPsfs are built.
     *  And it's simpler and much faster if we just always return true, rather than loop over the
     *  elements and check each one.
     */
    virtual bool isPersistable() const { return true; }

    // Factory used to read CoaddPsf from an InputArchive; defined only in the source file.
    class Factory;

protected:

    PTR(afw::detection::Psf::Image) doComputeKernelImage(
        afw::geom::Point2D const& ccdXY,
        afw::image::Color const& color
    ) const;

    // See afw::table::io::Persistable::getPersistenceName
    virtual std::string getPersistenceName() const;

    // See afw::table::io::Persistable::write
    virtual void write(OutputArchiveHandle & handle) const;

    // Used by persistence; the coadd Wcs is expected to be in the last record of the catalog.
    explicit CoaddPsf(afw::table::ExposureCatalog const & catalog);

private:

    afw::table::ExposureCatalog _catalog;
    CONST_PTR(afw::image::Wcs) _coaddWcs;
    afw::table::Key<double> _weightKey;
};

}}}


#endif
