// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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
#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/table/Exposure.h"
#include "lsst/afw/table/types.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/polygon/Polygon.h"
#include "lsst/afw/math/warpExposure.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 *  @brief CoaddPsf is the Psf derived to be used for non-PSF-matched Coadd images.
 *
 *  It incorporates the logic of James Jee's Stackfit algorithm for estimating the
 *  Psf of coadd by coadding the images of the Psf models of each input exposure.
 */
class CoaddPsf : public afw::table::io::PersistableFacade<CoaddPsf>, public ImagePsf {
public:

    /**
     * @brief Main constructors for CoaddPsf
     *
     * The ExposureCatalog contains info about each visit/ccd in Coadd; this must be provided to the
     * constructor, and cannot be changed.
     *
     * @param[in] catalog           ExposureCatalog containing the id, bbox, wcs, psf and weight for
     *                              each ccd/visit.  This is usually the same catalog as the "ccds"
     *                              catalog in the coadd Exposure's CoaddInputs.
     * @param[in] coaddWcs          Wcs for the coadd.
     * @param[in] weightFieldName   Field name that contains the weight of the exposure in the coadd;
     *                              defaults to "weight".
     * @param[in] warpingKernelName Name of warping kernel
     * @param[in] cacheSize         Warping kernel cache size
     */
    explicit CoaddPsf(
        afw::table::ExposureCatalog const & catalog,
        afw::image::Wcs const & coaddWcs,
        std::string const & weightFieldName = "weight",
        std::string const & warpingKernelName="lanczos3",
        int cacheSize=10000
    );

    /// Polymorphic deep copy.  Usually unnecessary, as Psfs are immutable.
    virtual PTR(afw::detection::Psf) clone() const;

    /**
     *  @brief Return the average of the positions of the stars that went into this Psf.
     *
     *  For CoaddPsf, this is calculated as the weighted average of the average positions
     *  of all the component Psfs.
     */
    virtual afw::geom::Point2D getAveragePosition() const { return _averagePosition; }

    /// Return the Wcs of the coadd (defines the coordinate system of the Psf).
    PTR(afw::image::Wcs const) getCoaddWcs() { return _coaddWcs; }

    /// Return the number of component Psfs in this CoaddPsf
    int getComponentCount() const;

    /// Return the Psf of the component image at index
    CONST_PTR(afw::detection::Psf) getPsf(int index);

    /// Return the Wcs of the component image at index
    CONST_PTR(afw::image::Wcs) getWcs(int index);

    /// Return the weight of the component image at index
    double getWeight(int index);

    /// Return the exposure ID of the component image at index
    afw::table::RecordId getId(int index);

    /// Return the bounding box (in component image Pixel coordinates) of the component image at index
    afw::geom::Box2I getBBox(int index);

    /// Return the valid Polygon (in component image Pixel coordinates) of the component image at index
    CONST_PTR(afw::geom::polygon::Polygon) getValidPolygon(int index);

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
        afw::geom::Point2D const & ccdXY,
        afw::image::Color const & color
    ) const;

    // See afw::table::io::Persistable::getPersistenceName
    virtual std::string getPersistenceName() const;

    // See afw::table::io::Persistable::getPythonModule
    virtual std::string getPythonModule() const;

    // See afw::table::io::Persistable::write
    virtual void write(OutputArchiveHandle & handle) const;

    // Used by persistence only
    explicit CoaddPsf(
        afw::table::ExposureCatalog const & catalog, ///< Unpersisted catalog
        PTR(afw::image::Wcs const) coaddWcs,         ///< WCS for the coadd
        afw::geom::Point2D const & averagePosition,  ///< Default position for accessors
        std::string const & warpingKernelName="lanczos3",    ///< Warping kernel name
        int cacheSize=10000                          ///< Kernel cache size
    );

private:

    afw::table::ExposureCatalog _catalog;
    CONST_PTR(afw::image::Wcs) _coaddWcs;
    afw::table::Key<double> _weightKey;
    afw::geom::Point2D _averagePosition;
    std::string _warpingKernelName;   // could be removed if we could get this from _warpingControl (#2949)
    CONST_PTR(afw::math::WarpingControl) _warpingControl;
};

}}} // namespace lsst::meas::algorithms

#endif
