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

#include "lsst/geom/AffineTransform.h"
#include "lsst/geom/Box.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/WarpedPsf.h"
#include "lsst/afw/math/warpExposure.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/Persistable.cc"

namespace lsst {
namespace afw {
namespace table {
namespace io {

template std::shared_ptr<meas::algorithms::WarpedPsf>
PersistableFacade<meas::algorithms::WarpedPsf>::dynamicCast(std::shared_ptr<Persistable> const &);

}  // namespace io
}  // namespace table
}  // namespace afw

namespace meas {
namespace algorithms {

namespace {

// TODO: make this routine externally callable and more generic using templates
//  (also useful in e.g. math/offsetImage.cc)
std::shared_ptr<afw::detection::Psf::Image> zeroPadImage(afw::detection::Psf::Image const &im, int xPad,
                                                         int yPad) {
    int nx = im.getWidth();
    int ny = im.getHeight();

    auto out = std::make_shared<afw::detection::Psf::Image>(nx + 2 * xPad, ny + 2 * yPad);
    out->setXY0(im.getX0() - xPad, im.getY0() - yPad);

    geom::Box2I box(geom::Point2I(xPad, yPad), geom::Extent2I(nx, ny));
    out->assign(im, box, afw::image::LOCAL);

    return out;
}

geom::Box2I computeBBoxFromTransform(geom::Box2I const bbox, afw::geom::TransformPoint2ToPoint2 const &t, geom::Point2D const &dest_position) {
    // This is the maximum scaling I can imagine we'd want - it's approximately what you'd
    // get from trying to coadd 2"-pixel data (e.g. 2MASS) onto a 0.01"-pixel grid (e.g.
    // from JWST).  Anything beyond that is probably a bug elsewhere in the code, and
    // catching that can save us from segfaults and catastrophic memory consumption.
    // static const double maxTransformCoeff = 200.0;
    // if (t.getLinear().getMatrix().lpNorm<Eigen::Infinity>() > maxTransformCoeff) {
    //     throw LSST_EXCEPT(pex::exceptions::RangeError, "Unexpectedly large transform passed to WarpedPsf");
    // }

    // Floating point version of input bbox (expanded to include entire pixels).
    geom::Box2D in_box_fp(bbox);
    // Floating point version of output bbox (to be populated as bbox of
    // transformed points.
    geom::Box2D out_box_fp;
    auto const in_corners = in_box_fp.getCorners();
    for (auto const & in_corner : in_corners) {
        auto out_corner = t.applyForward(in_corner);
        out_box_fp.include(out_corner);
    }

    // We want to guarantee that the output bbox has odd dimensions, so instead
    // of using the Box2I converting constructor directly, we start with (0, 0)
    // and dilate by the floating-point box's half-dimensions.
    geom::Extent2I out_half_dims = geom::floor(0.5*out_box_fp.getDimensions());
    geom::Box2I out_box;
    geom::Point2I out_center(dest_position);
    out_box.include(out_center);
    return out_box.dilatedBy(out_half_dims);
}

/**
 * @brief Alternate interface to afw::math::warpImage()
 * in which the caller does not need to precompute the output bounding box.
 *
 * This version takes an affine transform instead of an arbitrary xy transform.
 *
 * @param[in] im  Image to warp
 * @param[in] srcToDest  Affine transformation from source pixels to destination pixels in the forward
 *                  direction; the warping code only uses the inverse direction
 * @param[in] wc  Warping parameters
 *
 * The input image is assumed zero-padded.
 */


// A helper struct for logic and intermediate results used by both
// doComputeKernelImage and doComputeBBox.
//
// This linearizes the transform from destination to source coordinates, then
// computes both a source position that corresponds to a rounded version of the
// destination position and a transform from source to destination that also
// undoes the rounding.

}  // namespace

WarpedPsf::WarpedPsf(std::shared_ptr<afw::detection::Psf const> undistortedPsf,
                     std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const> distortion,
                     std::shared_ptr<afw::math::WarpingControl const> control)
        : ImagePsf(false),
          _undistortedPsf(undistortedPsf),
          _distortion(distortion),
          _warpingControl(control) {
    _init();
}

WarpedPsf::WarpedPsf(std::shared_ptr<afw::detection::Psf const> undistortedPsf,
                     std::shared_ptr<afw::geom::TransformPoint2ToPoint2 const> distortion,
                     std::string const &kernelName, unsigned int cache)
        : ImagePsf(false),
          _undistortedPsf(undistortedPsf),
          _distortion(distortion),
          _warpingControl(new afw::math::WarpingControl(kernelName, "", cache)) {
    _init();
}

void WarpedPsf::_init() {
    if (!_undistortedPsf) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "Undistorted Psf passed to WarpedPsf must not be None/NULL");
    }
    if (!_distortion) {
        throw LSST_EXCEPT(pex::exceptions::LogicError, "Transform passed to WarpedPsf must not be None/NULL");
    }
    if (!_warpingControl) {
        throw LSST_EXCEPT(pex::exceptions::LogicError,
                          "WarpingControl passed to WarpedPsf must not be None/NULL");
    }
}

geom::Point2D WarpedPsf::getAveragePosition() const {
    return _distortion->applyForward(_undistortedPsf->getAveragePosition());
}

std::shared_ptr<afw::detection::Psf> WarpedPsf::clone() const {
    return std::make_shared<WarpedPsf>(_undistortedPsf->clone(), _distortion, _warpingControl);
}

std::shared_ptr<afw::detection::Psf> WarpedPsf::resized(int width, int height) const {
    // For a given set of requested dimensions and distortion, it is not guaranteed that a
    // _undistortedPsf would exist to manifest those dimensions after distortion
    // Not possible to implement with member data currently in WarpedPsf
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

std::shared_ptr<afw::detection::Psf::Image> WarpedPsf::doComputeKernelImage(
        geom::Point2D const &position, afw::image::Color const &color) const {
    geom::Point2D rounded_dest_position(
            std::round(position.getX()),
            std::round(position.getY())
        );
    auto src_position = _distortion->applyInverse(rounded_dest_position);
    std::shared_ptr<Image> src_im = _undistortedPsf->computeImage(src_position, color);

    afw::math::SeparableKernel const &kernel = *_warpingControl->getWarpingKernel();
    geom::Point2I const &center = kernel.getCtr();
    int const xPad = std::max(center.getX(), kernel.getWidth() - center.getX());
    int const yPad = std::max(center.getY(), kernel.getHeight() - center.getY());
    auto padded_src_im = zeroPadImage(*src_im, xPad, yPad);

    auto ret = std::make_shared<afw::detection::Psf::Image>(computeBBoxFromTransform(padded_src_im->getBBox(), *_distortion, rounded_dest_position));

    afw::math::warpImage(
        *ret,
        *padded_src_im,
        *_distortion,
        *_warpingControl,
        0.0
    );
    ret->setXY0(-ret->getBBox().getDimensions().getX() / 2, -ret->getBBox().getDimensions().getY() / 2);

    double normFactor = 1.0;
    //
    // Normalize the output image to sum 1
    // FIXME defining a member function Image::getSum() would be convenient here and in other places
    //
    normFactor = 0.0;
    for (int y = 0; y != ret->getHeight(); ++y) {
        afw::detection::Psf::Image::x_iterator imEnd = ret->row_end(y);
        for (afw::detection::Psf::Image::x_iterator imPtr = ret->row_begin(y); imPtr != imEnd; imPtr++) {
            normFactor += *imPtr;
        }
    }
    if (normFactor == 0.0) {
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError, "psf image has sum 0");
    }
    *ret /= normFactor;
    return ret;
}

geom::Box2I WarpedPsf::doComputeBBox(geom::Point2D const &position, afw::image::Color const &color) const {
    throw LSST_EXCEPT(pex::exceptions::LogicError, "Not Implemented");
}

namespace {

struct PersistenceHelper {
    afw::table::Schema schema;
    afw::table::Key<int> psfIndex;
    afw::table::Key<int> transformIndex;
    afw::table::Key<int> controlIndex;

    static PersistenceHelper const &get() {
        static PersistenceHelper const instance;
        return instance;
    }

private:
    PersistenceHelper()
            : schema(),
              psfIndex(schema.addField<int>("psfIndex", "archive ID of nested Psf object")),
              transformIndex(schema.addField<int>("transformIndex", "archive ID of nested Transform object")),
              controlIndex(
                      schema.addField<int>("controlIndex", "archive ID of nested WarpingControl object")) {}
};

std::string _getPersistenceName() { return "WarpedPsf"; }

class : public afw::table::io::PersistableFactory {
public:
    virtual std::shared_ptr<afw::table::io::Persistable> read(
            afw::table::io::InputArchive const &archive,
            afw::table::io::CatalogVector const &catalogs) const {
        static PersistenceHelper const &keys = PersistenceHelper::get();
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        LSST_ARCHIVE_ASSERT(catalogs.front().size() == 1u);
        afw::table::BaseRecord const &record = catalogs.front().front();
        LSST_ARCHIVE_ASSERT(record.getSchema() == keys.schema);
        return std::make_shared<WarpedPsf>(
                archive.get<afw::detection::Psf>(record.get(keys.psfIndex)),
                archive.get<afw::geom::TransformPoint2ToPoint2>(record.get(keys.transformIndex)),
                archive.get<afw::math::WarpingControl>(record.get(keys.controlIndex)));
    }

    using afw::table::io::PersistableFactory::PersistableFactory;
} warpedPsfFactory(_getPersistenceName());

}  // namespace

std::string WarpedPsf::getPersistenceName() const { return _getPersistenceName(); }
std::string WarpedPsf::getPythonModule() const { return "lsst.meas.algorithms"; }

void WarpedPsf::write(OutputArchiveHandle &handle) const {
    static PersistenceHelper const &keys = PersistenceHelper::get();
    afw::table::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    std::shared_ptr<afw::table::BaseRecord> record = catalog.addNew();

    record->set(keys.psfIndex, handle.put(_undistortedPsf));
    record->set(keys.transformIndex, handle.put(_distortion));
    record->set(keys.controlIndex, handle.put(_warpingControl));

    handle.saveCatalog(catalog);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
