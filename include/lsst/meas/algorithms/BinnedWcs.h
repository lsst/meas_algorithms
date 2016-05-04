// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2014 LSST Corporation.
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

#ifndef LSST_MEAS_ALGORITHMS_BINNEDWCS_H
#define LSST_MEAS_ALGORITHMS_BINNEDWCS_H

#include <memory>

#include "lsst/pex/exceptions.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Point.h"


namespace lsst { namespace meas { namespace algorithms {

/// An exception that indicates the feature has not been implemented
LSST_EXCEPTION_TYPE(NotImplementedException, lsst::pex::exceptions::RuntimeError,
                    lsst::meas::algorithms::NotImplementedException);

class BinnedWcs : public afw::image::Wcs, public std::enable_shared_from_this<BinnedWcs> {
public:
    BinnedWcs(PTR(afw::image::Wcs) parent, unsigned int xBin, unsigned int yBin, afw::geom::Point2I xy0);
    virtual ~BinnedWcs() {}

    virtual PTR(afw::image::Wcs) clone() const {
        return PTR(afw::image::Wcs)(new BinnedWcs(_parent, _xBin, _yBin, _xy0));
    }
    PTR(afw::image::Wcs) upcast() { return shared_from_this(); }

    PTR(afw::image::Wcs) getParent() const { return _parent; }
    unsigned int getXBin() const { return _xBin; }
    unsigned int getYBin() const { return _yBin; }
    afw::geom::Point2I getXY0() const { return _xy0; }

    virtual bool hasDistortion() const { return _parent->hasDistortion(); }
    virtual bool isPersistable() const { return false; }

    virtual void flipImage(int flipLR, int flipTB, afw::geom::Extent2I dimensions) const {
        notImplemented();
    }
    virtual void rotateImageBy90(int nQuarter, afw::geom::Extent2I dimensions) const {
        notImplemented();
    }
    virtual PTR(daf::base::PropertyList) getFitsMetadata() const {
        notImplemented();
        return PTR(daf::base::PropertyList)(); // unreached
    }

    afw::geom::AffineTransform getBinnedToOriginal() const { return _binnedToOriginal; }
    afw::geom::AffineTransform getOriginalToBinned() const { return _originalToBinned; }

protected:
    virtual void pixelToSkyImpl(double pixel1, double pixel2, afw::geom::Angle skyTmp[2]) const;
    virtual afw::geom::Point2D skyToPixelImpl(afw::geom::Angle sky1, afw::geom::Angle sky2) const;

private:
    void notImplemented() const {
        throw LSST_EXCEPT(NotImplementedException, "Feature has not been implemented");
    }
    PTR(afw::image::Wcs) const _parent;
    unsigned int const _xBin, _yBin;
    afw::geom::Point2I const _xy0;
    afw::geom::AffineTransform const _binnedToOriginal, _originalToBinned;
};

}}} // namespace lsst::meas::algorithms

#endif  // LSST_MEAS_ALGORITHMS_BINNED_WCS_H
