#include "lsst/meas/algorithms/BinnedWcs.h"
#include "wcslib/wcs.h" // Unfortunate, but necessary due to parent class mixing implementation and interface

namespace lsst { namespace meas { namespace algorithms {

BinnedWcs::BinnedWcs(
    PTR(afw::image::Wcs) parent,
    unsigned int xBin,
    unsigned int yBin,
    afw::geom::Point2I xy0
    ) :
    afw::image::Wcs(*parent), // We have to inherit the implementation as well as the interface...
    _parent(parent), _xBin(xBin), _yBin(yBin), _xy0(xy0),
    _binnedToOriginal(afw::geom::AffineTransform::makeTranslation(afw::geom::Extent2D(_xy0))*
                      afw::geom::AffineTransform::makeScaling(_xBin, _yBin)),
    _originalToBinned(_binnedToOriginal.invert())
{
    // Correct CRPIX
    afw::geom::Point2D const crpix = _originalToBinned(parent->getPixelOrigin());
    _wcsInfo->crpix[0] = crpix.getX() + 1; // convert LSST --> FITS
    _wcsInfo->crpix[1] = crpix.getY() + 1; // convert LSST --> FITS
}


void BinnedWcs::pixelToSkyImpl(double x, double y, afw::geom::Angle skyTmp[2]) const
{
    PTR(afw::coord::Coord) coord = _parent->pixelToSky(_binnedToOriginal(afw::geom::Point2D(x, y)));
    skyTmp[0] = coord->getLongitude();
    skyTmp[1] = coord->getLatitude();
}

afw::geom::Point2D BinnedWcs::skyToPixelImpl(afw::geom::Angle sky1, afw::geom::Angle sky2) const
{
    return _originalToBinned(_parent->skyToPixel(*makeCorrectCoord(sky1, sky2)));
}

}}} // namespace lsst::meas::algorithms
