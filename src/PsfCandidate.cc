// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2008-2015 AURA/LSST.
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
 * @brief Implementation of code to determine spatial model of PSF
 *
 * @file
 *
 * @ingroup algorithms
 */
#include "lsst/afw/detection/Footprint.h"
#include "lsst/geom.h"
#include "lsst/afw/image/ImageAlgorithm.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/meas/algorithms/PsfCandidate.h"

namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/*
 * PsfCandidate's members
 */
template <typename PixelT>
int PsfCandidate<PixelT>::_border = 0;
template <typename PixelT>
int PsfCandidate<PixelT>::_defaultWidth = 21;
template <typename PixelT>
float PsfCandidate<PixelT>::_pixelThreshold = 0.0;
template <typename PixelT>
bool PsfCandidate<PixelT>::_doMaskBlends = true;

/************************************************************************************************************/
namespace {
template <typename T>  // functor used by makeImageFromMask to return inputMask
struct noop : public afw::image::pixelOp1<T> {
    T operator()(T x) const { return x; }
};

template <typename T>  // functor used by makeImageFromMask to return (inputMask & mask)
struct andMask : public afw::image::pixelOp1<T> {
    andMask(T mask) : _mask(mask) {}
    T operator()(T x) const { return (x & _mask); }

private:
    T _mask;
};

template <typename T>
andMask<T> makeAndMask(T val) {
    return andMask<T>(val);
}

/*
 * Return an Image initialized from a Mask (possibly modified by func)
 */
template <typename LhsT, typename RhsT>
std::shared_ptr<afw::image::Image<LhsT>> makeImageFromMask(
        afw::image::Mask<RhsT> const& rhs,                     ///< mask to process
        afw::image::pixelOp1<RhsT> const& func = noop<RhsT>()  ///< functor to call
) {
    std::shared_ptr<afw::image::Image<LhsT>> lhs =
            std::make_shared<afw::image::Image<LhsT>>(rhs.getDimensions());
    lhs->setXY0(rhs.getXY0());

    for (int y = 0; y != lhs->getHeight(); ++y) {
        typename afw::image::Image<RhsT>::const_x_iterator rhsPtr = rhs.row_begin(y);

        for (typename afw::image::Image<LhsT>::x_iterator lhsPtr = lhs->row_begin(y),
                                                          lhsEnd = lhs->row_end(y);
             lhsPtr != lhsEnd; ++rhsPtr, ++lhsPtr) {
            *lhsPtr = func(*rhsPtr);
        }
    }

    return lhs;
}

/// Return square of the distance between a point and a peak
double distanceSquared(double x, double y, afw::detection::PeakRecord const& peak) {
    return std::pow(peak.getIx() - x, 2) + std::pow(peak.getIy() - y, 2);
}

/// Functor to mask pixels in a blended candidate
///
/// We mask pixels in the footprint that are closest to a peak other than the central one.
/// For these pixels, we activate the "turnOn" bit mask and deactivate the "turnOff" bit mask.
template <typename MaskT>
class BlendedFunctor {
public:
    BlendedFunctor(afw::detection::PeakRecord const& central,  ///< Central peak
                   afw::detection::PeakCatalog const& peaks,   ///< Other peaks
                   afw::image::MaskPixel turnOff,              ///< Bit mask to deactivate
                   afw::image::MaskPixel turnOn                ///< Bit mask to activate
                   )
            : _central(central), _peaks(peaks), _turnOff(~turnOff), _turnOn(turnOn) {}

    /// Functor operation to mask pixels closest to a peak other than the central one.
    void operator()(geom::Point2I const& point, MaskT& val) {
        int x = point.getX();
        int y = point.getY();
        double const central = distanceSquared(x, y, _central);
        for (afw::detection::PeakCatalog::const_iterator iter = _peaks.begin(), end = _peaks.end();
             iter != end; ++iter) {
            double const dist2 = distanceSquared(x, y, *iter);
            if (dist2 < central) {
                val &= _turnOff;
                val |= _turnOn;
            }
        }
    }

private:
    afw::detection::PeakRecord const& _central;
    afw::detection::PeakCatalog const& _peaks;
    afw::image::MaskPixel const _turnOff;
    afw::image::MaskPixel const _turnOn;
};

}  // anonymous namespace

/// Extract an image of the candidate.
///
/// The MaskedImage is a deep copy of a sub-image of the original image.  No offsets are applied.
///
/// In the mask, the INTRP bit is set and DETECTED unset for any pixels that are not considered part of the
/// actual candidate.  You should consider that, for the output mask:
/// * INTRP means "ignore this pixel", i.e., it's contaminated
/// * DETECTED means "fit this pixel", i.e., it contains the object of interest
/// * Nothing means "do what you want"
///
/// Three schemes are used for masking pixels:
/// * Pixels closer to a peak in the source footprint other than the central peak are masked.  This deals with
///   sources blended with the actual candidate.
/// * Sources (identified from the DETECTED bit plane) unconnected to the central peak are masked, and this
///   mask is grown.  This deals with bright neighbouring sources.
/// * Pixels exceeding the pixelThreshold relative to the expected noise according to the variance plane are
///   masked if they are not in the central footprint.  This is only applied if the pixelThreshold is
///   positive.  This deals with faint neighbouring sources.
template <typename PixelT>
PTR(afw::image::MaskedImage<PixelT>)
PsfCandidate<PixelT>::extractImage(unsigned int width,  // Width of image
                                   unsigned int height  // Height of image
                                   ) const {
    geom::Point2I const cen(afw::image::positionToIndex(getXCenter()),
                            afw::image::positionToIndex(getYCenter()));
    geom::Point2I const llc(cen[0] - width / 2 - _parentExposure->getX0(),
                            cen[1] - height / 2 - _parentExposure->getY0());

    geom::BoxI bbox(llc, geom::ExtentI(width, height));

    PTR(MaskedImageT) image;
    try {
        MaskedImageT mimg = _parentExposure->getMaskedImage();
        image.reset(new MaskedImageT(mimg, bbox, afw::image::LOCAL, true));  // a deep copy
    } catch (pex::exceptions::LengthError& e) {
        LSST_EXCEPT_ADD(e, "Extracting image of PSF candidate");
        throw e;
    }

    //
    // Set INTRP and unset DETECTED for any pixels we don't want to deal with.
    //
    afw::image::MaskPixel const intrp =
            MaskedImageT::Mask::getPlaneBitMask("INTRP");  // mask bit for bad pixels
    afw::image::MaskPixel const detected = MaskedImageT::Mask::getPlaneBitMask("DETECTED");  // object pixels

    // Mask out blended objects
    if (getMaskBlends()) {
        CONST_PTR(afw::detection::Footprint) foot = getSource()->getFootprint();
        typedef afw::detection::PeakCatalog PeakCatalog;
        PeakCatalog const& peaks = foot->getPeaks();
        if (peaks.size() > 1) {
            // Mask all pixels in the footprint except for those closest to the central peak
            double best = std::numeric_limits<double>::infinity();
            PTR(afw::detection::PeakRecord) central;
            for (PeakCatalog::const_iterator iter = peaks.begin(), end = peaks.end(); iter != end; ++iter) {
                double const dist2 = distanceSquared(getXCenter(), getYCenter(), *iter);
                if (dist2 < best) {
                    best = dist2;
                    central = iter;
                }
            }
            assert(central);  // We must have found something

            PeakCatalog others(peaks.getTable());
            others.reserve(peaks.size() - 1);
            for (PeakCatalog::const_iterator iter = peaks.begin(), end = peaks.end(); iter != end; ++iter) {
                PTR(afw::detection::PeakRecord) ptr(iter);
                if (central != ptr) {
                    others.push_back(ptr);
                }
            }

            BlendedFunctor<typename MaskedImageT::Mask::Pixel> functor(*central, others, detected, intrp);
            foot->getSpans()->clippedTo(image->getBBox())->applyFunctor(functor, *image->getMask());
        }
    }

    /*
     * Mask any DETECTED pixels other than the one in the center of the object;
     * we grow the Footprint a bit first
     */
    typedef afw::detection::FootprintSet::FootprintList FootprintList;

    PTR(afw::image::Image<int>) mim = makeImageFromMask<int>(*image->getMask(), makeAndMask(detected));
    PTR(afw::detection::FootprintSet)
    fs = std::make_shared<afw::detection::FootprintSet>(*mim, afw::detection::Threshold(1));
    CONST_PTR(FootprintList) feet = fs->getFootprints();

    if (feet->size() > 1) {
        int const ngrow = 3;  // number of pixels to grow bad Footprints
        //
        // Go through Footprints looking for ones that don't contain cen
        //
        for (FootprintList::const_iterator fiter = feet->begin(); fiter != feet->end(); ++fiter) {
            PTR(afw::detection::Footprint) foot = *fiter;
            if (foot->contains(cen)) {
                continue;
            }

            // Dilate and clip to the image bounding box, incase the span grows outside the image
            auto bigSpan = foot->getSpans()->dilated(ngrow)->clippedTo(image->getBBox());
            bigSpan->clearMask(*image->getMask(), detected);
            bigSpan->setMask(*image->getMask(), intrp);
        }
    }

    // Mask high pixels unconnected to the center
    if (_pixelThreshold > 0.0) {
        CONST_PTR(afw::detection::FootprintSet)
        fpSet = std::make_shared<afw::detection::FootprintSet>(
                *image, afw::detection::Threshold(_pixelThreshold, afw::detection::Threshold::PIXEL_STDEV));
        for (FootprintList::const_iterator fpIter = fpSet->getFootprints()->begin();
             fpIter != fpSet->getFootprints()->end(); ++fpIter) {
            CONST_PTR(afw::detection::Footprint) fp = *fpIter;
            if (!fp->contains(cen)) {
                fp->getSpans()->clearMask(*image->getMask(), detected);
                fp->getSpans()->setMask(*image->getMask(), intrp);
            }
        }
    }

    return image;
}

/**
 * Return the %image at the position of the Source, without any sub-pixel shifts to put the centre of the
 * object in the centre of a pixel (for that, use getOffsetImage())
 *
 */
template <typename PixelT>
CONST_PTR(afw::image::MaskedImage<PixelT>)
PsfCandidate<PixelT>::getMaskedImage(int width, int height) const {
    if (!_image || (width != _image->getWidth() || height != _image->getHeight())) {
        _image = extractImage(width, height);
    }
    return _image;
}

/**
 * Return the %image at the position of the Source, without any sub-pixel shifts to put the centre of the
 * object in the centre of a pixel (for that, use getOffsetImage())
 *
 */
template <typename PixelT>
CONST_PTR(afw::image::MaskedImage<PixelT>)
PsfCandidate<PixelT>::getMaskedImage() const {
    int const width = getWidth() == 0 ? _defaultWidth : getWidth();
    int const height = getHeight() == 0 ? _defaultWidth : getHeight();
    return getMaskedImage(width, height);
}

template <typename PixelT>
int PsfCandidate<PixelT>::getBorderWidth() {
    return _border;
}

template <typename PixelT>
void PsfCandidate<PixelT>::setBorderWidth(int border) {
    _border = border;
}

template <typename PixelT>
void PsfCandidate<PixelT>::setPixelThreshold(float threshold) {
    _pixelThreshold = threshold;
}

template <typename PixelT>
float PsfCandidate<PixelT>::getPixelThreshold() {
    return _pixelThreshold;
}

template <typename PixelT>
void PsfCandidate<PixelT>::setMaskBlends(bool doMaskBlends) {
    _doMaskBlends = doMaskBlends;
}

template <typename PixelT>
bool PsfCandidate<PixelT>::getMaskBlends() {
    return _doMaskBlends;
}

/**
 * @brief Return an offset version of the image of the source.
 * The returned image has been offset to put the centre of the object in the centre of a pixel.
 *
 */
template <typename PixelT>
PTR(afw::image::MaskedImage<PixelT>)
PsfCandidate<PixelT>::getOffsetImage(std::string const algorithm,  // Warping algorithm to use
                                     unsigned int buffer           // Buffer for warping
                                     ) const {
    unsigned int const width = getWidth() == 0 ? _defaultWidth : getWidth();
    unsigned int const height = getHeight() == 0 ? _defaultWidth : getHeight();
    if (_offsetImage && static_cast<unsigned int>(_offsetImage->getWidth()) == width + 2 * buffer &&
        static_cast<unsigned int>(_offsetImage->getHeight()) == height + 2 * buffer) {
        return _offsetImage;
    }

    PTR(MaskedImageT) image = extractImage(width + 2 * buffer, height + 2 * buffer);

    double const xcen = getXCenter(), ycen = getYCenter();
    double const dx = afw::image::positionToIndex(xcen, true).second;
    double const dy = afw::image::positionToIndex(ycen, true).second;

    PTR(MaskedImageT) offset = afw::math::offsetImage(*image, -dx, -dy, algorithm);
    geom::Point2I llc(buffer, buffer);
    geom::Extent2I dims(width, height);
    geom::Box2I box(llc, dims);
    _offsetImage.reset(new MaskedImageT(*offset, box, afw::image::LOCAL, true));  // Deep copy

    return _offsetImage;
}

/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
typedef float Pixel;
// template class PsfCandidate<afw::image::MaskedImage<Pixel> >;
template class PsfCandidate<Pixel>;
/// \endcond

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
