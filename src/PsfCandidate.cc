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
 * @brief Implementation of code to determine spatial model of PSF
 *
 * @file
 *
 * @ingroup algorithms
 */
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/FootprintFunctor.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/ImageAlgorithm.h"
#include "lsst/afw/math/offsetImage.h"
#include "lsst/meas/algorithms/PsfCandidate.h"

namespace afwDetection = lsst::afw::detection;
namespace afwGeom      = lsst::afw::geom;
namespace afwImage     = lsst::afw::image;
namespace afwMath      = lsst::afw::math;
namespace measAlg      = lsst::meas::algorithms;

/************************************************************************************************************/
/*
 * PsfCandidate's members
 */
template <typename PixelT>
int measAlg::PsfCandidate<PixelT>::_border = 0;
template <typename PixelT>
int measAlg::PsfCandidate<PixelT>::_defaultWidth = 21;
template <typename PixelT>
float measAlg::PsfCandidate<PixelT>::_pixelThreshold = 0.0;
template <typename PixelT>
bool measAlg::PsfCandidate<PixelT>::_doMaskBlends = true;

/************************************************************************************************************/
namespace {
    template<typename T>                // functor used by makeImageFromMask to return inputMask
    struct noop : public afwImage::pixelOp1<T> {
        T operator()(T x) const { return x; }
    };

    template<typename T>                // functor used by makeImageFromMask to return (inputMask & mask)
    struct andMask : public afwImage::pixelOp1<T> {
        andMask(T mask) : _mask(mask) {}
        T operator()(T x) const { return (x & _mask); }
    private:
        T _mask;
    };

    template<typename T>
    andMask<T> makeAndMask(T val) {
        return andMask<T>(val);
    }

    /*
     * Return an Image initialized from a Mask (possibly modified by func)
     */
    template<typename LhsT, typename RhsT>
    typename afwImage::Image<LhsT>::Ptr
    makeImageFromMask(afwImage::Mask<RhsT> const& rhs,     ///< mask to process
                      afwImage::pixelOp1<RhsT> const& func=noop<RhsT>() ///< functor to call
                     )
    {
        typename afwImage::Image<LhsT>::Ptr lhs =
            boost::make_shared<afwImage::Image<LhsT> >(rhs.getDimensions());
        lhs->setXY0(rhs.getXY0());

        for (int y = 0; y != lhs->getHeight(); ++y) {
            typename afwImage::Image<RhsT>::const_x_iterator rhsPtr = rhs.row_begin(y);

            for (typename afwImage::Image<LhsT>::x_iterator lhsPtr = lhs->row_begin(y),
                     lhsEnd = lhs->row_end(y); lhsPtr != lhsEnd; ++rhsPtr, ++lhsPtr) {
                *lhsPtr = func(*rhsPtr);
            }
        }
        
        return lhs;
    }

    /// Return square of the distance between a point and a peak
    double distanceSquared(double x, double y, afwDetection::Peak const& peak) {
        return std::pow(peak.getIx() - x, 2) + std::pow(peak.getIy() - y, 2);
    }

    /// Functor to mask pixels in a blended candidate
    ///
    /// We mask pixels in the footprint that are closest to a peak other than the central one.
    /// For these pixels, we activate the "turnOn" bit mask and deactivate the "turnOff" bit mask.
    template <typename PixelT>
    class BlendedFunctor : public afwDetection::FootprintFunctor<afwImage::MaskedImage<PixelT> > {
    public:
        typedef typename afwImage::MaskedImage<PixelT> Image;
        typedef typename Image::Mask Mask;
        typedef typename afwDetection::FootprintFunctor<Image> Super;
        BlendedFunctor(
            Image const& image,             ///< Image; unused except for xy0, but required by superclass
            Mask & mask,                    ///< Mask to modify
            afwDetection::Peak const& central, ///< Central peak
            afwDetection::Footprint::PeakList const& peaks, ///< Other peaks
            afwImage::MaskPixel turnOff,                   ///< Bit mask to deactivate
            afwImage::MaskPixel turnOn                     ///< Bit mask to activate
            ) :
            Super(image),
            _central(central),
            _peaks(peaks),
            _mask(mask),
            _turnOff(~turnOff),
            _turnOn(turnOn)
            {}

        /// Functor operation to mask pixels closest to a peak other than the central one.
        virtual void operator()(typename Image::xy_locator loc, int x, int y) {
            double const central = distanceSquared(x, y, _central);
            int const xImage = x - getImage().getX0();
            int const yImage = y - getImage().getY0();
            for (afwDetection::Footprint::PeakList::const_iterator iter = _peaks.begin(), end = _peaks.end();
                 iter != end; ++iter) {
                double const dist2 = distanceSquared(x, y, **iter);
                if (dist2 < central) {
                    (_mask)(xImage, yImage) &= _turnOff;
                    (_mask)(xImage, yImage) |= _turnOn;
                    return;
                }
            }
        }

        using Super::getImage;

    private:
        afwDetection::Peak const& _central;
        afwDetection::Footprint::PeakList const& _peaks;
        Mask & _mask;
        afwImage::MaskPixel const _turnOff;
        afwImage::MaskPixel const _turnOn;
    };

} // anonymous namespace

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
PTR(afwImage::MaskedImage<PixelT>)
measAlg::PsfCandidate<PixelT>::extractImage(
    unsigned int width,                 // Width of image
    unsigned int height                 // Height of image
) const {
    afwGeom::Point2I const cen(afwImage::positionToIndex(getXCenter()),
                               afwImage::positionToIndex(getYCenter()));
    afwGeom::Point2I const llc(cen[0] - width/2 - _parentExposure->getX0(), 
                               cen[1] - height/2 - _parentExposure->getY0());
    
    afwGeom::BoxI bbox(llc, afwGeom::ExtentI(width, height));
        
    PTR(MaskedImageT) image;
    try {
        MaskedImageT mimg = _parentExposure->getMaskedImage();
        image.reset(new MaskedImageT(mimg, bbox, afwImage::LOCAL, true)); // a deep copy
    } catch(lsst::pex::exceptions::LengthErrorException &e) {
        LSST_EXCEPT_ADD(e, "Extracting image of PSF candidate");
        throw e;
    }

    //
    // Set INTRP and unset DETECTED for any pixels we don't want to deal with.
    //
    afwImage::MaskPixel const intrp = MaskedImageT::Mask::getPlaneBitMask("INTRP"); // mask bit for bad pixels
    afwImage::MaskPixel const detected = MaskedImageT::Mask::getPlaneBitMask("DETECTED"); // object pixels

    // Mask out blended objects
    if (getMaskBlends()) {
        CONST_PTR(afwDetection::Footprint) foot = getSource()->getFootprint();
        typedef afwDetection::Footprint::PeakList PeakList;
        PeakList const& peaks = foot->getPeaks();
        if (peaks.size() > 1) {
            // Mask all pixels in the footprint except for those closest to the central peak
            double best = std::numeric_limits<double>::infinity();
            CONST_PTR(afwDetection::Peak) central;
            for (PeakList::const_iterator iter = peaks.begin(), end = peaks.end(); iter != end; ++iter) {
                double const dist2 = distanceSquared(getXCenter(), getYCenter(), **iter);
                if (dist2 < best) {
                    best = dist2;
                    central = *iter;
                }
            }
            assert(central);                // We must have found something

            PeakList others;
            others.reserve(peaks.size() - 1);
            for (PeakList::const_iterator iter = peaks.begin(), end = peaks.end(); iter != end; ++iter) {
                if (central != *iter) {
                    others.push_back(*iter);
                }
            }

            BlendedFunctor<PixelT> functor(*image, *image->getMask(), *central, others, detected, intrp);
            functor.apply(*foot);
        }
    }

    /*
     * Mask any DETECTED pixels other than the one in the center of the object;
     * we grow the Footprint a bit first
     */
    typedef afwDetection::FootprintSet::FootprintList FootprintList;

    PTR(afwImage::Image<int>) mim = makeImageFromMask<int>(*image->getMask(), makeAndMask(detected));
    PTR(afwDetection::FootprintSet) fs =
        boost::make_shared<afwDetection::FootprintSet>(*mim, afwDetection::Threshold(1));
    CONST_PTR(FootprintList) feet = fs->getFootprints();

    if (feet->size() > 1) {
        int const ngrow = 3;            // number of pixels to grow bad Footprints
        //
        // Go through Footprints looking for ones that don't contain cen
        //
        for (FootprintList::const_iterator fiter = feet->begin(); fiter != feet->end(); ++fiter) {
            PTR(afwDetection::Footprint) foot = *fiter;
            if (foot->contains(cen)) {
                continue;
            }
            
            PTR(afwDetection::Footprint) bigfoot = afwDetection::growFootprint(foot, ngrow);
            afwDetection::clearMaskFromFootprint(image->getMask().get(), *bigfoot, detected);
            afwDetection::setMaskFromFootprint(image->getMask().get(), *bigfoot, intrp);
        }
    }

    // Mask high pixels unconnected to the center
    if (_pixelThreshold > 0.0) {
        CONST_PTR(afwDetection::FootprintSet) fpSet =
            boost::make_shared<afwDetection::FootprintSet>(*image,
                afwDetection::Threshold(_pixelThreshold, afwDetection::Threshold::PIXEL_STDEV));
        for (FootprintList::const_iterator fpIter = fpSet->getFootprints()->begin();
             fpIter != fpSet->getFootprints()->end(); ++fpIter) {
            CONST_PTR(afwDetection::Footprint) fp = *fpIter;
            if (!fp->contains(cen)) {
                afwDetection::clearMaskFromFootprint(image->getMask().get(), *fp, detected);
                afwDetection::setMaskFromFootprint(image->getMask().get(), *fp, intrp);
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
CONST_PTR(afwImage::MaskedImage<PixelT>)
measAlg::PsfCandidate<PixelT>::getMaskedImage(int width, int height) const {


    if (_haveImage && (width != _image->getWidth() || height != _image->getHeight())) {
        _haveImage = false;
    }

    if (!_haveImage) {
        _image = extractImage(width, height);
        _haveImage = true;
    }
    
    return _image;
}

/**
 * Return the %image at the position of the Source, without any sub-pixel shifts to put the centre of the
 * object in the centre of a pixel (for that, use getOffsetImage())
 *
 */
template <typename PixelT>
CONST_PTR(afwImage::MaskedImage<PixelT>) measAlg::PsfCandidate<PixelT>::getMaskedImage() const {

    int const width = getWidth() == 0 ? _defaultWidth : getWidth();
    int const height = getHeight() == 0 ? _defaultWidth : getHeight();

    return getMaskedImage(width, height);
    
}

/**
 * @brief Return an offset version of the image of the source.
 * The returned image has been offset to put the centre of the object in the centre of a pixel.
 *
 */
template <typename PixelT>
PTR(afwImage::MaskedImage<PixelT>)
measAlg::PsfCandidate<PixelT>::getOffsetImage(
    std::string const algorithm,        // Warping algorithm to use
    unsigned int buffer                 // Buffer for warping
) const {
    unsigned int const width = getWidth() == 0 ? _defaultWidth : getWidth();
    unsigned int const height = getHeight() == 0 ? _defaultWidth : getHeight();
    if (_offsetImage && static_cast<unsigned int>(_offsetImage->getWidth()) == width + 2*buffer && 
        static_cast<unsigned int>(_offsetImage->getHeight()) == height + 2*buffer) {
        return _offsetImage;
    }

    PTR(MaskedImageT) image = extractImage(width + 2*buffer, height + 2*buffer);

    double const xcen = getXCenter(), ycen = getYCenter();
    double const dx = afwImage::positionToIndex(xcen, true).second;
    double const dy = afwImage::positionToIndex(ycen, true).second;

    PTR(MaskedImageT) offset = afwMath::offsetImage(*image, -dx, -dy, algorithm);
    afwGeom::Point2I llc(buffer, buffer);
    afwGeom::Extent2I dims(width, height);
    afwGeom::Box2I box(llc, dims);
    _offsetImage.reset(new MaskedImageT(*offset, box, afwImage::LOCAL, true)); // Deep copy

    return _offsetImage;
}




/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
typedef float Pixel;
//template class measAlg::PsfCandidate<afwImage::MaskedImage<Pixel> >;
template class measAlg::PsfCandidate<Pixel>;
/// \endcond
