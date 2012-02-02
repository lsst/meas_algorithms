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
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/afw/geom/Box.h"
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
}

/// Extract an image of the candidate.
///
/// No offsets are applied.
/// The INTRP bit is set for any pixels that are detected but not part of the Source
template <typename PixelT>
PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::extractImage(
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

    /*
     * Set the INTRP bit for any DETECTED pixels other than the one in the center of the object;
     * we grow the Footprint a bit first
     */
    typedef afwDetection::FootprintSet<int>::FootprintList FootprintList;

    MaskPixel const detected = MaskedImageT::Mask::getPlaneBitMask("DETECTED");
    PTR(afwImage::Image<int>) mim = makeImageFromMask<int>(*image->getMask(), makeAndMask(detected));
    PTR(afwDetection::FootprintSet<int>) fs =
        afwDetection::makeFootprintSet<int, MaskPixel>(*mim, afwDetection::Threshold(1));
    CONST_PTR(FootprintList) feet = fs->getFootprints();

    if (feet->size() <= 1) {         // only one Footprint, presumably the one we want
        return image;
    }

    MaskPixel const intrp = MaskedImageT::Mask::getPlaneBitMask("INTRP"); // bit to set for bad pixels
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
        afwDetection::setMaskFromFootprint(image->getMask().get(), *bigfoot, intrp);
    }

    return image;
}


/**
 * Return the %image at the position of the Source, without any sub-pixel shifts to put the centre of the
 * object in the centre of a pixel (for that, use getOffsetImage())
 *
 */
template <typename PixelT>
CONST_PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::getImage(int width, int height) const {

    if (_haveImage && (width != _image->getWidth() || height != _image->getHeight())) {
        _haveImage = false;
    }

    if (!_haveImage) {
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
CONST_PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::getImage() const {

    int const width = getWidth() == 0 ? _defaultWidth : getWidth();
    int const height = getHeight() == 0 ? _defaultWidth : getHeight();

    return getImage(width, height);
    
}


/**
 * @brief Return an undistorted offset version of the image of the source.
 * The returned image has been offset to put the centre of the object in the centre of a pixel.
 *
 */
template <typename PixelT>
PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::getUndistOffsetImage(
    std::string const algorithm,        ///< Warping algorithm to use
    unsigned int offsetBuffer,          ///< Buffer for warping
    bool keepEdge                       ///< Keep a warping edge so image can be distorted again later.
                                                                                            ) const {

    // if we don't have a detector and distortion object, just give them a regular offset image
    if ((!_haveDetector) || (!_haveDistortion)) {
        return getOffsetImage(algorithm, offsetBuffer);
    }
    
    int width  = getWidth() == 0  ? _defaultWidth : getWidth();
    int height = getHeight() == 0 ? _defaultWidth : getHeight();
    
    if (_haveUndistOffsetImage &&
        (width != _undistOffsetImage->getWidth() || height != _undistOffsetImage->getHeight())) {
        _haveUndistOffsetImage = false;
    }


    // When getUndistImage() is called, it will need to compute an oversized image before distorting
    // the image will then be trimmed to the requested size.
    // However, we may intend to redistort things later, and that means we need to keep
    // that extra edge available to be trimmed-off later.
    double const xcen = getXCenter(), ycen = getYCenter();
    if (keepEdge) {
        int edge = std::abs(0.5*((height > width) ? height : width) *
                            (1.0-_distortion->computeMaxShear(*_detector)));
        width += 2*edge;
        height += 2*edge;
    }
    
    if (! _haveUndistOffsetImage) {
        // undistort
        //make it a bit bigger for the offset warp
        PTR(MaskedImageT) undistImgTmp = this->getUndistImage(width + 2*offsetBuffer,
                                                              height + 2*offsetBuffer);
            
        // offset subpixel shift
        double const dx = afwImage::positionToIndex(xcen, true).second;
        double const dy = afwImage::positionToIndex(ycen, true).second;
        PTR(MaskedImageT) undistOffsetImgTmp = afwMath::offsetImage(*undistImgTmp, -dx, -dy, algorithm);
        
        afwGeom::Point2I llc(offsetBuffer, offsetBuffer);
        afwGeom::Extent2I dims(width, height);
        afwGeom::Box2I box(llc, dims);
        _undistOffsetImage.reset(new MaskedImageT(*undistOffsetImgTmp, box, afwImage::LOCAL, true)); //Deep cp
        _haveUndistOffsetImage = true;
    }

    return _undistOffsetImage;
}


/**
 * @brief Return an undistorted version of the image of the source.
 *
 * Here, we mimic the original getImage() call which uses default parameters
 *
 */
template <typename PixelT>
PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::getUndistImage() const {
    return getUndistImage(getWidth(), getHeight());
}

/**
 * @brief Return the *undistorted* %image at the position of the Source,
 * without any sub-pixel shifts to put the centre of the
 * object in the centre of a pixel (for that, use getOffsetImage())
 *
 */
template <typename PixelT>
PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::getUndistImage(
                                                                                       int width,
                                                                                       int height
                                                                                      ) const {


    // if we don't have a detector and distortion object, just give them a regular image
    if ((!_haveDetector) || (!_haveDistortion)) {
        return this->extractImage(width, height);
    }
    
    // if we have it, return it
    if (_haveUndistImage &&
        (width == _undistImage->getWidth() && height == _undistImage->getHeight())) {
        return _undistImage;
        
    // ok.  do the work
    } else {

        // undistort
        
        int distBuffer = _distortion->getLanczosOrder() + 1;
        double const xcen = getXCenter(), ycen = getYCenter();
        afwGeom::Point2D pPixel(xcen, ycen);
        // use an Extent here so we can use /pixelSize() to convert to pixels (Point has no operator/())

        // add on the order of the lanczos kernel ... we're guaranteed to lose that much.
        int edge = distBuffer;

        edge += std::abs(0.5*((height > width) ? height : width) *
                         (1.0 - _distortion->computeMaxShear(*_detector)));
        int widthInit  = width + 2*edge;
        int heightInit = height + 2*edge;
        
        // get a raw image

        PTR(MaskedImageT) imgTmp = this->extractImage(widthInit, heightInit);

        // undistort it at pBoreSight *position*, with pix x,y pixel coordinate
        PTR(MaskedImageT) undistImgTmp =
            _distortion->undistort(pPixel, *imgTmp, *_parentExposure->getDetector());

        // trim the warping buffer
        afwGeom::Point2I llc(edge, edge);
        afwGeom::Extent2I dims(width, height);
        afwGeom::Box2I box(llc, dims);
        _undistImage.reset(new MaskedImageT(*undistImgTmp, box, afwImage::LOCAL, true)); // Deep copy
            
        // remember
        _haveUndistImage = true;
    }
    return _undistImage;
}


/**
 * @brief Return an offset version of the image of the source.
 * The returned image has been offset to put the centre of the object in the centre of a pixel.
 *
 */

template <typename PixelT>
PTR(afwImage::MaskedImage<PixelT,afwImage::MaskPixel,afwImage::VariancePixel>) measAlg::PsfCandidate<PixelT>::getOffsetImage(
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
