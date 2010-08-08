#ifndef MeasAlgoShapeletPsf_H
#define MeasAlgoShapeletPsf_H

/**
 * \file
 *
 * \brief A PSF class that describes the PSF in terms of its shapelet decomposition.
 *
 * \author Mike Jarvis
 */

#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"
//#include "lsst/meas/algorithms/PSF.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletPsfImpl;

    class ShapeletPsf : public afw::detection::Psf
    {
    public:
        typedef afw::detection::Psf Base;

        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::math::SpatialCellSet SpatialCellSet;
        typedef lsst::afw::image::MaskedImage<double> MaskedImage;
        typedef lsst::afw::image::Image<double> Image;
        typedef lsst::afw::image::Wcs Wcs;
        typedef lsst::afw::geom::Point2D Point;
        typedef lsst::afw::geom::Extent2I Extent;
        typedef lsst::afw::image::Color Color;
        typedef lsst::afw::math::Kernel Kernel;

        typedef boost::shared_ptr<ShapeletPsf> Ptr;
        typedef boost::shared_ptr<const ShapeletPsf> ConstPtr;

        /*!
         * \brief Construct the ShapeletPsf from a list of Sources
         *
         * Here is a list of parameters that should be in the policy object 
         * (along with some suggested values)
         *
         * shapeletOrder      int        10     The order of the shapelet measurements.
         * shapeletSigma      double     -1.    The sigma to use.  If <= 0, then determine from the data.
         * psfAperture        double     5.     The aperture radius in arcsec to use.
         * nStarsPerCell      int        5      The max number of stars to use per cell.
         * interpOrder        int        2      The order of the polynomial fit in (x,y)
         * interpNSigmaClip   double     3.     The number of sigma to use for outlier rejection.
         * pcaThresh          double     1.e-5  The theshold value for which principal components to keep.
         * colorTerm          string     "r-i"  ** Not implemented.  Need some way to define what color to use.
         *
         * The candidates in cellSet must be ShapeletPsfCandidate's.
         *
         * The cellSet that is input will be modified. 
         * The Candidates will have their Shapelet information filled in, and the outliers will be marked as BAD.
         * So the method getCellSet below returns something that is _NOT_ equivalent to the input cellSet here.
         *
         * Lupton and Jarvis discussed having the constructor take an optional Filter argument.
         * I think in this case, the filter info can be grabbed from the star candidates in cellSet.
         * They are Source's, and BaseSourceAttributes has a getFilterId() method.
         * However, if this doesn't work for some reason, then we should add such an argument here.
         */
        ShapeletPsf(
            const Policy& policy,            ///< The policy file for parameters to use.
            const SpatialCellSet& cellSet,   ///< The stars to be measured
            const MaskedImage& image,        ///< The image on which to measure the decompoisition.
            const Wcs& wcs                   ///< The wcs to use to convert from x,y to ra,dec.
        );

        /*!
         * \brief Destructor needs to delete pImpl
         */
        ~ShapeletPsf();

        /*!
         * \brief Copy constructor does a shallow copy.
         *
         * The copy shares the details with the original, so all changes to 
         * either one affect the other.
         */
        ShapeletPsf(const ShapeletPsf& rhs);

        /*!
         * \brief Make a clone of this
         */
        inline Base::Ptr clone() 
        { return boost::make_shared<ShapeletPsf>(*this); }

        /*!
         * \brief Get the cellSet of Psf candidates used for the interpolation.
         *
         * This isn't the same as the cellSet that was input in the constructor, because:
         *
         * 1) The candidates now have measured Shapelets
         * 2) Some candidates are marked bad if they were deemed to be outliers.
         */
        const SpatialCellSet& getCellSet() const;

        /*!
         * \brief Get the preferred size of the image.
         *
         * If you are flexible about the size of the image version of the
         * kernel, then in routines like getKernel (in the base class Psf), you
         * can omit the size argument, and just let the Psf class pick out a
         * good size to use.  The way the base class implements this is to call
         * this function to get a good size to use.
         * FIXME: This doesn't seem to be the case right now.  Should I get rid
         * of this method, or should be implement this as a virtual method
         * in the base class?
         *
         * In our case, the size is chosen to be 10 sigma in each direction,
         * where sigma is the shapelet scale used for measuring the Psf (which
         * is designed to be optimal for the average star in the input cellSet).
         */
        Extent getDefaultExtent() const;

        /*!
         * \brief Get the Kernel of the PSF at a given point.
         *
         * The position is given in pixel coordinates.
         *
         * The color term is currently ignored, but is provided for future 
         * implementation
         *
         */
        Kernel::ConstPtr doGetLocalKernel(
            const Point& ccdXY,     ///< position to interpolate to
            const Color& color      ///< color to interpolate to
        ) const;

        Kernel::Ptr doGetLocalKernel(
            const Point& ccdXY,     ///< position to interpolate to
            const Color& color      ///< color to interpolate to
        );

        /*!
         * \brief Get a general Kernel that varies across an image.
         *
         * [ Required virtual function from base class Psf ]
         *
         * The color term is currently ignored, but is provided for future 
         * implementation
         *
         * The width and height are optional.
         * If not specified, then the Kernel will choose an appropriate size 
         * automatically.
         */
        Kernel::ConstPtr doGetKernel(
            const Color& color    ///< color to interpolate to
        ) const;

        Kernel::Ptr doGetKernel(
            const Color& color    ///< color to interpolate to
        );

        /*! 
         * Use the base class implemenation of doComputeImage.
         * It might be slightly inefficient because it does color 
         * interpolation first.  Then computes the image of the 
         * kernel at a particular point.
         *
         * Probably faster to interpolate on both color and position.
         * Then compute the image from the local kernel.
         *
         * But there are weird things with the normalize parameter
         * that I think are poor design, so I don't particularly want to 
         * duplicate them here.  Specifically, the normalizePeak
         * argument optionally sets the center value of the image to 1.
         *
         * First, it is erroneously (or perhaps presumptuously) called the 
         * peak, whereas not all kernels will have the peak at the center
         * (e.g. comatic PSF).
         *
         * Second, this is a different convention than was proposed in 
         * Kernel where the bool argument optionally sets the _sum_ to 1.
         *
         * I think any normalization operation should be moved to being 
         * a method of the Image class, and let the caller subsequently call
         * that method rather than use the confusing boolean argument here.
         * i.e.
         * image = computeImage(color, pos, size);
         * image.normalizeCenterToUnity();
         *
         * Much clearer than:
         * image = computeImage(color, pos, size, true) 
         *   // true here means normalize center pixel to unity
         */
#if 0
        virtual Image::Ptr doComputeImage(
            const Color& color,       ///< color to interpolate to
            const Point& ccdXY,       ///< position to interpolate to
            const Extent& size,       ///< width/height of Kernel image
            bool normalizePeak        ///< normalize the image (so center is 1)?
        ) const;
#endif


    private :

        ShapeletPsfImpl* pImpl;

        // Not implemented.
        ShapeletPsf& operator=(const ShapeletPsf& rhs);

    };

}}}

#endif
