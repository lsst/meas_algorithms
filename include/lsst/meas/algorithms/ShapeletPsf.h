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
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Point.h"
// FIXME: SpatialCell needs to be after Image, since it forgot to 
// #include Image.h, and it uses stuff from it.
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletPsfImpl;

    class ShapeletPsf
    {
    public:
        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::math::SpatialCellSet SpatialCellSet;
        typedef lsst::afw::image::Image<double> Image;
        typedef lsst::afw::image::Wcs Wcs;
        typedef lsst::afw::geom::PointD PointD;
        typedef ShapeletKernel::Ptr KernelPtr;
        typedef LocalShapeletKernel::Ptr LocalKernelPtr;

        typedef double Color;  // A placeholder for some sophisticated Color type

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
         * The cellSet that is input will be modified. (i.e. That Ptr rather 
         * than ConstPtr is not an oversight.)  The Candidates will have their
         * Shapelet information filled in, and the outliers will be marked 
         * as BAD.
         *
         * Lupton and Jarvis discussed having the constructor take an optional Filter argument.
         * I think in this case, the filter info can be grabbed from the star candidates in cellSet.
         * They are Source's, and BaseSourceAttributes has a getFilterId() method.
         * However, if this doesn't work for some reason, then we should add such an argument here.
         */
        ShapeletPsf(
            const Policy& policy,             ///< The policy file for parameters to use.
            SpatialCellSet::Ptr cellSet,      ///< The stars to be measured
            Image::ConstPtr image,            ///< The image on which to measure the decompoisition.
            Wcs::Ptr wcs,                     ///< The wcs to use to convert from x,y to ra,dec.
            Image::ConstPtr weightImage=Image::ConstPtr() ///< If != 0, the image of weight values for pixels.
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
         * \brief op= does a shallow copy
         *
         * After the op=, this will have the same details as rhs, after which 
         * all changes to either one affect the other.
         */
        ShapeletPsf& operator=(const ShapeletPsf& rhs);

        /*!
         * \brief Get the Kernel of the PSF at a given point.
         *
         * The position is given in pixel coordinates.
         *
         * The color term is currently ignored, but is provided for future 
         * implementation
         *
         * The width and height are optional.
         * If not specified, then the Kernel will choose an appropriate size 
         * automatically.
         */
        LocalKernelPtr getLocalKernel(
            const PointD& pos,   ///< position to interpolate to
            Color color,         ///< color to interpolate to
            int width=0,         ///< width of Kernel image if you want to specify something particular.
            int height=0         ///< height of Kernel image if you want to specify something particular.
        ) const;

        /*!
         * \brief Get a general Kernel that varies across an image.
         *
         * The color term is currently ignored, but is provided for future 
         * implementation
         *
         * The width and height are optional.
         * If not specified, then the Kernel will choose an appropriate size 
         * automatically.
         */
        KernelPtr getKernel(
            Color color,     ///< color to interpolate to
            int width=0,     ///< width of Kernel image if you want to specify something particular.
            int height=0     ///< height of Kernel image if you want to specify something particular.
        ) const;

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
         * \brief Convenience funtion to skip having to go through a Kernel
         *
         * If other computeImage syntaxes emerge in Kernel.h, we can corresponding functions here too.
         * Lupton thinks scientists will find this syntax useful.
         */
        inline double computeImage(
            lsst::afw::image::Image<double> &image,   ///< image whose pixels are to be set (output)
            bool doNormalize,   ///< normalize the image (so sum is 1)?
            const PointD& pos,  ///< position to interpolate to
            Color color         ///< color to interpolate to
        ) const 
        { return getLocalKernel(pos,color)->computeImage(image,doNormalize); }


    private :

        ShapeletPsfImpl* pImpl;
    };

}}}

#endif
