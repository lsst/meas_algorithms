#ifndef MeasAlgoStarFinder_H
#define MeasAlgoStarFinder_H

/**
 * \file
 *
 * \brief A module for determining which objects are good PSF stars
 *
 * \author Mike Jarvis
 */

#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/math/SpatialCell.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class StarFinderImpl;

    class StarFinder 
    {
    public:

        typedef lsst::afw::image::MaskedImage<double> MaskedImage;
        //TODO: Does the image need to be templated?  Or is it ok to use double?
        typedef lsst::afw::image::Wcs Wcs;
        typedef lsst::afw::geom::PointD PointD;
        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::detection::Source Source;
        typedef lsst::afw::detection::SourceSet SourceSet;
        typedef lsst::afw::math::SpatialCellSet SpatialCellSet;


        /*!
         * \brief Read the parameters from the policy file.
         *
         * Here is a list of parameters that should be in the policy object
         * along with defaults for when it doesn't find them:
         * minSize      double  0.          Minimum size to use.
         * maxSize      double  1.e100      Maximum size to use.
         * isSizeLog    bool    false       Are sizes alread log(size)?
         * minMag       double  -1.e100     Minimum mag to use.
         * maxMag       double  1.e100      Maximum mag to use.
         * starFrac     double  0.5         What fraction of objects are likely stars?
         * startN       double  0.1         Fraction of objects to use in first pass.
         * fitOrder     int     1           Order of polynomial fit of size(x,y).
         * fitSigClip   double  4.0         nSigma to reject a star as an outlier.
         * fitStars     int     30          Do size(x,y) fit with fitStars brightest stars.
         * purity       double  0.05        Smaller = purer sample of stars, larger = more stars
         * cellSize     int     256         Size of SpatialCellSet cells (in pixels)
         * aperture     double  5.          Aperture size in arcsec.
         *
         * I don't know how to provide a default value for a policy.getDouble()
         * or similar call without an onerous try/catch block for each parameter.
         * So currently, these are actually all required to be in the policy object.  
         *
         * As far as selecting good values for the parameters, I recommend starting 
         * with the above default values and see how that goes.
         */
        StarFinder(const Policy& policy);

        /*!
         * \brief Destructor needs to delete pImpl
         */
        ~StarFinder();

        /*!
         * \brief Calculates a robust size measurement for a source.
         *
         * This measures the 2nd order shapelet decomposition of the object,
         * allowing the sigma to vary until the b11 term goes to zero.
         * When this happens, the sigma basically the best-fit Gaussian sigma.
         * Also, this value of sigma gives the best S/N properties for the
         * rest of the shapelet vector.
         *
         * This size measurement has proven to be a good one to use for 
         * the size-magnitude star finder algorithm, since the size is very
         * stable with stellar magnitude.
         *
         * \note This may not be necessary.  Lupton says that there is already
         * a similar measurement being done in the LSST stack.  Should check
         * if the existing values are ok for star finder algorithm.
         * (I'm not sure which entry in Source they would be...)
         * Therefore, this function may be a good candidate for a having
         * its action be specifyable by a Policy parameter.
         */
        double calculateSourceSize(
            const Source& source, 
            const MaskedImage& image,
            const Wcs& wcs) const;

        /*!
         * \brief Calculates a magnitude for a source.
         *
         * The star finder is written in terms of using magnitudes rather than
         * fluxes, wherease Source seems to store fluxes (specifically PetroFluxes).
         * So this just translates the flux into a magnitude.
         *
         * The normalization doesn't matter, so don't worry about the units
         * unless you specify the minMag and maxMag parameters to be something
         * other than (effectively) +-infinity.
         *
         * \note This function may also be a good candidate for a having
         * its action be specifyable by a Policy parameter.
         */
        double calculateSourceMagnitude(const Source& source) const;

        /*!
         * \brief Get the x and y values for the source
         *
         * This is put here as a method to make is clear which x,y values we 
         * are useing.  My guess is that XAstrom and YAstrom are the right ones,
         * but we might want something else, or even control it with the policy file.
         */
        double getSourceX(const Source& source) const;
        double getSourceY(const Source& source) const;

        /*!
         * \brief Find a set of stars from an input list of Sources
         *
         * This function uses an algorithm based on looking for a stellar locus
         * in the size-magnitude diagram that has constant size.  
         *
         * It basically construct a histogram of counts with respect to size,
         * starting with just the bright objects and then pushing down in magnitude
         * until the stellare peak starts to bleed into the galaxies, at 
         * which point it stops.
         *
         * The image and wcs parameters are only necessary because we currently
         * remeasure the size.  If we have access to a size in the Source that is
         * already good enough, then we can remove the dependency of this file
         * on Image nad Wcs.
         */
        SpatialCellSet::Ptr findStars(
            const SourceSet& allObj,    ///< The input list of objects to consider
            const MaskedImage& image,   ///< The image on which the sources are found
            const Wcs& wcs              ///< The wcs to convert from x,y to ra,dec
        ) const;

    private :

        /*!
         * \brief This class is not intended to be copied.
         */
        StarFinder(const StarFinder& rhs);
        void operator=(const StarFinder& rhs);

        StarFinderImpl* pImpl;
    };

}}}

#endif
