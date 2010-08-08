#ifndef MeasAlgoShapelet_H
#define MeasAlgoShapelet_H

/**
 * \file
 *
 * \brief Defines the Shapelet class
 *
 * \author Mike Jarvis
 */

#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/geom/Point.h"
#include "boost/shared_ptr.hpp"
#include "Eigen/Core"
#include <complex>

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletImpl;
    namespace shapelet { class BVec; }

    class Shapelet 
    {
        /*! 
         * \brief This class includes the basic functionality for shapelets
         *
         * There are a few definitions of shapelets out there.  
         * This class implements the shapelets defined in 
         * Bernstein & Jarvis (2002) (where they are called the
         * Gauss-Laguerre basis functions).
         * These are similar to the "polar shapelets" of Massey &
         * Refregier.  
         *
         * Shapelets are a complete basis set that can be used to 
         * describe any image.  However, they are most useful for 
         * describing things that are intrinsically similar to a 
         * 2-d Gaussian.  Since the PSF (for ground-based telescopes)
         * is dominated by a Gaussian component, this makes them
         * a good choice for describing the PSF.
         *
         * For specificity, the basis functions we use are:
         *
         * For p >= q:
         * psi_pq(x,y,sigma) = (pi p! q!)^-1/2 z^m exp(-r^2/2) K_pq(r^2)
         * where z = (x + Iy)/sigma
         *       m = p-q
         *       r = |z|^2
         *       K_pq(r^2) = (-)^q q! L_q^(m) (r^2)
         *       L_q^(m)(x) are Laguerre polynomials
         *
         * For any particular value of sigma, these define a complete
         * set of basis functions.
         *
         * The image is thus decomposed as:
         *
         * f(x,y) = Sum_pq b_pq psi_pq(x,y)
         *
         * As with any basis set, the decomposition is only exact 
         * if an infinite number of basis functions are used.
         * We truncate the expansion at p+q <= order
         * for some (specified) order.
         */

    public:

        typedef boost::shared_ptr<Shapelet> Ptr;
        typedef boost::shared_ptr<const Shapelet> ConstPtr;

        typedef Eigen::VectorXd ShapeletVector;
        typedef Eigen::MatrixXd ShapeletCovariance;

        typedef lsst::afw::detection::Source Source;
        typedef lsst::afw::image::MaskedImage<double> MaskedImage;
        typedef lsst::afw::image::Wcs Wcs;
        typedef lsst::afw::geom::PointD PointD;
        typedef lsst::afw::image::MaskPixel MaskPixel;

        /*!
         * \brief Basic constructor requires order and sigma.
         *
         * order defines how many shapelet coefficients to measure.
         * e.g. order = 4 includes:
         * b00, b10, b20, b11, b30, b21, b40, b31, b22
         * (bii are real, the others are complex.)
         * 
         * sigma is the scale size of Gaussian factor in the shapelet
         * decomposition, measured in arcsec.
         */
        Shapelet(
            int order,      ///< the maximum value of p+q to measure.
            double sigma    ///< the scale size of the shapelet decomposition (arcsec).
        );

        /*!
         * \brief A constructor from a vector of values
         *
         * The input vector should have (order+1)*(order+2)/2 elements.
         *
         * Caveat: The input vector is real, not complex.  
         * See the comment for getValues for the expected order of values.
         */
        Shapelet(
            int order,      ///< the maximum value of p+q to measure.
            double sigma,   ///< the scale size of the shapelet decomposition (arcsec).
            const ShapeletVector& vector ///< the shapelet vector
        );

        /*!
         * \brief A constructor from a vector of values, with covariance
         *
         * The input vector should have (order+1)*(order+2)/2 elements.
         * The covariance matrix is of the input vector of values.
         * So it should be square, symmetric, and have the same size in
         * each dimension as the input vector.
         *
         * Caveat: The input vector is real, not complex.  
         * See the comment for getValues for the expected order of values.
         */
        Shapelet(
            int order,      ///< the maximum value of p+q to measure.
            double sigma,   ///< the scale size of the shapelet decomposition (arcsec).
            const ShapeletVector& vector, ///< the shapelet vector
            const ShapeletCovariance& cov ///< the covariance matrix
        );

        /*!
         * \brief Destructor needs to delete pImpl
         */
        ~Shapelet();

        /*!
         * \brief Copy constructor does a deep copy
         */
        Shapelet(const Shapelet& rhs);

        /*! 
         * \brief op= does a deep copy
         */
        Shapelet& operator=(const Shapelet& rhs);

        /*!
         * \brief get the order of the shapelet
         */
        int getOrder() const;

        /*!
         * \brief get the scale size of the shapelet
         */
        double getSigma() const;

        /*!
         * \brief the size of the shapelet vector
         */
        int size() const;

        /*!
         * \brief get the values as a vector.
         *
         * The order of values are:
         * b00
         * real(b10)  imag(b10)
         * real(b20)  imag(b20)  b11
         * real(b30)  imag(b30)  real(b21)  imag(b21)
         * real(b40)  imag(b40)  real(b31)  imag(b31)  b22
         * etc.
         *
         * where bpq corresponds to the coefficients $b_{pq}$ as defined in
         * Bernstein and Jarvis (2002).  
         */
        const ShapeletVector& getValues() const;

        /*!
         * \brief does the shapelet have a covariance matrix stored?
         */
        bool hasCovariance() const;

        /*!
         * \brief get the covariance matrix
         */
        boost::shared_ptr<const ShapeletCovariance> getCovariance() const;

        /*!
         * \brief set a new value of sigma
         */
        void setSigma(double sigma);

        /*!
         * \brief Get a complex b_pq value 
         *
         * Get the coefficient b_pq for the shapelet decomposition.
         */
        std::complex<double> getPQ(int p, int q);

        /*!
         * \brief Evaluate f(x,y) 
         *
         * Evaluate the intensity as a function of chip position (x,y)
         *
         * This is the approximation to the true f(x,y) defined by 
         * the decomposition Sum_pq b_pq psi_pq(x,y,sigma).
         */
        double evaluateAt(const PointD& pos);
        double evaluateAt(double x, double y);

        /*!
         * \brief measure shapelet decomposition of an image
         *
         * The initial estimate of the centroid is taken from the source.  
         * The initial estimate of sigma is the value given by 
         * this->getSigma().
         *
         * The normal operation would be to allow both values to shift to 
         * find the best values for the image.  However, there are reasons
         * to keep each fixed.
         *
         * First, if we have very accurate positions given in the source 
         * object, and we want the PSF to include the centroid shift that 
         * is present in the seeing, then this would be a good reason to 
         * keep it fixed.
         *
         * Second, interpolation of the PSF is more straightforward
         * if all the stars are measured with the saem sigma.
         * So the normal practice is to let sigma float for each
         * star, then find the mean of all the stars in an image,
         * and finally remeasure each star with this constant
         * value of sigma.
         *
         * So, there are two boolean parameters that govern what 
         * should be kept fixed during the calculation.
         *
         * isCentroidFixed says whether to keep the centoid fixed.  
         * If it is false, then the centoid is allowed to vary until
         * b10 = 0 (our definition of being centroided).
         *
         * isSigmaFixed says whether to keep the sigma fixed.
         * If it is false, then sigma is allowed to vary until 
         * b11 = 0 (which corresponds to when the measured shapelet 
         * components have the maximum signal-to-noise).
         *
         * The measurement is done using a circular aperture,
         * whose radius is given by the paremeter aperture.
         * This aperture (given in arcsec) defines a circle in
         * _sky_ coordinates, not chip coordinates.
         *
         * The return value is true if the measurement is successful,
         * and falso if not.
         */
        bool measureFromImage(
            const Source& source,       ///< The source of which to measure the decomposition.
            const PointD& pos,          ///< The position of the object
            bool isCentroidFixed,       ///< Is sigma fixed? or should it be allowed to vary?
            bool isSigmaFixed,          ///< Is sigma fixed? or should it be allowed to vary?
            double aperture,            ///< The aperture size in arcsec for the measurement.
            const MaskedImage& image,   ///< The image on which to measure the decompoisition.
            const Wcs& wcs,             ///< The wcs to use to convert from x,y to ra,dec.
            const MaskPixel okmask=0    ///< The mask values that are ok to use
        );

        /*!
         * \brief A constructor that is only used by various Shapelet implementations.
         *
         * This is used by ShapeletInterpolation for example.
         * It is more convenient to leave this public, but you can't actually use
         * it, since BVec isn't defined publicly.
         */
        Shapelet(const shapelet::BVec& bvec);

        /*!
         * \brief View the shapelet as a BVec.
         *
         * This is also only used by various Shapelet implementations, but it
         * is convenient to leave it public rather than define a bunch of friend 
         * classes.
         */
        const shapelet::BVec& viewAsBVec() const;
        shapelet::BVec& viewAsBVec();

    private :

        ShapeletImpl* pImpl;
    };

    /*!
     * \brief a helper function to deal with the fact that Wcs doesn't directly return a Jacobian.
     *
     * This function calculates the local Jacobian of the distortion pattern
     * at a given location.
     *
     * The input position is given in chip coordinates, since that's what 
     * we generally have when we want to calculate the local distortion.
     *
     * The return matrix is:
     *
     * J = ( du/dx   du/dy )
     *     ( dv/dx   dv/dy )
     *
     * where (u,v) are the sky coordinates, and (x,y) are the chip coordinates.
     */
    Eigen::Matrix2d getJacobian(
        const lsst::afw::image::Wcs& wcs,       ///< The wcs defining the distortion function
        const lsst::afw::geom::PointD& pos     ///< The position in chip coordinates at which to find the Jacobian
    );

}}}

#endif
