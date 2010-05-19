// -*- LSST-C++ -*-
/*!
 * \brief Implementation of PSF code
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <typeinfo>
#include <cmath>
#include "lsst/afw/image/ImagePca.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"

/************************************************************************************************************/

namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

PSF::PSF(int const width,               // desired width of Image realisations of the kernel
         int const height               // desired height of Image realisations of the kernel; default: width
        ) :  lsst::daf::data::LsstBase(typeid(this)),
             _kernel(afwMath::Kernel::Ptr()),
             _width(width), _height(height == 0 ? width : height) {}

PSF::PSF(lsst::afw::math::Kernel::Ptr kernel ///< The Kernel corresponding to this PSF
        ) : lsst::daf::data::LsstBase(typeid(this)),
            _kernel(kernel),
            _width(kernel.get()  == NULL ? 0 : kernel->getWidth()),
            _height(kernel.get() == NULL ? 0 : kernel->getHeight()) {}

/// PSF's destructor; declared pure virtual, but we still need an implementation
PSF::~PSF() {}

///
/// Set the PSF's kernel
///
void PSF::setKernel(lsst::afw::math::Kernel::Ptr kernel) {
    _kernel = kernel;
}

///
/// Return the PSF's kernel
///
afwMath::Kernel::Ptr PSF::getKernel() {
    return _kernel;
}

///
/// Return the PSF's kernel
///
boost::shared_ptr<const afwMath::Kernel> PSF::getKernel() const {
    return boost::shared_ptr<const afwMath::Kernel>(_kernel);
}

/**
 * Return an Image of the the PSF at the point (x, y), setting the sum of all the PSF's pixels to 1.0
 *
 * The specified position is a floating point number, and the resulting image will
 * have a PSF with the correct fractional position, with the centre within pixel (width/2, height/2)
 * Specifically, fractional positions in [0, 0.5] will appear above/to the right of the center,
 * and fractional positions in (0.5, 1] will appear below/to the left (0.9999 is almost back at middle)
 *
 * @note If a fractional position is specified, the central pixel value may not be 1.0
 *
 * @note This is a virtual function; we expect that derived classes will do something
 * more useful than returning a NULL pointer
 */
afwImage::Image<PSF::Pixel>::Ptr PSF::getImage(double const, ///< column position in parent %image
                                                double const  ///< row position in parent %image
                                               ) const {
    return afwImage::Image<PSF::Pixel>::Ptr();
}

/************************************************************************************************************/
/*
 * Register a factory object by name;  if the factory's NULL, return the named factory
 */
PsfFactoryBase& PSF::_registry(std::string const& name, PsfFactoryBase* factory) {
    static std::map<std::string const, PsfFactoryBase *> psfRegistry;

    std::map<std::string const, PsfFactoryBase *>::iterator el = psfRegistry.find(name);

    if (el == psfRegistry.end()) {      // failed to find name
        if (factory) {
            psfRegistry[name] = factory;
        } else {
            throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                              "Unable to lookup Psf variety \"" + name + "\"");
        }
    } else {
        if (!factory) {
            factory = (*el).second;
        } else {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException,
                              "Psf variety \"" + name + "\" is already declared");
        }
    }

    return *factory;
}

/**
 * Declare a PsfFactory for a variety "name"
 *
 * @throws lsst::pex::exceptions::InvalidParameterException if name is already declared
 */
void PSF::declare(std::string name,          ///< name of variety
                  PsfFactoryBase* factory ///< Factory to make this sort of PSF
                 ) {
    (void)_registry(name, factory);
}

/**
 * Return the named PsfFactory
 *
 * @throws lsst::pex::exceptions::NotFoundException if name can't be found
 */
PsfFactoryBase& PSF::lookup(std::string name ///< desired variety
                                 ) {
    return _registry(name, NULL);
}

/************************************************************************************************************/
/**
 * Return a Psf of the requested variety
 *
 * @throws std::runtime_error if name can't be found
 */
PSF::Ptr createPSF(std::string const& name,       ///< desired variety
                   int width,                     ///< Number of columns in realisations of PSF
                   int height,                    ///< Number of rows in realisations of PSF
                   double p0,                     ///< PSF's 1st parameter
                   double p1,                     ///< PSF's 2nd parameter
                   double p2                      ///< PSF's 3rd parameter
            ) {
    return PSF::lookup(name).create(width, height, p0, p1, p2);
}

/**
 * Return a Psf of the requested variety
 *
 * @throws std::runtime_error if name can't be found
 */
PSF::Ptr createPSF(std::string const& name,             ///< desired variety
                   lsst::afw::math::Kernel::Ptr kernel ///< Kernel specifying the PSF
                  ) {
    return PSF::lookup(name).create(kernel);
}


    

/**
 * @brief Constructor for PsfAttributes
 *
 */
PsfAttributes::PsfAttributes(
                             PSF::Ptr psf, ///< The psf whose attributes we want
                             int const iX, ///< the x position in the frame we want the attributes at
                             int const iY  ///< the y position in the frame we want the attributes at
                            )
{
    // N.b. (iX, iY) are ints so that we know this image is centered in the central pixel of _psfImage
    _psfImage = psf->getImage(iX, iY);
}

namespace {

/*
 * Return an estimate of <r> == <sqrt(x^2 + y^2)> for an image (i.e. sum(I*r)/sum(I))
 *
 * For a Gaussian N(0, alpha^2),  <r> = sqrt(pi/2) alpha
 */
template<typename ImageT>
double
computeFirstMoment(ImageT const& image,        // the data to process
                   float const xCen, float const yCen // centre of object
                  )
{
    double sum = 0.0;
    double norm = 0.0;
    for (int iY = 0; iY != image->getHeight(); ++iY) {
        int iX = 0;
        for (afwImage::Image<double>::x_iterator ptr = image->row_begin(iY),
                                                 end = image->row_end(iY); ptr != end; ++ptr, ++iX) {
            double const x = iX - xCen;
            double const y = iY - yCen;
            double const r = std::sqrt( x*x + y*y );
            double const m = (*ptr)*r;
            norm += *ptr;
            sum += m;
        }
    }
    
    std::string errmsg("");
    if (sum < 0.0) {
        errmsg = "sum(I*r) is negative.  ";
    }
    if (norm <= 0.0) {
        errmsg += "sum(I) is <= 0.";
    }
    if (errmsg != "") {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException, errmsg);
    }
    
    return sum/norm;
}

/*
 * Return an estimate of <r^2> == <x^2 + y^2> for an image (i.e. sum(I*r^2)/sum(I))
 *
 * For a Gaussian N(0, alpha^2),  <r^2> = 2 alpha^2
 */
template<typename ImageT>
double
computeSecondMoment(ImageT const& image,        // the data to process
                    float const xCen, float const yCen // centre of object
                   )
{
    double sum = 0.0;
    double norm = 0.0;
    for (int iY = 0; iY != image->getHeight(); ++iY) {
        int iX = 0;
        for (afwImage::Image<double>::x_iterator ptr = image->row_begin(iY),
                                                 end = image->row_end(iY); ptr != end; ++ptr, ++iX) {
            double const x = iX - xCen;
            double const y = iY - yCen;
            double const r2 = x*x + y*y;
            double const m = (*ptr)*r2;
            norm += *ptr;
            sum += m;
        }
    }
    
    std::string errmsg("");
    if (sum < 0.0) {
        errmsg = "sum(I*r*r) is negative.  ";
    }
    if (norm <= 0.0) {
        errmsg += "sum(I) is <= 0.";
    }
    if (errmsg != "") {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException, errmsg);
    }

    return sum/norm;
}
    
/*****************************************************************************/
/*
 * Calculate weighted moments of an object up to 2nd order
 */
template<typename ImageT>
std::pair<bool, double>
calcmom(ImageT const& image,                // the image data
        float const xCen, float const yCen, // centre of object
        double const w11                    // weights
       )
{
    assert(w11 >= 0);                   /* i.e. it was set */
    if (fabs(w11) > 1e6) {
        return std::make_pair(false, std::numeric_limits<double>::quiet_NaN());
    }

    double sum = 0, sumrr = 0.0;

    for (int i = 0; i <= image.getHeight(); ++i) {
        float const y = i - yCen;
        float const y2 = y*y;
        
        typename ImageT::x_iterator ptr = image.row_begin(i);
        for (int j = 0; j <= image.getWidth(); ++j, ++ptr) {
            float const x = j - xCen;
            float const x2 = x*x;
            float const expon = (x2 + y2)*w11;
            
            if (expon <= 14.0) {
                float const weight = exp(-0.5*expon);
                float const tmod = *ptr;
                float const ymod = tmod*weight;
                sum += ymod;
                sumrr += (x2 + y2)*ymod;
            }
        }
    }

    if (sum <= 0 || sumrr < 0) {
        return std::make_pair(false, std::numeric_limits<double>::quiet_NaN());
    }

    return std::make_pair(true, 0.5*sumrr/sum); // 0.5:  1-D moment
}

/*
 * Return an estimate of <r^2> == <x^2 + y^2> for an image using adaptive moments
 *
 * For a Gaussian N(0, alpha^2),  <r^2> = 2 alpha^2
 *
 * This is basically the SdssShape code simplified for a circularly symmetrical case.  I don't want to call
 * the shape code here as this class may well be moving to afw along with Psf
 */
template<typename ImageT>
double
computeSecondMomentAdaptive(ImageT const& image,        // the data to process
                            float const xCen, float const yCen // centre of object
                      )
{
    int const MAXIT = 100;              // \todo from Policy XXX
    float const TOL = 0.0001;
    double sigma11W = 1.5;              // quadratic moment of the weighting function
    double w11 = -1;                    // current weight for moments; always set when iter == 0
    float sigma11_ow_old = 1e6;         // previous version of sigma11_ow
   
    bool unweighted = false;            // do we need to use an unweighted moment?
    int iter = 0;                       // iteration number
    for (; iter < MAXIT; ++iter) {
        assert(sigma11W > std::numeric_limits<float>::epsilon());
        w11 = 1/sigma11W;

        std::pair<bool, double> moments = calcmom(*image, xCen, yCen, w11);
        
        if (not moments.first) {
            unweighted = true;
            break;
        }
/*
 * Did we converge?
 */
        float const sigma11_ow = moments.second; // quadratic moments of weight*object

        if (iter > 0 && fabs(sigma11_ow/sigma11_ow_old - 1.0) < TOL) {
            break;                              /* yes; we converged */
        }

        sigma11_ow_old = sigma11_ow;
/*
 * Didn't converge, calculate new values for weighting function
 *
 * The product of two Gaussians is a Gaussian, the inverse-variances add
 *
 * We know sigma11_ow and sigma11W, the variances of the weighted object
 * and of the weights themselves.  We can estimate the object's variance as
 *   1/sigma11_ow - 1/sigma11W
 * and, as we want to find a set of weights with the _same_ covariance as the
 * object we take this to be the an estimate of our correct weights.
 *
 * N.b. This assumes that the object is roughly Gaussian.
 * Consider the object:
 *   O == delta(x + p) + delta(x - p)
 * the covariance of the weighted object is equal to that of the unweighted
 * object, and this prescription fails badly.  If we detect this, we set
 * unweighted, and calculate the UNweighted moments
 * instead.
 */
        {
            if (sigma11_ow <= 0) {
                unweighted = true;
                break;
            }
         
            float const ow11 =  1/sigma11_ow;
            float const n11 = ow11 - w11; // inverse of new sigma11_ow

            if (n11 <= 0) {             // product-of-Gaussians assumption failed
                unweighted = true;
                break;
            }
      
            sigma11W = 1/n11;
        }
    }
/*
 * Problems; try calculating the un-weighted moments
 */
    if (iter == MAXIT || unweighted) {
        w11 = 0;
        std::pair<bool, double> moments = calcmom(*image, xCen, yCen, w11);

        if (moments.first) {
            sigma11W = moments.second;  // estimate of object moment
        } else {
            sigma11W = 1/12.0;          // a single pixel
        }
    }

    return 2*sigma11W;
}

}
    
/**
 * @brief Compute the 'sigma' value for an equivalent gaussian psf.
 *
 */
double PsfAttributes::computeGaussianWidth(PsfAttributes::Method how) {
    double const xCen = _psfImage->getWidth()/2;
    double const yCen = _psfImage->getHeight()/2;

    switch (how) {
      case PsfAttributes::ADAPTIVE:
        return ::sqrt(0.5*computeSecondMomentAdaptive(_psfImage, xCen, yCen));
      case PsfAttributes::FIRST_MOMENT:
        return ::sqrt(2.0/M_PI)*computeFirstMoment(_psfImage, xCen, yCen);
      case PsfAttributes::SECOND_MOMENT:
        return ::sqrt(0.5*computeSecondMoment(_psfImage, xCen, yCen));
      case PsfAttributes::BICKERTON:
        double sum = 0.0;
        double norm = 0.0;
        for (int iY = 0; iY != _psfImage->getHeight(); ++iY) {
            int iX = 0;
            afwImage::Image<double>::x_iterator end = _psfImage->row_end(iY);
            for (afwImage::Image<double>::x_iterator ptr = _psfImage->row_begin(iY); ptr != end; ++ptr, ++iX) {
                double const x = iX - xCen;
                double const y = iY - yCen;
                double const r = std::sqrt( x*x + y*y );
                double const m = (*ptr)*r;
                norm += (*ptr)*(*ptr);
                sum += m*m;
            }
        }
        return sqrt(sum/norm);
    }
}
    
/**
 * @brief Compute the effective area of the psf ( (sum(I))^2 / sum(I^2) )
 *
 */
double PsfAttributes::computeEffectiveArea() {
    
    double sum = 0.0;
    double sumsqr = 0.0;
    for (int iY = 0; iY != _psfImage->getHeight(); ++iY) {
        afwImage::Image<double>::x_iterator end = _psfImage->row_end(iY);
        for (afwImage::Image<double>::x_iterator ptr = _psfImage->row_begin(iY); ptr != end; ++ptr) {
            sum += *ptr;
            sumsqr += (*ptr)*(*ptr);
        }
    }
    return sum*sum/sumsqr;
}

    
    
}}}


//template lsst::meas::algorithms::PsfAttributes::PsfAttributes<lsst::meas::algorithms::details::dgPSF>(lsst::meas::algorithms::PSF::Ptr psf, double const, double const);

