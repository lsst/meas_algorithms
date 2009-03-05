/*!
 * \brief Implementation of PSF code
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <typeinfo>
#include <cmath>
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"

/************************************************************************************************************/
/*
 * Include concrete implementations
 */
#include "PSFImpl.h"
#include "dgPSF.h"

namespace afwImage = lsst::afw::image;

namespace lsst {
namespace meas {
namespace algorithms {

PSF::PSF(int const width,               // desired width of Image realisations of the kernel
         int const height               // desired height of Image realisations of the kernel; default: width
        ) :  lsst::daf::data::LsstBase(typeid(this)),
        _kernel(lsst::afw::math::Kernel::PtrT()),
        _width(width), _height(height == 0 ? width : height) {}


PSF::PSF(lsst::afw::math::Kernel::PtrT kernel ///< The Kernel corresponding to this PSF
        ) : lsst::daf::data::LsstBase(typeid(this)),
            _kernel(kernel),
            _width(kernel->getWidth()), _height(kernel->getHeight()) {
    ;
}

/// PSF's destructor; declared pure virtual, but we still need an implementation
PSF::~PSF() {}

///
/// Set the PSF's kernel
///
void PSF::setKernel(lsst::afw::math::Kernel::PtrT kernel) {
    _kernel = kernel;
}

///
/// Return the PSF's kernel
///
lsst::afw::math::Kernel::PtrT PSF::getKernel() {
    return _kernel;
}

///
/// Return the PSF's kernel
///
boost::shared_ptr<const lsst::afw::math::Kernel> PSF::getKernel() const {
    return boost::shared_ptr<const lsst::afw::math::Kernel>(_kernel);
}

/**
 * Return an Image of the the PSF at the point (x, y), setting the PSF's peak value to 1.0
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
lsst::afw::image::Image<PSF::PixelT>::Ptr PSF::getImage(double const x, ///< column position in parent %image
                                                        double const y  ///< row position in parent %image
                   ) const {
    return lsst::afw::image::Image<PSF::PixelT>::Ptr();
}
    
/************************************************************************************************************/
/**
 * @brief The mapping between type names (e.g. "DGPSF") and an enum (DGPSF)
 */
std::map<std::string, psfType>* PSF::_psfTypes = NULL;

/**
 * @brief Register a (name, enum) pair.
 *
 * This routine should only be called by createPSF
 */
void PSF::registerType(std::string const&name, psfType type) {
    if (_psfTypes == NULL) {
        _psfTypes = new(std::map<std::string, psfType>);
    }

    (*_psfTypes)[name] = type;
}

/**
 * @brief Return the typename for this PSF
 *
 * Names are registered using registerType
 */
psfType PSF::lookupType(std::string const& name ///< Name of this type of PSF
                       ) {
    assert (_psfTypes != NULL);
    
    std::map<std::string, psfType>::const_iterator i = _psfTypes->find(name);
    if (i == _psfTypes->end()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException,
                          (boost::format("Unknown psf algorithm: %s") % name).str());
    }

    return i->second;
}

/**
 * @brief A factory function to return a PSF of the specified type, given as a string.
 */
PSF* createPSF(std::string const& type,           ///< desired type
               int width,                         ///< Number of columns in realisations of PSF
               int height,                        ///< Number of rows in realisations of PSF
               double p0,                         ///< PSF's 1st parameter
               double p1,                         ///< PSF's 2nd parameter
               double p2                          ///< PSF's 3rd parameter
              ) {
    switch (PSF::lookupType(type)) {
      case DGPSF:
        return new dgPSF(width, height, p0, p1, p2);
      default:
        throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException, 
                          (boost::format("PSF of type %d is not implemented") % type).str());
    }
    // NOTREACHED
}

/************************************************************************************************************/
/*
 * PsfCandidate's members
 */
/**
 * Return the %image at the position of the Source
 */
template <typename ImageT>
typename ImageT::ConstPtr lsst::meas::algorithms::PsfCandidate<ImageT>::getImage() const {
    if (!_haveImage) {
        int const width = getWidth() == 0 ? 15 : getWidth();
        int const height = getHeight() == 0 ? 15 : getHeight();
    
        std::pair<int, double> xCen = afwImage::positionToIndex(getXCenter(), true); // true => return the std::pair
        std::pair<int, double> yCen = afwImage::positionToIndex(getYCenter(), true);

        afwImage::PointI center(xCen.first, yCen.first); // integral part
        afwImage::BBox bbox(center - afwImage::PointI(width/2, height/2), width, height);
        bbox.shift(-_parentImage->getX0(), -_parentImage->getY0());
        
        typename ImageT::Ptr patch(new ImageT(*_parentImage, bbox, false)); // a shallow copy

        std::string algorithmName = "lanczos5";
        _image = lsst::afw::math::offsetImage(algorithmName, *patch, -xCen.second, -yCen.second);
    }
    
    return _image;
}

    

//
// Explicit instantiations
//
// \cond
    template class PsfCandidate<lsst::afw::image::MaskedImage<float> >;
// \endcond
            
}}}
