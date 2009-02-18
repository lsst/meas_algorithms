#if !defined(LSST_DETECTION_PSF_H)
#define LSST_DETECTION_PSF_H
//!
// Describe an image's PSF
//
#include "boost/shared_ptr.hpp"
#include "lsst/daf/data.h"
#include <lsst/pex/exceptions.h>
#include "lsst/afw/math.h"

namespace lsst { namespace meas { namespace algorithms {

/**
 * Types of supported PSFs
 */
typedef int psfType;

/*!
 * \brief Represent an image's PSF
 */
class PSF : public lsst::daf::data::LsstBase {
public:
    typedef boost::shared_ptr<PSF> Ptr; ///< shared_ptr to a PSF
    typedef boost::shared_ptr<const PSF> ConstPtr; ///< shared_ptr to a const PSF

    explicit PSF(lsst::afw::math::Kernel::PtrT kernel=lsst::afw::math::Kernel::PtrT());
    virtual ~PSF() = 0;
    ///
    /// Convolve an image with a Kernel
    ///
    template <typename ImageT>
    void convolve(ImageT& convolvedImage,          ///< convolved image
                  ImageT const& inImage,           ///< image to convolve
                  bool doNormalize=true,           ///< if True, normalize the kernel, else use "as is"
                  int edgeBit=-1        ///< mask bit to indicate pixel includes edge-extended data;
                  ///< if negative (default) then no bit is set; only relevant for MaskedImages
                 ) const {
        if (!getKernel() || getKernel()->getWidth() <= 0 || getKernel()->getHeight() <= 0) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                              "PSF does not have a realisation that can be used for convolution");            
        }
        lsst::afw::math::convolve(convolvedImage, inImage, *getKernel(), doNormalize, edgeBit);        
    }

    ///< Evaluate the PSF at (dx, dy)
    ///
    /// This routine merely calls doGetValue, but here we can provide default values
    /// for the virtual functions that do the real work
    ///
    double getValue(double const dx=0,  ///< column position (relative to centre of PSF)
                    double const dy=0   ///< row position (relative to centre of PSF)
                   ) const {
        return doGetValue(dx, dy);
    }

    void setKernel(lsst::afw::math::Kernel::PtrT kernel);
    lsst::afw::math::Kernel::PtrT getKernel();
    boost::shared_ptr<const lsst::afw::math::Kernel> getKernel() const;

    static psfType lookupType(std::string const& name);
protected:
    static void registerType(std::string const& name, psfType type);
private:
    virtual double doGetValue(double const dx, double const dy) const = 0;
    static std::map<std::string, psfType>* _psfTypes;

    lsst::afw::math::Kernel::PtrT _kernel; // Kernel that corresponds to the PSF
};

/************************************************************************************************************/
/**
 * Factory functions to return a PSF
 */
PSF *createPSF(std::string const& type, int size=0, double=0, double=0, double=0);
}}}
#endif
