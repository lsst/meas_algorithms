#if !defined(LSST_DETECTION_PSF_H)
#define LSST_DETECTION_PSF_H
//!
// Describe an image's PSF
//
#include "boost/shared_ptr.hpp"
#include "lsst/daf/data.h"
#include "lsst/afw/math.h"

namespace lsst { namespace meas { namespace algorithms {

/*!
 * \brief Represent an image's PSF
 */
class PSF : public lsst::daf::data::LsstBase {
public:
    typedef boost::shared_ptr<PSF> Ptr;

    explicit PSF(lsst::afw::math::Kernel::PtrT kernel=lsst::afw::math::Kernel::PtrT());
    virtual ~PSF() = 0;

    void setA(double const A) { _A = A; } ///< Set the central amplitude
    double getA() const { return _A; }  ///< Get the central amplitude
    ///
    /// Convolve an image with a Kernel
    ///
    template <typename ImageT>
    void convolve(ImageT& convolvedImage,          ///< convolved image
                  ImageT const& inImage,           ///< image to convolve
                  bool doNormalize=true,           ///< if True, normalize the kernel, else use "as is"
                  int edgeBit=-1        ///< mask bit to indicate pixel includes edge-extended data;
                  ///< if negative (default) then no bit is set; only relevant for MaskedImages
                 ) {
        lsst::afw::math::convolve(convolvedImage, inImage, *getKernel(), doNormalize, edgeBit);        
    }

    virtual double getValue(double const dx=0, double const dy=0) const = 0; ///< Evaluate the PSF at (x, y)

    void setKernel(lsst::afw::math::Kernel::PtrT kernel);
    lsst::afw::math::Kernel::PtrT getKernel();
    boost::shared_ptr<const lsst::afw::math::Kernel> getKernel() const;
private:
    lsst::afw::math::Kernel::PtrT _kernel; // Kernel that corresponds to the PSF
    double _A;                          // Central amplitude of PSF
};

/*!
 * \brief Represent a PSF as a circularly symmetrical double Gaussian
 *
 * @note This is just a place holder
 */
class dgPSF : public PSF {
public:
    typedef boost::shared_ptr<dgPSF> Ptr;

    explicit dgPSF() : PSF() {}
    explicit dgPSF(
                   double sigma1,         ///< Width of inner Gaussian
                   double sigma2=1,       ///< Width of outer Gaussian (1 as 0 upsets DoubleGaussianFunction2)
                   double b=0,            ///< Central amplitude of outer Gaussian (inner amplitude == 1)
                   int size=0             ///< Kernel should have dimensions (size*size)
                  );
    double getValue(double const dx=0, double const dy=0) const; ///< Evaluate the PSF at (x, y)
private:
    double _sigma1;                     ///< Width of inner Gaussian
    double _sigma2;                     ///< Width of outer Gaussian
    double _b;                          ///< Central amplitude of outer Gaussian (inner amplitude == 1)
};

}}}

#endif
