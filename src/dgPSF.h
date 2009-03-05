#if !defined(LSST_DETECTION_DGPSF_H)
#define LSST_DETECTION_DGPSF_H
//!
// Describe an image's PSF
//
#include "lsst/meas/algorithms/PSF.h"
#include "PSFImpl.h"

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * \brief Represent a PSF as a circularly symmetrical double Gaussian
 */
class dgPSF : public PSF {
public:
    typedef boost::shared_ptr<dgPSF> Ptr;

    /**
     * @brief constructors for a dgPSF
     *
     * Parameters:
     */
    explicit dgPSF(int width, int height, double sigma1, double sigma2=1, double b=0);

    lsst::afw::image::Image<PSF::PixelT>::Ptr getImage(double const x, double const y) const;
private:
    double doGetValue(double const dx, double const dy) const;

    double _sigma1;                     ///< Width of inner Gaussian
    double _sigma2;                     ///< Width of outer Gaussian
    double _b;                          ///< Central amplitude of outer Gaussian (inner amplitude == 1)
};

}}}
#endif
