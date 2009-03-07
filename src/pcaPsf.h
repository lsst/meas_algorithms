#if !defined(LSST_DETECTION_PCAPSF_H)
#define LSST_DETECTION_PCAPSF_H
//!
// Describe an image's PSF
//
#include "lsst/meas/algorithms/PSF.h"
#include "PSFImpl.h"

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * \brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class pcaPSF : public PSF {
public:
    typedef boost::shared_ptr<pcaPSF> Ptr;

    /**
     * @brief constructors for a pcaPSF
     *
     * Parameters:
     */
    explicit pcaPSF(lsst::afw::math::Kernel::PtrT kernel);

    lsst::afw::image::Image<PSF::PixelT>::Ptr getImage(double const x, double const y) const;
private:
    double doGetValue(double const dx, double const dy, int const xPositionInImage, int const yPositionInImage) const;
};

}}}
#endif
