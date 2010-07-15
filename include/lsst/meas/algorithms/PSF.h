// -*- LSST-C++ -*-
#if !defined(LSST_DETECTION_PSF_H)
#define LSST_DETECTION_PSF_H
//!
// Describe an image's PSF
//
#include <string>
#include "lsst/base.h"
#include "lsst/daf/data.h"

namespace lsst {
namespace afw {
    namespace detection {
        class Psf;
    }
    namespace image {
        template<typename T> class Image;
    }
}
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * @class PsfAttributes
 *
 * A class to contain various attributes of the Psf
 * - most notably, a width (1-D RMS size) to be used to
 *   make a single gaussian psf for fast convolution.
 */
class PsfAttributes {
public:
    enum Method { ADAPTIVE_MOMENT,      ///< Calculate width using adaptive Gaussian weights
                  FIRST_MOMENT,         ///< Calculate width using \<r>
                  SECOND_MOMENT,        ///< Calculate width using \<r^2>
                  NOISE_EQUIVALENT,     ///< Calculate width as sqrt(n_eff/(4 pi))
                  BICKERTON             ///< Weight \<r^2> by I^2 to avoid negative fluxes
    };

    PsfAttributes(PTR(lsst::afw::detection::Psf) psf, int const iX, int const iY);
    
    double computeGaussianWidth(Method how=ADAPTIVE_MOMENT);
    double computeEffectiveArea();
    
private:
#define SPACE                           /* Macro PTR is fixed in base 3.1.3 */
    PTR(lsst::afw::image::Image<double> SPACE ) _psfImage;
};

    
}}}
#endif
