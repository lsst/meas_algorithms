#if !defined(LSST_DETECTION_SGPSF_H)
#define LSST_DETECTION_SGPSF_H
//!
// Describe an image's PSF
//
#include "lsst/meas/algorithms/PSF.h"
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/void_cast.hpp"

// Forward declarations

namespace lsst { namespace meas { namespace algorithms {
    class sgPsf;
}}}

namespace boost {
namespace serialization {
    template <class Archive>
    void save_construct_data(
        Archive& ar, lsst::meas::algorithms::sgPsf const* p,
        unsigned int const file_version);
}}

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * \brief Represent a PSF as a circularly symmetrical double Gaussian
 */
class sgPsf : public PSF {
public:
    typedef boost::shared_ptr<sgPsf> Ptr;

    /**
     * @brief constructors for a sgPsf
     *
     * Parameters:
     */
    explicit sgPsf(int width, int height, double sigma, double=0, double=0);

    lsst::afw::image::Image<PSF::Pixel>::Ptr getImage(double const x, double const y) const;
private:
    double doGetValue(double const dx, double const dy, int xPositionInImage, int yPositionInImage) const;

    double _sigma;                     ///< Width of Gaussian

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<sgPsf, PSF>(
            static_cast<sgPsf*>(0), static_cast<PSF*>(0));
    }
    template <class Archive>
    friend void boost::serialization::save_construct_data(
            Archive& ar, sgPsf const* p, unsigned int const file_version);
};

}}} // namespace lsst::meas::algorithms

namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::sgPsf const* p,
    unsigned int const) {
    int width = p->getKernel()->getWidth();
    int height = p->getKernel()->getHeight();
    ar << make_nvp("width", width);
    ar << make_nvp("height", height);
    ar << make_nvp("sigma", p->_sigma);
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::sgPsf* p,
    unsigned int const) {
    int width;
    int height;
    double sigma;
    ar >> make_nvp("width", width);
    ar >> make_nvp("height", height);
    ar >> make_nvp("sigma", sigma);
    ::new(p) lsst::meas::algorithms::sgPsf(width, height, sigma);
}

}}


#endif
