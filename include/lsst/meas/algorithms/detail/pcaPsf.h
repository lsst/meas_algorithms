#if !defined(LSST_DETECTION_PCAPSF_H)
#define LSST_DETECTION_PCAPSF_H
//!
// Describe an image's PSF
//
#include "lsst/meas/algorithms/PSF.h"
#include "PsfImpl.h"

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * \brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class pcaPsf : public PSF {
public:
    typedef boost::shared_ptr<pcaPsf> Ptr;

    /**
     * @brief constructors for a pcaPsf
     *
     * Parameters:
     */
    explicit pcaPsf(lsst::afw::math::Kernel::PtrT kernel);

    lsst::afw::image::Image<PSF::PixelT>::Ptr getImage(double const x, double const y) const;
private:
    double doGetValue(double const dx, double const dy, int const xPositionInImage, int const yPositionInImage) const;

private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, unsigned int const version) {
        boost::serialization::void_cast_register<pcaPsf, PSF>(
            static_cast<pcaPsf*>(0), static_cast<PSF*>(0));
    };
};

}}}

namespace boost {
namespace serialization {

using boost::serialization::make_nvp;

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::pcaPsf const* p,
    unsigned int const file_version) {
    boost::shared_ptr<const lsst::afw::math::Kernel> kernel = p->getKernel();
    ar << make_nvp("kernel", kernel);
};

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::pcaPsf* p,
    unsigned int const file_version) {
    lsst::afw::math::Kernel::PtrT kernel;
    ar >> make_nvp("kernel", kernel);
    ::new(p) lsst::meas::algorithms::pcaPsf(kernel);
};

}}

#endif
