#if !defined(LSST_DETECTION_PCAPSF_H)
#define LSST_DETECTION_PCAPSF_H
//!
// Describe an image's PSF
//
#include "lsst/base.h"
#include "lsst/afw/detection/Psf.h"

namespace lsst { namespace meas { namespace algorithms {
            
/*!
 * \brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class pcaPsf : public lsst::afw::detection::KernelPsf {
public:
    typedef PTR(pcaPsf) Ptr;
    typedef CONST_PTR(pcaPsf) ConstPtr;

    /**
     * @brief constructors for a pcaPsf
     *
     * Parameters:
     */
    explicit pcaPsf(PTR(Kernel) kernel);
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive&, unsigned int const) {
        boost::serialization::void_cast_register<pcaPsf, PSF>(
            static_cast<pcaPsf*>(0), static_cast<PSF*>(0));
    }
};

}}}

namespace boost {
namespace serialization {

template <class Archive>
inline void save_construct_data(
    Archive& ar, lsst::meas::algorithms::pcaPsf const* p,
    unsigned int const) {
    lsst::afw::math::Kernel const* kernel = p->getKernel().get();
    ar << make_nvp("kernel", kernel);
}

template <class Archive>
inline void load_construct_data(
    Archive& ar, lsst::meas::algorithms::pcaPsf* p,
    unsigned int const) {
    lsst::afw::math::Kernel* kernel;
    ar >> make_nvp("kernel", kernel);
    ::new(p) lsst::meas::algorithms::pcaPsf(PTR(kernel));
}

}}

#endif
