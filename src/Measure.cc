#include "lsst/pex/logging/Trace.h"
#include "lsst/detection/Measure.h"

namespace image = lsst::afw::image;
namespace afwDetection = lsst::afw::detection;

template<typename MaskedImageT>
lsst::detection::Measure<MaskedImageT>::Measure(MaskedImageT &img) 
    : lsst::daf::data::LsstBase(typeid(this)) {
    _img = img;
}

template<typename MaskedImageT>
void lsst::detection::Measure<MaskedImageT>::measureSource(afwDetection::Source::Ptr pDia,
                                                           afwDetection::Footprint const &fp,
                                                           float background
                                                    ) {
    float xCentroid = _img.getX0() /* + measureFunc.getXCentroid() */;
    float yCentroid = _img.getY0() /* + measureFunc.getYCentroid() */;
    float flux = 0;

    pDia->setColc(xCentroid);
    pDia->setRowc(yCentroid);
    pDia->setFlux(flux);
}

template<typename MaskedImageT>
void lsst::detection::Measure<MaskedImageT>::measureSource(afwDetection::Source::Ptr pDia,
                                                           afwDetection::Footprint::Ptr fpPtr,
                                                           float background
                                                          ) {
    measureSource(pDia, *fpPtr, background);
}

//
// Explicit instantiations
//
template class lsst::detection::Measure<image::MaskedImage<float> >;
//
// Why do we need double images?
//
#if 1
template class lsst::detection::Measure<image::MaskedImage<double> >;
#endif
