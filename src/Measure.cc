#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Measure.h"

namespace image = lsst::afw::image;
namespace detection = lsst::afw::detection;
namespace algorithms = lsst::meas::algorithms;

template<typename MaskedImageT>
algorithms::Measure<MaskedImageT>::Measure(MaskedImageT &img) 
    : lsst::daf::data::LsstBase(typeid(this)) {
    _img = img;
}

template<typename MaskedImageT>
void algorithms::Measure<MaskedImageT>::measureSource(detection::Source::Ptr pDia,
                                                      detection::Footprint const &fp,
                                                      float background
                                                     ) {
    FootprintCentroid<MaskedImageT> centroid(fp, _img);
    centroid.apply();

    float xCentroid = _img.getX0() + centroid.getX();
    float yCentroid = _img.getY0() + centroid.getY();
    
    pDia->setColc(xCentroid);
    pDia->setRowc(yCentroid);
    pDia->setFlux(centroid.getSum());
}

template<typename MaskedImageT>
void algorithms::Measure<MaskedImageT>::measureSource(detection::Source::Ptr pDia,
                                                      detection::Footprint::Ptr fpPtr,
                                                      float background
                                                     ) {
    measureSource(pDia, *fpPtr, background);
}

//
// Explicit instantiations
//
template class algorithms::Measure<image::MaskedImage<float> >;
//
// Why do we need double images?
//
#if 1
template class algorithms::Measure<image::MaskedImage<double> >;
#endif
