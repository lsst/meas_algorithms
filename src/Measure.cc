/// \file

#include "lsst/pex/logging/Trace.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/SecondMoments.h"
#include "lsst/meas/algorithms/Ellipse.h"

namespace image = lsst::afw::image;
namespace detection = lsst::afw::detection;
namespace algorithms = lsst::meas::algorithms;

/************************************************************************************************************/
/**
 * \brief Set some fields in a Source from foot (which was found in mimage)
 */
template<typename MaskedImageT>
void algorithms::measureSource(lsst::afw::detection::Source::Ptr src,    ///< the Source to receive results
                               MaskedImageT& mimage,      ///< image wherein Footprint dwells
                               lsst::afw::detection::Footprint const& foot, ///< Footprint to measure
                               float background                ///< background level to subtract
                              ) {
    FootprintCentroid<MaskedImageT> centroid(mimage);
    centroid.apply(foot);
    FootprintSecondMoments<MaskedImageT> secondMoments(mimage, centroid.getX(), centroid.getY());
    secondMoments.apply(foot);
    
    Ellipse ellipse(centroid.getX(), 
    				centroid.getY(),
    				secondMoments.getXx(), 
    				secondMoments.getYy(), 
    				secondMoments.getXy());
    
    
    //set some properties of the source
    //src->setXFlux(centroid.getX());
    //src->setYFlux(centroid.getY());
    //src->setFwhmA(ellipse.getMajorAxis());
    //src->setFwhmB(ellipse.getMinorAxis());
    //src->setFwhmTheta(ellipse.getTheta());
}

//
// Explicit instantiations
//
template void algorithms::measureSource(detection::Source::Ptr src, image::MaskedImage<float>& mimage,
                                        detection::Footprint const &foot, float background);
//
// Why do we need double images?
//    
#if 1
template void algorithms::measureSource(detection::Source::Ptr src, image::MaskedImage<double>& mimage,
                                        detection::Footprint const &foot, float background);
#endif
