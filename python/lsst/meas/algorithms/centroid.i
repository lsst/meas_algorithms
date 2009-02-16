// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Centroid.h"
%}

%include "lsst/meas/algorithms/Centroid.h"

%template(measureCentroidF) lsst::meas::algorithms::measureCentroid<lsst::afw::image::Image<float> >;
%template(createmeasureCentroid) lsst::meas::algorithms::createmeasureCentroid<lsst::afw::image::Image<float> >;
