// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Centroid.h"
%}

%include "lsst/meas/algorithms/Centroid.h"

%template(measureCentroidF) lsst::meas::algorithms::measureCentroid<lsst::afw::image::Image<float> >;
%template(createMeasureCentroid) lsst::meas::algorithms::createMeasureCentroid<lsst::afw::image::Image<float> >;
