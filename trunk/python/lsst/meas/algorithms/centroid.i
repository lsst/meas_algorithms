// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Centroid.h"
%}

%include "lsst/meas/algorithms/Centroid.h"

%template(MeasureCentroidF) lsst::meas::algorithms::MeasureCentroid<lsst::afw::image::Image<float> >;
%template(createMeasureCentroid) lsst::meas::algorithms::createMeasureCentroid<lsst::afw::image::Image<float> >;
