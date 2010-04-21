// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Centroid.h"
%}

%include "lsst/meas/algorithms/Centroid.h"

%createMeasureProperty(MeasureCentroid, F, lsst::afw::image::Image<float>);
%createMeasureProperty(MeasureCentroid, I, lsst::afw::image::Image<int>);
