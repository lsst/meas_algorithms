// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Shape.h"
%}

%include "lsst/meas/algorithms/Shape.h"

%template(measureShapeF) lsst::meas::algorithms::measureShape<lsst::afw::image::Image<float> >;
%template(createmeasureShape) lsst::meas::algorithms::createmeasureShape<lsst::afw::image::Image<float> >;
