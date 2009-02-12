// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Shape.h"
%}

%include "lsst/meas/algorithms/Shape.h"

%template(ShapeFinderF) lsst::meas::algorithms::ShapeFinder<lsst::afw::image::Image<float> >;
%template(createShapeFinder) lsst::meas::algorithms::createShapeFinder<lsst::afw::image::Image<float> >;
