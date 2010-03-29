// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Shape.h"
%}

%include "lsst/meas/algorithms/Shape.h"

%createMeasureProperty(MeasureShape, F, lsst::afw::image::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>);
