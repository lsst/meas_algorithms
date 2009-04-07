// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Shape.h"
%}

%include "lsst/meas/algorithms/Shape.h"

%template(MeasureShapeF) lsst::meas::algorithms::MeasureShape<lsst::afw::image::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> >;
%template(createMeasureShape) lsst::meas::algorithms::createMeasureShape<lsst::afw::image::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> >;
