// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Photometry.h"
%}

%include "lsst/meas/algorithms/Photometry.h"

%define %declarePhotometry(PIXTYPE, SUFFIX)
    %createMeasureProperty(MeasurePhotometry, SUFFIX, lsst::afw::image::MaskedImage<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>);
%enddef

%declarePhotometry(float, F)
 // %declarePhotometry(double, D)
