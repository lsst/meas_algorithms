// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Photometry.h"
%}

%include "lsst/meas/algorithms/Photometry.h"

%define %declarePhotometry(PIXTYPE, SUFFIX)
    %template(MeasurePhotometry ## SUFFIX) lsst::meas::algorithms::MeasurePhotometry<lsst::afw::image::MaskedImage<PIXTYPE> >;
    %template(createMeasurePhotometry) lsst::meas::algorithms::createMeasurePhotometry<lsst::afw::image::MaskedImage<PIXTYPE> >;
%enddef

%declarePhotometry(float, F)
 // %declarePhotometry(double, D)
