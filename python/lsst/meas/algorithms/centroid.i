// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Centroid.h"
%}

%include "lsst/meas/algorithms/Centroid.h"

%define %createMeasureProperty(WHAT, TYPEID, IMAGE_T)
    %template() lsst::meas::algorithms::MeasureProperty< lsst::meas::algorithms::WHAT< IMAGE_T >, IMAGE_T >;
    %template(_##WHAT##TYPEID) lsst::meas::algorithms::WHAT<IMAGE_T >;
    %template(create##WHAT) lsst::meas::algorithms::create##WHAT<IMAGE_T >;
%enddef
 /*
  * Because of the version of createMeasureCentroid that doesn't take an image, we get swig warnings (#302)
  * about ambiguities.  The float create function is created, which is good for backward compatibility
  */
%createMeasureProperty(MeasureCentroid, F, lsst::afw::image::Image<float>);
%createMeasureProperty(MeasureCentroid, I, lsst::afw::image::Image<int>);
