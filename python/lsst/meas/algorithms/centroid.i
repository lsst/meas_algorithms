// -*- lsst-c++ -*-

%{
#   include "lsst/meas/algorithms/Centroid.h"
%}

%template(xyAndError) std::pair<double, double>;

%include "lsst/meas/algorithms/Centroid.h"

%template(CentroiderF) lsst::meas::algorithms::Centroider<lsst::afw::image::Image<float> >;
%template(make_Centroider) lsst::meas::algorithms::make_Centroider<lsst::afw::image::Image<float> >;
