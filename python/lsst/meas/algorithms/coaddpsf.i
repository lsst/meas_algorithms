// -*- lsst-C++ -*-

%{
#include "lsst/afw/math/Kernel.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
%}
%shared_ptr(lsst::meas::algorithms::CoaddPsf);

%include "lsst/meas/algorithms/CoaddPsf.h"
%import "lsst/afw/table/Exposure.i"
%lsst_persistable(lsst::meas::algorithms::CoaddPsf);
