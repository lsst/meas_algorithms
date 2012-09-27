// -*- lsst-C++ -*-

%{
#include "lsst/afw/math/Kernel.h"
#include "lsst/meas/algorithms/CoaddPsf.h"
%}
%shared_ptr(lsst::meas::algorithms::CoaddPsf);
%shared_ptr(lsst::meas::algorithms::CoaddPsfKernel);

%include "lsst/meas/algorithms/CoaddPsf.h"

//%lsst_persistable(lsst::afw::detection::Psf);
