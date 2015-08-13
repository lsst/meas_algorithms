// -*- lsst-C++ -*-

%{
#include "lsst/meas/algorithms/CoaddBoundedField.h"
%}
%shared_ptr(lsst::afw::math::BoundedField)
%shared_ptr(lsst::afw::image::Wcs)
%declareTablePersistable(CoaddBoundedField, lsst::meas::algorithms::CoaddBoundedField);

%ignore std::vector<lsst::meas::algorithms::CoaddBoundedFieldElement>::vector(size_type);
%ignore std::vector<lsst::meas::algorithms::CoaddBoundedFieldElement>::resize(size_type);
%template(CoaddBoundedFieldElementVector) std::vector<lsst::meas::algorithms::CoaddBoundedFieldElement>;

%include "lsst/meas/algorithms/CoaddBoundedField.h"

%pythoncode %{
CoaddBoundedField.Element = CoaddBoundedFieldElement
CoaddBoundedField.ElementVector = CoaddBoundedFieldElementVector
%}

%castShared(lsst::meas::algorithms::CoaddBoundedField, lsst::afw::math::BoundedField)
