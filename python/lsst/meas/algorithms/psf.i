// -*- lsst-C++ -*-

SWIG_SHARED_PTR_DERIVED(PSFPtrT, lsst::daf::data::LsstBase, lsst::meas::algorithms::PSF);
//
// Must go Before the %include
//
%define %PsfCandidatePtr(NAME, PIXEL_TYPES...)
SWIG_SHARED_PTR_DERIVED(PsfCandidate##NAME,
                        lsst::afw::math::SpatialCellCandidate,
                        lsst::meas::algorithms::PsfCandidate<lsst::afw::image::MaskedImage<PIXEL_TYPES> >);
%enddef
//
// Must go After the %include
//
%define %PsfCandidate(NAME, PIXEL_TYPES...)
%template(PsfCandidate##NAME) lsst::meas::algorithms::PsfCandidate<lsst::afw::image::MaskedImage<PIXEL_TYPES> >;
%template(makePsfCandidate) lsst::meas::algorithms::makePsfCandidate<lsst::afw::image::MaskedImage<PIXEL_TYPES> >;
%enddef

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%PsfCandidatePtr(F, float,  lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel);

%include "lsst/meas/algorithms/PSF.h"
%include "lsst/meas/algorithms/SpatialModelPsf.h"

%PsfCandidate(F, float,  lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel);
