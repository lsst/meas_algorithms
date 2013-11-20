// -*- lsst-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

//
// We need this macro so as to avoid having commas in the 2nd argument to SWIG_SHARED_PTR_DERIVED,
// which confuses the swig parser.  It's also convenient
//
%define %MASKEDIMAGE(PIXTYPE)
lsst::afw::image::MaskedImage<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
%enddef
%define %EXPOSURE(PIXTYPE)
    lsst::afw::image::Exposure<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
%enddef

//
// Must go Before the %include
//
%define %PsfCandidatePtr(TYPE)
%shared_ptr(lsst::meas::algorithms::PsfCandidate<TYPE>);
%enddef
//
// Must go After the %include
//
%define %PsfCandidate(NAME, TYPE)
%template(PsfCandidate##NAME) lsst::meas::algorithms::PsfCandidate<TYPE>;
%template(makePsfCandidate) lsst::meas::algorithms::makePsfCandidate<TYPE>;

//
// When swig sees a SpatialCellImageCandidates it doesn't know about PsfCandidates; all it knows is that it has
// a SpatialCellImageCandidate, and SpatialCellCandidates don't know about e.g. getSource().
//
// We therefore provide a cast to PsfCandidate<TYPE> and swig can go from there;  In fact,
// we can cast all the way from the ultimate base class, so let's do that.
//
%inline %{
    PTR(lsst::meas::algorithms::PsfCandidate<TYPE>)
        cast_PsfCandidate##NAME(PTR(lsst::afw::math::SpatialCellCandidate) candidate) {
        return boost::dynamic_pointer_cast<lsst::meas::algorithms::PsfCandidate<TYPE> >(candidate);
    }
%}

%enddef

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


%ignore PsfFactoryBase;

%include "lsst/meas/algorithms/PSF.h"
%include "lsst/meas/algorithms/PsfCandidate.h"
%include "lsst/meas/algorithms/SpatialModelPsf.h"


%PsfCandidatePtr(float);
%PsfCandidate(F, float);

 //
// N.b. Swig won't will be able to resolve the overload for *FromPsfCandidates
// if you define another image type (there are no dependent parameters); so you'll have to
// append some type marker (e.g. "I") to the name
//
%template(pair_Psf_vector_double) std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, std::vector<double> >;
%template(pair_vector_double_KernelList) std::pair<std::vector<double>, lsst::afw::math::KernelList>;
%template(pair_bool_double) std::pair<bool, double>;
%template(pair_Kernel_double_double) std::pair<lsst::afw::math::Kernel::Ptr, std::pair<double, double> >;

%template(createKernelFromPsfCandidates) lsst::meas::algorithms::createKernelFromPsfCandidates<float>;
%template(fitSpatialKernelFromPsfCandidates) lsst::meas::algorithms::fitSpatialKernelFromPsfCandidates<float>;
%template(countPsfCandidates) lsst::meas::algorithms::countPsfCandidates<float>;
%template(subtractPsf) lsst::meas::algorithms::subtractPsf<%MASKEDIMAGE(float)>;
%template(fitKernelParamsToImage) lsst::meas::algorithms::fitKernelParamsToImage<%MASKEDIMAGE(float)>;
%template(fitKernelToImage) lsst::meas::algorithms::fitKernelToImage<%MASKEDIMAGE(float)>;

%{
#include "lsst/meas/algorithms/SingleGaussianPsf.h"
#include "lsst/meas/algorithms/PcaPsf.h"
%}

%import "lsst/afw/table/io/ioLib.i"

%declareTablePersistable(ImagePsf, lsst::meas::algorithms::ImagePsf);
%declareTablePersistable(KernelPsf, lsst::meas::algorithms::KernelPsf);
%declareTablePersistable(SingleGaussianPsf, lsst::meas::algorithms::SingleGaussianPsf);
%declareTablePersistable(DoubleGaussianPsf, lsst::meas::algorithms::DoubleGaussianPsf);
%declareTablePersistable(PcaPsf, lsst::meas::algorithms::PcaPsf);

%include "lsst/meas/algorithms/ImagePsf.h"
%include "lsst/meas/algorithms/KernelPsf.h"
%include "lsst/meas/algorithms/SingleGaussianPsf.h"
%include "lsst/meas/algorithms/DoubleGaussianPsf.h"
%include "lsst/meas/algorithms/PcaPsf.h"

%lsst_persistable(lsst::meas::algorithms::ImagePsf);
%lsst_persistable(lsst::meas::algorithms::KernelPsf);
%lsst_persistable(lsst::meas::algorithms::SingleGaussianPsf);
%lsst_persistable(lsst::meas::algorithms::DoubleGaussianPsf);
%lsst_persistable(lsst::meas::algorithms::PcaPsf);

%include "lsst/meas/algorithms/WarpedPsf.i"
