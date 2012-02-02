// -*- lsst-c++ -*-

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
 
%define meas_algorithmsLib_DOCSTRING
"
Python bindings for meas/algorithms module
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.algorithms",docstring=meas_algorithmsLib_DOCSTRING) algorithmsLib

// Suppress swig complaints
// I had trouble getting %warnfilter to work; hence the pragmas
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#   include <exception>
#   include <list>
#   include <map>
#   include <boost/cstdint.hpp>
#   include <boost/shared_ptr.hpp>
#   include "lsst/pex/logging.h"
#   include "lsst/pex/logging/BlockTimingLog.h"
#   include "lsst/pex/logging/DualLog.h"
#   include "lsst/pex/logging/ScreenLog.h"
#   include "lsst/afw.h"
#   include "lsst/afw/detection/Peak.h"
#   include "lsst/afw/detection/Psf.h"
#   include "lsst/afw/detection/AperturePhotometry.h"
#   include "lsst/meas/algorithms/Flags.h"
#   include "lsst/meas/algorithms/CR.h"
#   include "lsst/meas/algorithms/Interp.h"
#   include "lsst/meas/algorithms/PSF.h"
#   include "lsst/meas/algorithms/PsfCandidate.h"
#   include "lsst/meas/algorithms/SpatialModelPsf.h"
#   include "lsst/meas/algorithms/MeasureQuantity.h"
#   include "lsst/meas/algorithms/Measure.h"
#   include "lsst/meas/algorithms/Shapelet.h"
#   include "lsst/meas/algorithms/ShapeletInterpolation.h"
#   include "lsst/meas/algorithms/ShapeletKernel.h"
#   include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
#   include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
#   include "lsst/meas/algorithms/ShapeletPsf.h"
#   include "lsst/meas/algorithms/detail/SincPhotometry.h"
#   include "lsst/meas/algorithms/PhotometryControl.h"
#   include "lsst/meas/algorithms/AstrometryControl.h"
#   include "lsst/meas/algorithms/ShapeControl.h"

#   define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_ALGORITHMS_NUMPY_ARRAY_API
#   include "numpy/arrayobject.h"
#   include "lsst/ndarray/python.h"
#   include "lsst/ndarray/python/eigen.h"

#ifdef __clang__
#pragma clang diagnostic ignored "-Warray-bounds"
#endif
%}

namespace lsst { namespace meas { namespace algorithms { namespace interp {} } } }

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"                  // PTR(); should be in p_lsstSwig.i
%include "lsst/pex/config.h"            // LSST_CONTROL_FIELD.
%include "lsst/daf/base/persistenceMacros.i"

%lsst_exceptions();

%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"

%pythoncode %{
def version(HeadURL = r"$HeadURL$"):
    """Return a version given a HeadURL string; default: afw's version"""
    return guessSvnVersion(HeadURL)
%}

/************************************************************************************************************/

%include "psf.i"
%include "lsst/meas/algorithms/CR.h"

/************************************************************************************************************/

%declareNumPyConverters(lsst::meas::algorithms::Shapelet::ShapeletVector)
%declareNumPyConverters(lsst::meas::algorithms::Shapelet::ShapeletCovariance)

%shared_ptr(lsst::meas::algorithms::Shapelet)
%shared_ptr(lsst::meas::algorithms::ShapeletInterpolation)
%shared_ptr(lsst::meas::algorithms::LocalShapeletKernel);
%shared_ptr(lsst::meas::algorithms::ShapeletKernel);
%shared_ptr(lsst::meas::algorithms::ShapeletPsfCandidate);
%shared_ptr(lsst::meas::algorithms::ShapeletPsf);
%shared_vec(lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList);
%shared_ptr(std::vector<lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList>);

%include "lsst/meas/algorithms/Flags.h"
%include "lsst/meas/algorithms/Shapelet.h" // causes tons of numpy warnings; due to Eigen?
%include "lsst/meas/algorithms/ShapeletInterpolation.h"
%include "lsst/meas/algorithms/ShapeletKernel.h"
%include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
%include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
%include "lsst/meas/algorithms/ShapeletPsf.h"

%include "lsst/meas/algorithms/ExposurePatch.h"
%include "lsst/meas/algorithms/Algorithm.h"
%include "lsst/meas/algorithms/MeasureQuantity.h"
%include "lsst/meas/algorithms/Measure.h"

%extend lsst::meas::algorithms::MeasureQuantity {
%pythoncode %{
    def addAlgorithms(self, iterable):
        for item in iterable:
            self.addAlgorithm(item)
%}
}

%template(PhotometryControl) lsst::meas::algorithms::AlgorithmControl<lsst::afw::detection::Photometry>;
%template(AstrometryControl) lsst::meas::algorithms::AlgorithmControl<lsst::afw::detection::Astrometry>;
%template(ShapeControl) lsst::meas::algorithms::AlgorithmControl<lsst::afw::detection::Shape>;

%include "lsst/meas/algorithms/PhotometryControl.h"
%include "lsst/meas/algorithms/AstrometryControl.h"
%include "lsst/meas/algorithms/ShapeControl.h"

/************************************************************************************************************/

%shared_ptr(lsst::meas::algorithms::Defect);
%shared_vec(lsst::meas::algorithms::Defect::Ptr);
%shared_ptr(std::vector<lsst::meas::algorithms::Defect::Ptr>);

%include "lsst/meas/algorithms/Interp.h"

/************************************************************************************************************/

%define %Exposure(PIXTYPE)
    lsst::afw::image::Exposure<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
%enddef

%define %MeasureQuantity(MEASUREMENT, PIXTYPE)
    lsst::meas::algorithms::MeasureQuantity<lsst::afw::detection::MEASUREMENT, %Exposure(PIXTYPE)>
%enddef

%define %MeasureQuantityAstrometry(PIXTYPE)
    %MeasureQuantity(Astrometry, PIXTYPE)
%enddef
%define %MeasureQuantityPhotometry(PIXTYPE)
    %MeasureQuantity(Photometry, PIXTYPE)
%enddef
%define %MeasureQuantityShape(PIXTYPE)
    %MeasureQuantity(Shape, PIXTYPE)
%enddef

%define %MeasureQuantityPtrs(MEASUREMENT, PIXTYPE)
    %shared_ptr(%MeasureQuantity##MEASUREMENT(PIXTYPE));
    %shared_ptr(lsst::meas::algorithms::Measure##MEASUREMENT<%Exposure(PIXTYPE)>);
%enddef

%define %Algorithm(MEASUREMENT, PIXTYPE)
    lsst::meas::algorithms::Algorithm<lsst::afw::detection::MEASUREMENT, %Exposure(PIXTYPE)>
%enddef

%define %AlgorithmAstrometry(PIXTYPE)
    %Algorithm(Astrometry, PIXTYPE)
%enddef
%define %AlgorithmPhotometry(PIXTYPE)
    %Algorithm(Photometry, PIXTYPE)
%enddef
%define %AlgorithmShape(PIXTYPE)
    %Algorithm(Shape, PIXTYPE)
%enddef

%define %AlgorithmPtrs(MEASUREMENT, PIXTYPE)
    %shared_ptr(%Algorithm##MEASUREMENT(PIXTYPE))
%enddef

%define %MeasureSources(PIXTYPE)
    %shared_ptr(lsst::meas::algorithms::MeasureSources<%Exposure(PIXTYPE)>);

    %MeasureQuantityPtrs(Astrometry, PIXTYPE);
    %MeasureQuantityPtrs(Photometry, PIXTYPE);
    %MeasureQuantityPtrs(Shape, PIXTYPE);

    %AlgorithmPtrs(Astrometry, PIXTYPE);
    %AlgorithmPtrs(Photometry, PIXTYPE);
    %AlgorithmPtrs(Shape, PIXTYPE);

    %shared_ptr(lsst::meas::algorithms::ExposurePatch<%Exposure(PIXTYPE)>);
%enddef

%MeasureSources(float);
%MeasureSources(double);

%include "lsst/meas/algorithms/ExposurePatch.h"
%include "lsst/meas/algorithms/MeasureQuantity.h"
%include "lsst/meas/algorithms/Measure.h"
%include "lsst/meas/algorithms/Algorithm.h"

/************************************************************************************************************/
/*
 * Now %template declarations
 */
%define %MeasureAlgorithm(SUFFIX, MEASUREMENT, PIXTYPE)
    %template(MeasureQuantity##MEASUREMENT##SUFFIX) %MeasureQuantity##MEASUREMENT(PIXTYPE);
    %template(Measure##MEASUREMENT##SUFFIX) lsst::meas::algorithms::Measure##MEASUREMENT<%Exposure(PIXTYPE)>;
    %template(makeMeasure##MEASUREMENT) lsst::meas::algorithms::makeMeasure##MEASUREMENT<%Exposure(PIXTYPE)>;
    %template(Algorithm##MEASUREMENT##SUFFIX) %Algorithm##MEASUREMENT(PIXTYPE);
%enddef

%define %ExposurePatch(SUFFIX, PIXTYPE)
    %template(ExposurePatch##SUFFIX) lsst::meas::algorithms::ExposurePatch<%Exposure(PIXTYPE)>;
    %template(makeExposurePatch) lsst::meas::algorithms::makeExposurePatch<%Exposure(PIXTYPE)>;
    %template(ExposureList##SUFFIX) std::vector<CONST_PTR(%Exposure(PIXTYPE))>;
%template(ExposurePatchList##SUFFIX) std::vector<typename lsst::meas::algorithms::ExposurePatch<%Exposure(PIXTYPE)>::ConstPtr>;

%typemap(in) std::vector<CONST_PTR(%Exposure(PIXTYPE))> const {
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  size_t size = PySequence_Size($input);
  std::cout << "Converting sequence of " << size << std::endl;
  $1 = std::vector<CONST_PTR(%Exposure(PIXTYPE))>(size);
  for (i = 0; i < size; ++i) {
      PyObject* obj = PySequence_GetItem($input, i);
      CONST_PTR(%Exposure(PIXTYPE)) exp;
      if ((SWIG_ConvertPtr(obj, (void **) &exp, SWIGTYPE_p_Exposure##SUFFIX, 1)) == -1) return NULL;
      $1[i] = exp;
  }
}

%enddef

%define %instantiate_templates(SUFFIX, PIXTYPE, UTILITIES)
#if UTILITIES
    %template(findCosmicRays) lsst::meas::algorithms::findCosmicRays<
                                  lsst::afw::image::MaskedImage<PIXTYPE,
                                                                lsst::afw::image::MaskPixel,
                                                                lsst::afw::image::VariancePixel> >;
    %template(interpolateOverDefects) lsst::meas::algorithms::interpolateOverDefects<
                                          lsst::afw::image::MaskedImage<PIXTYPE,
                                                                        lsst::afw::image::MaskPixel,
                                                                        lsst::afw::image::VariancePixel> >;
#endif

    %template(MeasureSources##SUFFIX) lsst::meas::algorithms::MeasureSources<%Exposure(PIXTYPE)>;
    %template(makeMeasureSources) lsst::meas::algorithms::makeMeasureSources<%Exposure(PIXTYPE)>;

    %MeasureAlgorithm(SUFFIX, Astrometry, PIXTYPE);
    %MeasureAlgorithm(SUFFIX, Photometry, PIXTYPE);
    %MeasureAlgorithm(SUFFIX, Shape, PIXTYPE);

    %ExposurePatch(SUFFIX, PIXTYPE);

%enddef

%instantiate_templates(F, float, 1)
%instantiate_templates(D, double, 0)

%include "lsst/meas/algorithms/detail/SincPhotometry.h";
%template(getCoeffImage) lsst::meas::algorithms::detail::getCoeffImage<float>;
%rename(computeGaussLeakage) lsst::meas::algorithms::detail::computeGaussLeakage;

%template(DefectListT) std::vector<lsst::meas::algorithms::Defect::Ptr>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
