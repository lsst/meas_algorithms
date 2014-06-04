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
#   include "lsst/afw/geom/ellipses.h"
#   include "lsst/meas/algorithms.h"

#   define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_ALGORITHMS_NUMPY_ARRAY_API
#   include "numpy/arrayobject.h"
#   include "ndarray/swig.h"
#   include "ndarray/swig/eigen.h"

#ifdef __clang__
#pragma clang diagnostic ignored "-Warray-bounds"
#endif
%}

namespace lsst { namespace meas { namespace algorithms { namespace interp {} namespace photometry {} } } }

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

%include "ndarray.i"

/************************************************************************************************************/

%include "psf.i"
%include "coaddpsf.i"
%include "lsst/meas/algorithms/CR.h"

/************************************************************************************************************/

%declareNumPyConverters(lsst::meas::algorithms::Shapelet::ShapeletVector)
%declareNumPyConverters(lsst::meas::algorithms::Shapelet::ShapeletCovariance)

%shared_ptr(lsst::meas::algorithms::Shapelet)
%shared_ptr(lsst::meas::algorithms::ShapeletInterpolation)
%shared_ptr(lsst::meas::algorithms::LocalShapeletKernel);
%shared_ptr(lsst::meas::algorithms::ShapeletKernel);
%shared_ptr(lsst::meas::algorithms::ShapeletPsfCandidate);
%shared_vec(lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList);
%shared_ptr(std::vector<lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList>);

%include "lsst/meas/algorithms/Shapelet.h" // causes tons of numpy warnings; due to Eigen?
%include "lsst/meas/algorithms/ShapeletInterpolation.h"
%include "lsst/meas/algorithms/ShapeletKernel.h"
%include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
%include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"


%shared_ptr(lsst::meas::algorithms::Algorithm)
%shared_ptr(lsst::meas::algorithms::AlgorithmControl)

%include "lsst/meas/algorithms/Algorithm.h"

%shared_ptr(lsst::meas::algorithms::CentroidAlgorithm)
%shared_ptr(lsst::meas::algorithms::CentroidControl)
%shared_ptr(lsst::meas::algorithms::GaussianCentroidControl)
%shared_ptr(lsst::meas::algorithms::NaiveCentroidControl)
%shared_ptr(lsst::meas::algorithms::SdssCentroidControl)
%shared_ptr(lsst::meas::algorithms::RecordCentroidControl)

%shared_ptr(lsst::meas::algorithms::ShapeAlgorithm)
%shared_ptr(lsst::meas::algorithms::ShapeControl)
%shared_ptr(lsst::meas::algorithms::SdssShapeControl)

%shared_ptr(lsst::meas::algorithms::FluxAlgorithm)
%shared_ptr(lsst::meas::algorithms::FluxControl)
%shared_ptr(lsst::meas::algorithms::ApertureFluxControl)
%shared_ptr(lsst::meas::algorithms::EllipticalApertureFluxControl)
%shared_ptr(lsst::meas::algorithms::PeakLikelihoodFluxControl)
%shared_ptr(lsst::meas::algorithms::GaussianFluxControl)
%shared_ptr(lsst::meas::algorithms::NaiveFluxControl)
%shared_ptr(lsst::meas::algorithms::PsfFluxControl)
%shared_ptr(lsst::meas::algorithms::SincFluxControl)

%shared_ptr(lsst::meas::algorithms::CorrectFluxesControl)
%shared_ptr(lsst::meas::algorithms::ClassificationControl)
%shared_ptr(lsst::meas::algorithms::PixelFlagControl)
%shared_ptr(lsst::meas::algorithms::SkyCoordControl)
%shared_ptr(lsst::meas::algorithms::JacobianControl)
%shared_ptr(lsst::meas::algorithms::FocalPlaneControl)

%include "lsst/meas/algorithms/FluxControl.h"
%include "lsst/meas/algorithms/CentroidControl.h"
%include "lsst/meas/algorithms/ShapeControl.h"
%include "lsst/meas/algorithms/Classification.h"
%include "lsst/meas/algorithms/Jacobian.h"
%include "lsst/meas/algorithms/FocalPlane.h"
%include "lsst/meas/algorithms/PixelFlags.h"
%include "lsst/meas/algorithms/SkyCoord.h"
%include "lsst/meas/algorithms/RecordCentroid.h"
%include "lsst/meas/algorithms/SdssShapeControl.h"
%include "lsst/meas/algorithms/GaussianFluxControl.h"
%include "lsst/meas/algorithms/CorrectFluxes.h"
%include "lsst/meas/algorithms/Photometry.h"
%returnCopy(lsst::meas::algorithms::MeasureSources::getAlgorithms)
%returnSelf(lsst::meas::algorithms::MeasureSourcesBuilder::setCentroider)
%returnSelf(lsst::meas::algorithms::MeasureSourcesBuilder::addAlgorithm)
%include "lsst/meas/algorithms/Measure.h"

%extend lsst::meas::algorithms::MeasureSources {
%template(apply) apply<float>;
%template(apply) apply<double>;
%template(applyForced) applyForced<float>;
%template(applyForced) applyForced<double>;
%template(applyWithPeak) applyWithPeak<float>;
%template(applyWithPeak) applyWithPeak<double>;
%template(applyWithCoord) applyWithCoord<float>;
%template(applyWithCoord) applyWithCoord<double>;
%template(applyWithPixel) applyWithPixel<float>;
%template(applyWithPixel) applyWithPixel<double>;
}
%extend lsst::meas::algorithms::MeasureSourcesBuilder {
%pythoncode %{
    def addAlgorithms(self, iterable):
        for item in iterable:
            self.addAlgorithm(item)
        return self
%}
}

/************************************************************************************************************/

%shared_ptr(lsst::meas::algorithms::Defect);
%shared_vec(lsst::meas::algorithms::Defect::Ptr);
%shared_ptr(std::vector<lsst::meas::algorithms::Defect::Ptr>);

%include "lsst/meas/algorithms/Interp.h"

/************************************************************************************************************/

%define %Exposure(PIXTYPE)
    lsst::afw::image::Exposure<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
%enddef

/************************************************************************************************************/
/*
 * Now %template declarations
 */

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

%define %instantiate_templates(SUFFIX, PIXTYPE)
    %template(findCosmicRays) lsst::meas::algorithms::findCosmicRays<
                                  lsst::afw::image::MaskedImage<PIXTYPE,
                                                                lsst::afw::image::MaskPixel,
                                                                lsst::afw::image::VariancePixel> >;
    %template(interpolateOverDefects) lsst::meas::algorithms::interpolateOverDefects<
                                          lsst::afw::image::MaskedImage<PIXTYPE,
                                                                        lsst::afw::image::MaskPixel,
                                                                        lsst::afw::image::VariancePixel> >;
%enddef

%instantiate_templates(F, float)

%template(DefectListT) std::vector<lsst::meas::algorithms::Defect::Ptr>;

%init %{
    import_array();
%}

%template(SincCoeffsF) lsst::meas::algorithms::photometry::SincCoeffs<float>;
%template(SincCoeffsD) lsst::meas::algorithms::photometry::SincCoeffs<double>;
%template(calculateSincApertureFlux) lsst::meas::algorithms::photometry::calculateSincApertureFlux<
    lsst::afw::image::MaskedImage<float> >;


%shared_ptr(lsst::meas::algorithms::BinnedWcs)
%include "lsst/meas/algorithms/BinnedWcs.h"

%include "lsst/meas/algorithms/CoaddBoundedField.i"
