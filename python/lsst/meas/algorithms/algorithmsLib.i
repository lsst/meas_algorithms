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

%include "lsst/meas/algorithms/algorithms_fwd.i"

//---------- Warning suppression ----------------------------------------------------------------------------

%{
#   pragma clang diagnostic ignored "-Warray-bounds" // PyTupleObject has an array declared as [1]
%}
// FIXME: should resolve these warnings some other way
#pragma SWIG nowarn=362                 // operator=  ignored

//---------- Dependencies that don't need to be seen by downstream imports ----------------------------------

%import "lsst/afw/geom/Box.i"
%import "lsst/afw/geom/ellipses/ellipses_fwd.i"
%import "lsst/afw/image/Defect.i"
%import "lsst/afw/image/Image.i"
%import "lsst/afw/image/MaskedImage.i"
%import "lsst/afw/image/Exposure.i"
%import "lsst/afw/math/spatialCell.i"
%import "lsst/afw/math/kernel.i"
%import "lsst/afw/table/Source.i"
%import "lsst/afw/detection/Psf.i"
%import "lsst/afw/detection/Footprint.i"

namespace lsst { namespace pex { namespace policy {
class Policy;
}}}
%shared_ptr(lsst::pex::policy::Policy);

%pythoncode %{
import lsst.afw.detection.detectionLib
%}

//---------- ndarray and Eigen NumPy conversion typemaps ----------------------------------------------------

%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_ALGORITHMS_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "ndarray.i"

%declareNumPyConverters(lsst::meas::algorithms::Shapelet::ShapeletVector)
%declareNumPyConverters(lsst::meas::algorithms::Shapelet::ShapeletCovariance)

//---------- meas_algorithms classes and functions ----------------------------------------------------------

%include "lsst/meas/algorithms/PsfCandidate.i"
%include "lsst/meas/algorithms/psfs.i"
%include "lsst/meas/algorithms/PsfAttributes.i"
%include "lsst/meas/algorithms/CR.i"
%include "lsst/meas/algorithms/Interp.i"
%include "lsst/meas/algorithms/Algorithm.i"
%include "lsst/meas/algorithms/CentroidControl.i"
%include "lsst/meas/algorithms/FluxControl.i"
%include "lsst/meas/algorithms/ShapeControl.i"
%include "lsst/meas/algorithms/miscAlgorithms.i"
%include "lsst/meas/algorithms/Measure.i"
%include "lsst/meas/algorithms/SincPhotometry.i"
%include "lsst/meas/algorithms/Shapelet.i"

/************************************************************************************************************/

%define %Exposure(PIXTYPE)
    lsst::afw::image::Exposure<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
%enddef

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
