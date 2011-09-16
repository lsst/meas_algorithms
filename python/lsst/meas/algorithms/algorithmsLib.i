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
#   include "lsst/pex/logging/Log.h"
#   include "lsst/pex/logging/BlockTimingLog.h"    
#   include "lsst/pex/logging/ScreenLog.h"
#   include "lsst/pex/logging/DualLog.h"
#   include "lsst/pex/logging/Debug.h"
#   include "lsst/afw.h"
#   include "lsst/afw/detection/Psf.h"
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

#   define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_ALGORITHMS_NUMPY_ARRAY_API
#   include "numpy/arrayobject.h"
#   include "lsst/ndarray/python.h"
#   include "lsst/ndarray/python/eigen.h"
%}

%inline %{
namespace boost { namespace filesystem { } }
namespace lsst { namespace afw {
        namespace detection { }
        namespace image { }
        namespace math { }
} }
namespace lsst { namespace meas { namespace algorithms { namespace interp {} } } }
namespace lsst { namespace daf { namespace data { } } }
    
using namespace lsst;
using namespace lsst::afw::image;
using namespace lsst::afw::detection;
using namespace lsst::meas::algorithms;
using namespace lsst::meas::algorithms::interp;
using namespace lsst::daf::data;
%}

%ignore boost::noncopyable;
namespace boost {
    class noncopyable {};
}

%init %{
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"                  // PTR(); should be in p_lsstSwig.i
%include "lsst/daf/base/persistenceMacros.i"

%lsst_exceptions();

%import "lsst/daf/data/dataLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"

%include "lsst/afw/image/lsstImageTypes.i"     // Image/Mask types and typedefs

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

SWIG_SHARED_PTR(ShapeletPtrT, lsst::meas::algorithms::Shapelet)
SWIG_SHARED_PTR(ShapeletInterpolationPtrT, lsst::meas::algorithms::ShapeletInterpolation)
SWIG_SHARED_PTR_DERIVED(LocalShapeletKernelPtrT, lsst::afw::math::AnalyticKernel,
    lsst::meas::algorithms::LocalShapeletKernel);
SWIG_SHARED_PTR_DERIVED(ShapeletKernelPtrT, lsst::afw::math::AnalyticKernel,
    lsst::meas::algorithms::ShapeletKernel);
SWIG_SHARED_PTR_DERIVED(ShapeletPsfCandidateT, lsst::afw::math::SpatialCellCandidate,
   lsst::meas::algorithms::ShapeletPsfCandidate);
SWIG_SHARED_PTR_DERIVED(ShapeletPsfPtrT, lsst::afw::detection::Psf, lsst::meas::algorithms::ShapeletPsf);
SWIG_SHARED_PTR(PsfCandidateListF,
    std::vector<lsst::meas::algorithms::SizeMagnitudeStarSelector::PsfCandidateList>);

%include "lsst/meas/algorithms/Shapelet.h" // causes tons of numpy warnings; due to Eigen?
%include "lsst/meas/algorithms/ShapeletInterpolation.h"
%include "lsst/meas/algorithms/ShapeletKernel.h"
%include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
%include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
%include "lsst/meas/algorithms/ShapeletPsf.h"

/************************************************************************************************************/

SWIG_SHARED_PTR(DefectPtrT, lsst::meas::algorithms::Defect);
SWIG_SHARED_PTR(DefectListT,  std::vector<lsst::meas::algorithms::Defect::Ptr>);

%include "lsst/meas/algorithms/Interp.h"

/************************************************************************************************************/
//
// We need this macro so as to avoid having commas in the 2nd argument to SWIG_SHARED_PTR_DERIVED,
// which confuses the swig parser.
//
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

%define %MeasureQuantityPtrs(SUFFIX, MEASUREMENT, PIXTYPE)
    SWIG_SHARED_PTR(MeasureQuantity##MEASUREMENT##SUFFIX, %MeasureQuantity##MEASUREMENT(PIXTYPE));
    SWIG_SHARED_PTR_DERIVED(Measure##MEASUREMENT##SUFFIX,
                            %MeasureQuantity##MEASUREMENT(PIXTYPE),
                            lsst::meas::algorithms::Measure##MEASUREMENT<%Exposure(PIXTYPE)>
        )
%enddef

%define %MeasureSources(SUFFIX, PIXTYPE)
    SWIG_SHARED_PTR(MeasureSources##SUFFIX,
                    lsst::meas::algorithms::MeasureSources<%Exposure(PIXTYPE)>);

    %MeasureQuantityPtrs(SUFFIX, Astrometry, PIXTYPE);
    %MeasureQuantityPtrs(SUFFIX, Photometry, PIXTYPE);
    %MeasureQuantityPtrs(SUFFIX, Shape, PIXTYPE);
%enddef

%MeasureSources(F, float);
%MeasureSources(I, int);

%include "lsst/meas/algorithms/MeasureQuantity.h"
%include "lsst/meas/algorithms/Measure.h"

/************************************************************************************************************/
/*
 * Now %template declarations
 */
%define %MeasureAlgorithm(SUFFIX, MEASUREMENT, PIXTYPE)
    %template(MeasureQuantity##MEASUREMENT##SUFFIX) %MeasureQuantity##MEASUREMENT(PIXTYPE);
    %template(Measure##MEASUREMENT##SUFFIX) lsst::meas::algorithms::Measure##MEASUREMENT<%Exposure(PIXTYPE)>;
    %template(makeMeasure##MEASUREMENT) lsst::meas::algorithms::makeMeasure##MEASUREMENT<%Exposure(PIXTYPE)>;
%enddef

%define %instantiate_templates(SUFFIX, PIXTYPE, UTILITIES)
#if UTILITIES
    %template(findCosmicRays) findCosmicRays<lsst::afw::image::MaskedImage<PIXTYPE,
                                                                           lsst::afw::image::MaskPixel> >;
    %template(interpolateOverDefects) interpolateOverDefects<lsst::afw::image::MaskedImage<PIXTYPE> >;
#endif

    %template(MeasureSources##SUFFIX) lsst::meas::algorithms::MeasureSources<%Exposure(PIXTYPE)>;
    %template(makeMeasureSources) lsst::meas::algorithms::makeMeasureSources<%Exposure(PIXTYPE)>;

    %MeasureAlgorithm(SUFFIX, Astrometry, PIXTYPE);
    %MeasureAlgorithm(SUFFIX, Photometry, PIXTYPE);
    %MeasureAlgorithm(SUFFIX, Shape, PIXTYPE);
%enddef

%instantiate_templates(F, float, 1)
%instantiate_templates(I, int, 0)

%include "lsst/meas/algorithms/detail/SincPhotometry.h";
%template(getCoeffImage) lsst::meas::algorithms::detail::getCoeffImage<float>;
%rename(computeGaussLeakage) lsst::meas::algorithms::detail::computeGaussLeakage;

%template(DefectListT) std::vector<lsst::meas::algorithms::Defect::Ptr>;

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
