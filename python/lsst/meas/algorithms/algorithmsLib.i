// -*- lsst-c++ -*-
%define meas_algorithmsLib_DOCSTRING
"
Python bindings for meas/algorithms module
"
%enddef

%feature("autodoc", "1");
%module(package="meas_algorithms",docstring=meas_algorithmsLib_DOCSTRING) algorithmsLib

// Suppress swig complaints
// I had trouble getting %warnfilter to work; hence the pragmas
 //#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#   include <exception>
#   include <list>
#   include <map>
#   include <boost/cstdint.hpp>
#   include <boost/shared_ptr.hpp>
#   include "lsst/afw/image.h"
#   include "lsst/afw/detection.h"
#   include "lsst/meas/algorithms/CR.h"
#   include "lsst/meas/algorithms/Interp.h"
#   include "lsst/meas/algorithms/PSF.h"
#   include "lsst/meas/algorithms/Measure.h"
%}

%inline %{
namespace boost { namespace filesystem { } }
namespace lsst { namespace afw {
        namespace detection { }
        namespace image { }
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
%include "lsst/daf/base/persistenceMacros.i"

%import "lsst/daf/data/dataLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"

%lsst_exceptions();

%include "lsst/afw/image/lsstImageTypes.i"     // Image/Mask types and typedefs

%pythoncode %{
def version(HeadURL = r"$HeadURL$"):
    """Return a version given a HeadURL string; default: afw's version"""
    return guessSvnVersion(HeadURL)

%}

%import "lsst/afw/math/mathLib.i"

/************************************************************************************************************/

SWIG_SHARED_PTR_DERIVED(PSFPtrT, lsst::daf::data::LsstBase, lsst::meas::algorithms::PSF);

%include "lsst/meas/algorithms/PSF.h"
%include "lsst/meas/algorithms/CR.h"

/************************************************************************************************************/

SWIG_SHARED_PTR(DefectPtrT, lsst::meas::algorithms::Defect);
SWIG_SHARED_PTR(DefectListT,  std::vector<lsst::meas::algorithms::Defect::Ptr>);

%include "lsst/meas/algorithms/Interp.h"

/************************************************************************************************************/

%include "lsst/meas/algorithms/Measure.h"

/************************************************************************************************************/
/*
 * Now %template declarations
 */
%define %instantiate_templates(TYPE)
    %template(convolve) lsst::meas::algorithms::PSF::convolve<lsst::afw::image::Image<TYPE> >;
    %template(findCosmicRays) findCosmicRays<lsst::afw::image::MaskedImage<TYPE, lsst::afw::image::MaskPixel> >;
    %template(interpolateOverDefects) interpolateOverDefects<lsst::afw::image::MaskedImage<TYPE> >;
    %template(measureSource) lsst::meas::algorithms::measureSource<lsst::afw::image::MaskedImage<TYPE> >;
%enddef

%instantiate_templates(float)
%instantiate_templates(double)

%template(DefectListT) std::vector<lsst::meas::algorithms::Defect::Ptr>;

/************************************************************************************************************/

%template(xyAndError) std::pair<double, double>;

%include "centroid.i"
%include "shape.i"

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
