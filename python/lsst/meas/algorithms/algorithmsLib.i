// -*- lsst-c++ -*-
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
#   include "lsst/pex/logging/ScreenLog.h"
#   include "lsst/pex/logging/DualLog.h"
#   include "lsst/pex/logging/Debug.h"
#   include "lsst/afw.h"
#   include "lsst/afw/geom.h"
#   include "lsst/meas/algorithms/CR.h"
#   include "lsst/meas/algorithms/Interp.h"
#   include "lsst/meas/algorithms/PSF.h"
#   include "lsst/meas/algorithms/SpatialModelPsf.h"
#   include "lsst/meas/algorithms/Measure.h"
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
%include "lsst/daf/base/persistenceMacros.i"

%lsst_exceptions();

%import "lsst/daf/data/dataLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"

%include "lsst/afw/image/lsstImageTypes.i"     // Image/Mask types and typedefs

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://taxelrod@lsstarchive.ncsa.uiuc.edu/DMS/meas/algorithms/trunk/python/lsst/meas/algorithms/algorithmsLib.i $"):
    """Return a version given a HeadURL string; default: afw's version"""
    return guessSvnVersion(HeadURL)
%}

/************************************************************************************************************/

%include "psf.i"
%include "lsst/meas/algorithms/CR.h"

%lsst_persistable(lsst::meas::algorithms::PSF);

/************************************************************************************************************/

SWIG_SHARED_PTR(DefectPtrT, lsst::meas::algorithms::Defect);
SWIG_SHARED_PTR(DefectListT,  std::vector<lsst::meas::algorithms::Defect::Ptr>);

%include "lsst/meas/algorithms/Interp.h"

/************************************************************************************************************/

SWIG_SHARED_PTR(MeasureSourcesF,
       lsst::meas::algorithms::MeasureSources<lsst::afw::image::Exposure<float, lsst::afw::image::MaskPixel, float> >);

%include "lsst/meas/algorithms/Measure.h"

/************************************************************************************************************/
/*
 * Now %template declarations
 */
%define %instantiate_templates(NAME, TYPE)
    %template(convolve) lsst::meas::algorithms::PSF::convolve<lsst::afw::image::Image<TYPE> >;
    %template(convolve) lsst::meas::algorithms::PSF::convolve<lsst::afw::image::MaskedImage<TYPE> >;
    %template(findCosmicRays) findCosmicRays<lsst::afw::image::MaskedImage<TYPE, lsst::afw::image::MaskPixel> >;
    %template(interpolateOverDefects) interpolateOverDefects<lsst::afw::image::MaskedImage<TYPE> >;
    %template(MeasureSources ## NAME)
        lsst::meas::algorithms::MeasureSources<lsst::afw::image::Exposure<TYPE, lsst::afw::image::MaskPixel, float> >;
    %template(makeMeasureSources) lsst::meas::algorithms::makeMeasureSources<lsst::afw::image::Exposure<TYPE> >;
%enddef

%instantiate_templates(F, float)

%template(DefectListT) std::vector<lsst::meas::algorithms::Defect::Ptr>;

/************************************************************************************************************/

%template(xyAndError) std::pair<double, double>;

%include "photometry.i"
%include "centroid.i"
%include "shape.i"

/******************************************************************************/
// Local Variables: ***
// eval: (setq indent-tabs-mode nil) ***
// End: ***
