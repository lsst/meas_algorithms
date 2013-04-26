// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

/*
 * The ImageMoments has a lot of interfaces that are tricky to expose to Python, so there are a lot
 * of semi-painful workarounds here:
 *
 *  - Results are returned via inner classes, but Swig can't handle inner classes, so we use typedefs
 *    here and "#ifndef SWIG" blocks in the .h file to trick it into thinking they aren't inner classes.
 *    We want them to remain inner classes on the C++ side so they can still have access to names
 *    defineded in the ImageMoments scope.
 *
 *  - Those same result structs contain Eigen members (some of which are in shared_ptrs) that aren't
 *    handled properly by the the typemaps in ndarray, so we provide custom wrappers for those.
 *
 *  - The convertRawMoments and correctWeightedMoments have optional Eigen pointer output arguments,
 *    used to compute derivatives, which are also not handled by the ndarray typemaps, so we change
 *    the signature to return a tuple of (originalResult, derivative) instead.
 */

%{
#include "lsst/meas/algorithms/ImageMoments.h"

// typedefs to help with tricking Swig about inner classes
namespace lsst { namespace meas { namespace algorithms {

typedef ImageMoments::EllipseResult ImageMomentsEllipseResult;
typedef ImageMoments::RawResult ImageMomentsRawResult;

}}} // namespace lsst::meas::algorithms

using lsst::meas::algorithms::ImageMoments;

// convenience function to return two Eigen objects as a pair of ndarrays; used below
template <typename VectorT, typename MatrixT>
    PyObject * makeVectorMatrixTuple(VectorT const & vector, MatrixT const & matrix) {
    PyObject * result = PyTuple_New(2);
    PyTuple_SET_ITEM(result, 0, ndarray::toPyObject(vector));
    PyTuple_SET_ITEM(result, 1, ndarray::toPyObject(matrix));
    return result;
}

%}

%declareNumPyConverters(lsst::meas::algorithms::ImageMoments::VectorQ)
%declareNumPyConverters(lsst::meas::algorithms::ImageMoments::VectorM)

namespace lsst { namespace meas { namespace algorithms {

class ImageMomentsEllipseResult;
class ImageMomentsRawResult;

}}} // namespace lsst::meas::algorithms

// n.b. We want to replace the C++ signatures for these methods with more Python-friendly ones
//      which means we have to put the %extend block before the %ignore commands, and put the
//      %ignore commands before the %include.
%extend lsst::meas::algorithms::ImageMoments {
    static PyObject * convertRawMoments(ImageMoments::VectorQ const & q) {
        ImageMoments::MatrixMQ dm_dq;
        ImageMoments::VectorM m = ImageMoments::convertRawMoments(q, &dm_dq);
        return makeVectorMatrixTuple(m, dm_dq);
    }
    static PyObject * correctWeightedMoments(
        lsst::afw::geom::ellipses::Quadrupole const & weight,
        ImageMoments::VectorM const & m
    ) {
        ImageMoments::MatrixM dc_dm;
        ImageMoments::VectorM c = ImageMoments::correctWeightedMoments(weight, m, &dc_dm);
        return makeVectorMatrixTuple(c, dc_dm);
    }
}

%ignore lsst::meas::algorithms::ImageMoments::convertRawMoments;
%ignore lsst::meas::algorithms::ImageMoments::correctWeightedMoments;

%include "lsst/meas/algorithms/ImageMoments.h"

// These %extend member functions can be used for both RawResults and 
%define %wrapEigenGetters()
    ndarray::Array<double,1,1> getMoments(PyObject ** PYTHON_SELF) {
        return ndarray::external(
            self->moments.data(),
            ndarray::makeVector(static_cast<int>(self->moments.size())),
            ndarray::ROW_MAJOR,
            ndarray::PyPtr(*PYTHON_SELF, true)
        );
    }
    PyObject * getCovariance() {
        if (self->covariance) {
            ndarray::Array<double,2,-2> r = ndarray::external(
                self->covariance->data(),
                ndarray::makeVector(
                    static_cast<int>(self->covariance->rows()),
                    static_cast<int>(self->covariance->cols())
                ),
                ndarray::COLUMN_MAJOR,
                self->covariance
            );
            return ndarray::toPyObject(r);
        } else {
            Py_RETURN_NONE;
        }
    }
%enddef

%extend lsst::meas::algorithms::ImageMomentsEllipseResult {
    %wrapEigenGetters()
    %pythoncode %{
        moments = property(getMoments)
        covariance = property(getCovariance)
        quadrupole = property(getQuadrupole)
        centroid = property(getCentroid)
        ellipse = property(getEllipse)
    %}
}

%extend lsst::meas::algorithms::ImageMomentsRawResult {
    %wrapEigenGetters()
    %pythoncode %{
        moments = property(getMoments)
        covariance = property(getCovariance)
    %}
}

%define %instantiateImageMoments(IM)
%template(measureRaw) lsst::meas::algorithms::ImageMoments::measureRaw< lsst::afw::image::IM >;
%template(measureEllipse) lsst::meas::algorithms::ImageMoments::measureEllipse< lsst::afw::image::IM >;
%enddef

%instantiateImageMoments(Image<float>);
%instantiateImageMoments(Image<double>);
%instantiateImageMoments(MaskedImage<float>);
%instantiateImageMoments(MaskedImage<double>);
