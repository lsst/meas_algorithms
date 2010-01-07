// -*- lsst-c++ -*-
%define testLib_DOCSTRING
"
Various swigged-up C++ classes for testing
"
%enddef

%feature("autodoc", "1");
%module(package="testLib", docstring=testLib_DOCSTRING) testLib

%pythonnondynamic;
%naturalvar;  // use const reference typemaps

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%inline %{
#include "SillyCentroid.cc"
%}
