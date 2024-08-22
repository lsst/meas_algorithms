// March 21, 2023: heliohypy (heliocentric hypothesis code for python)
// Implementation of make_tracklets for python.

#include "lsst/meas/algorithms/glintlib.h"
#include "cmath"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

void fill_struct(longpair & out, longpair const& in) {
    out.i1 = in.i1;
    out.i2 = in.i2;
}


void fill_struct(glint_trail & out, glint_trail const& in) {
    out.x  = in.x;
    out.y  = in.y;
    out.length  = in.length;
    out.PA  = in.PA;
    out.linrms  = in.linrms;
    out.eqrms  = in.eqrms;
    out.magmean  = in.magmean;
    out.magrms  = in.magrms;
    out.stepsize  = in.stepsize;
    out.qc1  = in.qc1;
    out.npt  = in.npt;
    out.flashnum  = in.flashnum;
};

void fill_struct(point3d_index & out, point3d_index const& in) {
    out.x  = in.x;
    out.y  = in.y;
    out.z  = in.z;
    out.index  = in.index;
}; 


template<typename T>
std::vector<T> ndarray_to_vec(py::array_t<T> py_vec) {
    std::vector<T> vec = {};

    // Get a reference to the py_array data
    auto data_ref = py_vec.unchecked();

    // Place the numpy data into the c++ type
    for (long int i = 0; i < data_ref.size(); i++) {
        T data_out;
        auto &data_in = data_ref[i];

        fill_struct(data_out, data_in);

        vec.push_back(data_out);
    }

    return vec;
}

template<typename T>
py::array vec_to_ndarray(std::vector<T> const& vec) {
    // Allocate a structured numpy array of type T
    auto py_vec = py::array_t<T>(vec.size());

    // Get a mutable reference to the ndarray data
    auto data_ref = py_vec.mutable_unchecked();

    // Place vector data into numpy array
    for (long int i = 0; i < data_ref.size(); i++) {
        auto data_in = vec[i];
        auto &data_out = data_ref[i];

        fill_struct(data_out, data_in);
    }
    return py_vec;
}


// findGlints: July 17, 2024:
// Minimalist wrapper, only handles the python <-> C++ translations.
// All the interesting algorithmic stuff happens in functions called from solarsyst_dyn_geo01.
// This is the pixel x,y version.

std::tuple<py::array, py::array>findGlints(
    FindGlintsConfig config,
    py::array_t<point3d_index> py_detvec
  ) {
  cout << "C++ wrapper for find_glints_xypix\n";
  
  std::vector <point3d_index> detvec = ndarray_to_vec(py_detvec);
  int status = 0;
  std::vector <glint_trail> trailvec;
  std::vector <longpair> trail2det;
     
  status = find_glints_xypix(detvec, config, trailvec, trail2det);
  if(status!=0) {
    cerr << "ERROR: find_glints_xypix returned failure status " << status << "\n";
    auto py_empty = vec_to_ndarray<point3d_index>({});
    return(std::make_tuple(py_empty, py_empty));
  }
      
  auto py_detout1 = vec_to_ndarray<glint_trail>(trailvec);
  auto py_detout2 = vec_to_ndarray<longpair>(trail2det);

  return(std::make_tuple(py_detout1, py_detout2));
}


template <typename S>
py::array_t<S> create_recarray(size_t n) {
    return py::array_t<S>(n);
}
#define NDARRAY_FACTORY(S) m.def("create_" #S, &create_recarray<S>);

PYBIND11_MODULE(glintpy, m) {
    m.doc() = "pybind11 I/O test"; // optional module docstring
    
    PYBIND11_NUMPY_DTYPE(longpair, i1, i2);
    PYBIND11_NUMPY_DTYPE(point3d_index, x, y, z, index);
    PYBIND11_NUMPY_DTYPE(glint_trail, x, y, length, PA, linrms, eqrms, magmean, magrms, stepsize, qc1, npt, flashnum);
		     
    NDARRAY_FACTORY(longpair)
    NDARRAY_FACTORY(point3d_index)
    NDARRAY_FACTORY(glint_trail)


    // Config class for FindGlints
    py::class_<FindGlintsConfig>(m, "FindGlintsConfig")
      .def(py::init<>())
      .def_readwrite("minpoints", &FindGlintsConfig::minpoints) 
      .def_readwrite("maxgcr", &FindGlintsConfig::maxgcr)
      .def_readwrite("maxrange", &FindGlintsConfig::maxrange)
      .def_readwrite("centerknown", &FindGlintsConfig::centerknown)
      .def_readwrite("incenRA", &FindGlintsConfig::incenRA)
      .def_readwrite("incenDec", &FindGlintsConfig::incenDec)
      .def_readwrite("freq_downscale", &FindGlintsConfig::freq_downscale)
      .def_readwrite("freq_upscale", &FindGlintsConfig::freq_upscale)
      .def_readwrite("max_phase_err", &FindGlintsConfig::max_phase_err);

    m.def("findGlints", &findGlints, "Identify glint trails produced by space junk, using pixel x,y coordinates");
}


