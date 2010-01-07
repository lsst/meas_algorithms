/*
 * Used In Most Places
 */
#include <cmath>
#include <vector>

#if !defined(NO_AFW)
#include "lsst/afw/image/Image.h"
#else
#include "Image.h"
#endif

namespace afwImage = lsst::afw::image;

struct FittedModel {
    enum { PEAK = 0, SKY, X0, Y0, SIGMA, NPARAM };
    
    enum {
        BAD_GUESS = -11,
        TOO_FEW = -12,
        CHI_SQUARED = -13,
        RANGE = -14,
        BAD_WIDTH = -15,
        LOST = -16,
        DIAGONAL = -17,
        BAD_A = -18,
        CONVERGE = 1,
        ITERATE = 2,
        ALMOST = 3,
        POOR = 4,
    };

    FittedModel(int status_, std::vector<double> params_, int iter_=0, double flamd_=0, double chnew_=0) :
        status(status_), params(params_), iter(iter_), flamd(flamd_), chnew(chnew_) { }
    int status;
    std::vector<double> params;
    int iter;
    double flamd;
    double chnew;
};

template<typename PixelT>
FittedModel twodg(afwImage::Image<PixelT> const& im, double x0, double y0);

