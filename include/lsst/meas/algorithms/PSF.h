// -*- LSST-C++ -*-

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
 
#if !defined(LSST_MEAS_ALGORITHMS_PSF_H)
#define LSST_MEAS_ALGORITHMS_PSF_H
//!
// Describe an image's PSF
//
#include <string>
#include "lsst/base.h"

namespace lsst {
namespace afw {
    namespace detection {
        class Psf;
    }
    namespace geom {
        template<typename T, int N> class Point;
        typedef Point<int,2> Point2I;
    }
    namespace image {
        template<typename T> class Image;
    }
}
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * @class PsfAttributes
 *
 * A class to contain various attributes of the Psf
 * - most notably, a width (1-D RMS size) to be used to
 *   make a single gaussian psf for fast convolution.
 */
class PsfAttributes {
public:
    enum Method { ADAPTIVE_MOMENT,      ///< Calculate width using adaptive Gaussian weights
                  FIRST_MOMENT,         ///< Calculate width using \<r>
                  SECOND_MOMENT,        ///< Calculate width using \<r^2>
                  NOISE_EQUIVALENT,     ///< Calculate width as sqrt(n_eff/(4 pi))
                  BICKERTON             ///< Weight \<r^2> by I^2 to avoid negative fluxes
    };

    PsfAttributes(CONST_PTR(lsst::afw::detection::Psf) psf, int const iX, int const iY);
    PsfAttributes(CONST_PTR(lsst::afw::detection::Psf) psf, lsst::afw::geom::Point2I const& cen);
    
    double computeGaussianWidth(Method how=ADAPTIVE_MOMENT);
    double computeEffectiveArea();
    
private:
    PTR(lsst::afw::image::Image<double>) _psfImage;
};

    
}}}
#endif
