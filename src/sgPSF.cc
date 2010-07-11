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
 
/*!
 * Represent a PSF as a circularly symmetrical single Gaussian
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <cmath>
#include "lsst/pex/exceptions.h"
#include "lsst/meas/algorithms/detail/sgPsf.h"
#include "lsst/afw/image/ImageUtils.h"

namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/**
 * Constructor for a sgPsf
 */
sgPsf::sgPsf(int width,                         ///< Number of columns in realisations of PSF
             int height,                        ///< Number of rows in realisations of PSF
             double sigma,                       ///< Width of Gaussian
             double,        ///< needed to match PSF.h definition but unused
             double         ///< needed to match PSF.h definition but unused
            ) :
    PSF(width, height),
    _sigma(sigma) {

    if (_sigma <= 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::DomainErrorException,
                          (boost::format("sigma may not be 0: %g") % _sigma).str());
    }
    
    if (width > 0) {
        lsst::afw::math::GaussianFunction1<double> sg(_sigma);
        setKernel(lsst::afw::math::SeparableKernel::Ptr(new lsst::afw::math::SeparableKernel(width, height, sg, sg)));
    }
}

/// \brief Evaluate the PSF at (dx, dy) (relative to the centre), taking the central amplitude to be 1.0
double sgPsf::doGetValue(double const dx,            ///< Desired column (relative to centre of PSF)
                         double const dy,            ///< Desired row (relative to centre of PSF)
                         int,                        ///< Desired column position in image (think "CCD")
                         int                         ///< Desired row position in image (think "CCD")
                        ) const {
    double const r2 = dx*dx + dy*dy;
    double const psf = exp(-r2/(2*_sigma*_sigma));
    
    return psf;
}

/*
 * Return an Image of the the PSF at the point (x, y), setting the sum of all the PSF's pixels to 1.0
 *
 * The specified position is a floating point number, and the resulting image will
 * have a PSF with the correct fractional position, with the centre within pixel (width/2, height/2)
 * Specifically, fractional positions in [0, 0.5] will appear above/to the right of the center,
 * and fractional positions in (0.5, 1] will appear below/to the left (0.9999 is almost back at middle)
 */
lsst::afw::image::Image<PSF::Pixel>::Ptr sgPsf::getImage(double const x, ///< column posn in parent %image
                                                         double const y  ///< row posn in parent %image
                                                        ) const {
    PSF::Image::Ptr image(new PSF::Image(getWidth(), getHeight()));

    // "ir" : (integer, residual)
    std::pair<int, double> const ir_dx = lsst::afw::image::positionToIndex(x, true);
    std::pair<int, double> const ir_dy = lsst::afw::image::positionToIndex(y, true);

    image->setXY0(ir_dx.first - getWidth()/2  + (ir_dx.second <= 0.5 ? 0 : 1),
                  ir_dy.first - getHeight()/2 + (ir_dy.second <= 0.5 ? 0 : 1));

    double const dx = ir_dx.second;     // fractional part of position
    double const dy = ir_dy.second;

    int const xcen = static_cast<int>(getWidth()/2);
    int const ycen = static_cast<int>(getHeight()/2);

    double sum = 0;
    for (int iy = 0; iy != image->getHeight(); ++iy) {
        PSF::Image::x_iterator row = image->row_begin(iy);
        for (int ix = 0; ix != image->getWidth(); ++ix) {
            PSF::Pixel val = getValue(ix - dx - xcen, iy - dy - ycen);

            row[ix] = val;
            sum += val;
        }
    }

    *image /= sum;

    return image;                                                    
}

//
// We need to make an instance here so as to register it
//
// \cond
namespace {
    volatile bool isInstance =
        PSF::registerMe<sgPsf, boost::tuple<int, int, double, double, double> >("SingleGaussian");
}

// \endcond
}}}
