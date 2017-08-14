// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008-2011 LSST Corporation.
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
 
#include <iostream>
#include <cmath>

#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/afw/geom/Angle.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ImagePsf

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop
#include "boost/test/floating_point_comparison.hpp"

using namespace lsst::meas::algorithms;
using namespace lsst::afw::detection;
using namespace lsst::afw::image;
using namespace lsst::afw::geom;
using namespace lsst::afw::geom::ellipses;

class TestGaussianPsf : public lsst::meas::algorithms::ImagePsf {
public:

    TestGaussianPsf(int size, Quadrupole const & ellipse) :
        ImagePsf(true), _size(size), _ellipse(ellipse)
    {
        assert(size % 2);
    }

    TestGaussianPsf(int size, double radius) :
        ImagePsf(true), _size(size), _ellipse(1.0, 1.0, 0.0)
    {
        assert(size % 2);
        _ellipse.scale(radius);
    }

    PTR(Psf) clone() const {
        return std::make_shared<TestGaussianPsf>(*this);
    }

    PTR(Psf) resized(int width, int height) const {
        throw LSST_EXCEPT(lsst::pex::exceptions::LogicError, "Not Implemented");
    }

private:
    virtual PTR(Image) doComputeKernelImage(
        Point2D const & position, Color const & color
    ) const {
        PTR(Image) result(new Image(_size, _size));
        result->setXY0(-_size / 2, -_size / 2);
        Eigen::Matrix2d cInv = _ellipse.getMatrix().inverse();
        double norm = std::sqrt(cInv.determinant()) / (2.0 * M_PI);
        for (int row = 0; row < _size; ++row) {
            for (int col = 0; col < _size; ++col) {
                Eigen::Vector2d xy(col + result->getX0(), row + result->getY0());
                double z = norm * std::exp(-0.5 * xy.dot(cInv * xy));
                (*result)(col, row) = z;
            }
        }
        return result;
    }

    virtual Box2I doComputeBBox(
        Point2D const & position, Color const & color
    ) const {
        return Box2I(Point2I(-_size/2, -_size/2), Extent2I(_size, _size));
    }

    int _size;
    Quadrupole _ellipse;
};

void checkShape(int size, Quadrupole const & ellipse, double tol) {
    TestGaussianPsf psf(size, ellipse);
    Quadrupole shape = psf.computeShape();
    BOOST_CHECK_CLOSE(ellipse.getIxx(), shape.getIxx(), tol);
    BOOST_CHECK_CLOSE(ellipse.getIyy(), shape.getIyy(), tol);
    BOOST_CHECK_CLOSE(ellipse.getIxy(), shape.getIxy(), tol);
}

void checkApertureFlux(int size, double sigma, double radius, double tol) {
    TestGaussianPsf psf(size, sigma);
    double flux = psf.computeApertureFlux(radius);
    double check = 1.0 - std::exp(-0.5*(radius*radius)/(sigma*sigma));
    BOOST_CHECK_CLOSE(flux, check, tol);
}

// n.b. tolerance is in %, so values are closer than they appear

BOOST_AUTO_TEST_CASE(PsfShape) {
    checkShape(25, Quadrupole(7.0, 7.0, 0.0), 1E-2);
    checkShape(25, Quadrupole(7.0, 5.0, 1.0), 1E-2);
    checkShape(25, Quadrupole(5.0, 7.0, -1.0), 1E-2);
}

BOOST_AUTO_TEST_CASE(PsfApertureFlux) {
    checkApertureFlux(25, 5.0, 3.0, 1E-2);
    checkApertureFlux(25, 5.0, 5.0, 1E-2);
    checkApertureFlux(25, 5.0, 10.0, 1E-2);
}
