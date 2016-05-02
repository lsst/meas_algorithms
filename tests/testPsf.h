// -*- LSST-C++ -*-
/*
 * LSST Data Management System
 * Copyright 2013 LSST Corporation.
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

#ifndef LSST_MEAS_ALGORITHMS_TESTS_testPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_TESTS_testPsf_h_INCLUDED

#include <boost/make_shared.hpp>
#include "lsst/pex/exceptions.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/image/Color.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/detection/Psf.h"

namespace test {
namespace foo {
namespace bar {

/// A Psf for testing positions involved in calculations
///
/// The TestPsf ensures that all positions on which calculations are based
/// are within the bounding box of an image.
class TestPsf : public lsst::afw::detection::Psf
{
public:
    template <typename ImageT>
    explicit TestPsf(CONST_PTR(ImageT) image, double value=1.0)
        : lsst::afw::detection::Psf(), _bbox(image->getBBox()), _value(value) {}

    virtual ~TestPsf() {}

    virtual PTR(lsst::afw::detection::Psf) clone() const {
        return PTR(TestPsf)(new TestPsf(_bbox, _value));
    }

    virtual PTR(Image) doComputeImage(lsst::afw::geom::Point2D const& position,
                                      lsst::afw::image::Color const&) const {
        assertPosition(position);
        return makeImage();
    }
    virtual PTR(Image) doComputeKernelImage(lsst::afw::geom::Point2D const& position,
                                            lsst::afw::image::Color const&) const {
        assertPosition(position);
        return makeImage();
    }
    virtual double doComputeApertureFlux(double radius,
                                         lsst::afw::geom::Point2D const& position,
                                         lsst::afw::image::Color const&) const {
        assertPosition(position);
        return _value;
    }
    virtual lsst::afw::geom::ellipses::Quadrupole doComputeShape(lsst::afw::geom::Point2D const& position,
                                                                 lsst::afw::image::Color const&) const {
        assertPosition(position);
        return lsst::afw::geom::ellipses::Quadrupole(0, 0, 0);
    }

protected:
    PTR(Image) makeImage() const {
        PTR(Image) image = std::make_shared<Image>(1, 1);
        *image = _value;
        return image;
    }

    void assertPosition(lsst::afw::geom::Point2D const& position) const {
        if (!_bbox.contains(lsst::afw::geom::Point2I(position))) {
            std::cout << "Position " << position << " outside bbox " << _bbox << std::endl;
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterError, "Position outside BBox");
        }
    }

    explicit TestPsf(lsst::afw::geom::Box2I const& bbox, double value=1.0) :
        lsst::afw::detection::Psf(), _bbox(bbox), _value(value) {}

private:

    lsst::afw::geom::Box2I _bbox;        ///< Bounds for PSF
    double _value;                       ///< Value of PSF
};

template <typename ImageT>
PTR(TestPsf) makeTestPsf(CONST_PTR(ImageT) image, double value=1.0) {
    return std::make_shared<TestPsf>(image, value);
}

}}} // namespace test::foo::bar

#endif
