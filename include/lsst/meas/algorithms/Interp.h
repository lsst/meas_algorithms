// -*- LSST-C++ -*-

/*
 * LSST Data Management System
 * Copyright 2008-2015 LSST Corporation.
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

#if !defined(LSST_MEAS_ALGORITHMS_INTERP_H)
#define LSST_MEAS_ALGORITHMS_INTERP_H
//!
// Interpolate over defects in a MaskedImage
//
#include <limits>
#include <vector>

#include "lsst/geom/Box.h"
#include "lsst/afw/image/Defect.h"
#include "lsst/afw/image/MaskedImage.h"

namespace lsst {
namespace afw {
namespace detection {

class Psf;
}
}  // namespace afw

namespace meas {
namespace algorithms {
namespace interp {
/**
 * LPC coefficients for sigma = 1, S/N = infty
 */
double const lpc_1_c1 = 0.7737;
double const lpc_1_c2 = -0.2737;
/**
 * LPC coefficients for sigma = 1/sqrt(2), S/N = infty. These are the coeffs
 * to use when interpolating at 45degrees to the row/column
 */
double const lpc_1s2_c1 = 0.7358;
double const lpc_1s2_c2 = -0.2358;
/*
 * Used to debias min(x, y)
 */
double const min2GaussianBias = -0.5641895835;  ///< Mean value of the minimum of two N(0,1) variates

template <typename MaskedImageT>
std::pair<bool, typename MaskedImageT::Image::Pixel> singlePixel(int x, int y, MaskedImageT const &image,
                                                                 bool horizontal, double minval);
}  // namespace interp

/**
 * @brief Encapsulate information about a bad portion of a detector
 */
class Defect : public afw::image::DefectBase {
public:
    typedef std::shared_ptr<Defect> Ptr;  //!< shared pointer to Defect

    enum DefectPosition {
        LEFT = 1,         //!< defect is at left boundary
        NEAR_LEFT,        //!< defect is near left boundary
        WIDE_LEFT,        //!< defect is wide at left boundary
        MIDDLE,           //!< defect is in middle of frame
        WIDE_NEAR_LEFT,   //!< defect is near left, and wide
        WIDE,             //!< defect is in middle, and wide
        WIDE_NEAR_RIGHT,  //!< defect is near right, and wide
        NEAR_RIGHT,       //!< defect is near right boundary
        WIDE_RIGHT,       //!< defect is wide at right boundary
        RIGHT             //!< defect is at right boundary
    };

    enum { WIDE_DEFECT = 11 };  //!< minimum width of a WIDE defect

    explicit Defect(const geom::BoxI &bbox = geom::BoxI()  //!< Region's bounding box
                    )
            : afw::image::DefectBase(bbox), _pos(static_cast<DefectPosition>(0)), _type(0) {}
    virtual ~Defect() {}

    void classify(DefectPosition pos,  //!< Position of defect in chip
                  unsigned int type    //!< Type of defect
    ) {
        _pos = pos;
        _type = type;
    }

    unsigned int getType() const { return _type; }  //!< Return the defect's interpolation type
    DefectPosition getPos() const { return _pos; }  //!< Return the position of the defect
private:
    DefectPosition _pos;  //!< Position of defect
    unsigned int _type;   //!< Type of defect
};

template <typename MaskedImageT>
void interpolateOverDefects(MaskedImageT &image, afw::detection::Psf const &psf,
                            std::vector<Defect::Ptr> &badList, double fallbackValue = 0.0,
                            bool useFallbackValueAtEdge = false);

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif
