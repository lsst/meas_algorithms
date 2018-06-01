// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_IMAGEPCA_H)
#define LSST_MEAS_ALGORITHMS_IMAGEPCA_H

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2012 LSST Corporation.
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

/**
 * @file
 *
 * @brief Class for doing PCA on PSF stars.
 *
 * @ingroup algorithms
 */
#include <memory>
#include <utility>
#include <vector>

#include "lsst/afw.h"

namespace lsst {
namespace meas {
namespace algorithms {

template <typename ImageT>
class PsfImagePca : public afw::image::ImagePca<ImageT> {
    typedef typename afw::image::ImagePca<ImageT> Super;  ///< Base class
public:
    /// Ctor
    explicit PsfImagePca(bool constantWeight = true, int border = 3)
            : Super(constantWeight), _border(border) {}

    /// Generate eigenimages that are normalised and background-subtracted
    ///
    /// The background subtraction ensures PSF variation doesn't couple with small background errors.
    virtual void analyze();

private:
    int const _border;  ///< Border width for background subtraction
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif
