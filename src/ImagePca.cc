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
 * @brief Implementation of code to determine spatial model of PSF
 *
 * @file
 *
 * @ingroup algorithms
 */

#include "lsst/afw.h"
#include "lsst/meas/algorithms/ImagePca.h"

namespace lsst {
namespace meas {
namespace algorithms {

template <typename ImageT>
void PsfImagePca<ImageT>::analyze() {
    Super::analyze();

    typename Super::ImageList const &eImageList = this->getEigenImages();
    typename Super::ImageList::const_iterator iter = eImageList.begin(), end = eImageList.end();
    for (size_t i = 0; iter != end; ++i, ++iter) {
        PTR(ImageT) eImage = *iter;

        /*
         * Normalise eigenImages to have a maximum of 1.0.  For n > 0 they
         * (should) have mean == 0, so we can't use that to normalize
         */
        afw::math::Statistics stats = afw::math::makeStatistics(*eImage, (afw::math::MIN | afw::math::MAX));
        double const min = stats.getValue(afw::math::MIN);
        double const max = stats.getValue(afw::math::MAX);

        double const extreme = (fabs(min) > max) ? min : max;
        if (extreme != 0.0) {
            *eImage /= extreme;
        }

        /*
         * Estimate and subtract the mean background level from the i > 0
         * eigen images; if we don't do that then PSF variation can get mixed
         * with subtle variations in the background and potentially amplify
         * them disasterously.
         *
         * It is not at all clear that doing this is a good idea; it'd be
         * better to get the sky level right in the first place.
         */
        if (i > 0 && _border > 0) { /* not the zeroth KL component */
            int const height = eImage->getHeight();
            int const width = eImage->getWidth();
            double background;
            if (2 * _border >= std::min(height, width)) {
                // _Border consumes the entire image
                background = afw::math::makeStatistics(*afw::image::GetImage<ImageT>::getImage(eImage),
                                                       afw::math::MEDIAN)
                                     .getValue();
            } else {
                // Use the median of the edge pixels

                // If ImageT is a MaskedImage, unpack the Image
                std::shared_ptr<typename afw::image::GetImage<ImageT>::type> eImageIm =
                        afw::image::GetImage<ImageT>::getImage(eImage);

                int const nEdge = width * height - (width - 2 * _border) * (height - 2 * _border);
                std::vector<double> edgePixels(nEdge);

                std::vector<double>::iterator bi = edgePixels.begin();

                typedef typename afw::image::GetImage<ImageT>::type::x_iterator imIter;
                int y = 0;
                for (; y != _border; ++y) {  // Bottom border of eImage
                    for (imIter ptr = eImageIm->row_begin(y), end = eImageIm->row_end(y); ptr != end;
                         ++ptr, ++bi) {
                        *bi = *ptr;
                    }
                }
                for (; y != height - _border; ++y) {  // Left and right borders of eImage
                    for (imIter ptr = eImageIm->row_begin(y), end = eImageIm->x_at(_border, y); ptr != end;
                         ++ptr, ++bi) {
                        *bi = *ptr;
                    }
                    for (imIter ptr = eImageIm->x_at(width - _border, y), end = eImageIm->row_end(y);
                         ptr != end; ++ptr, ++bi) {
                        *bi = *ptr;
                    }
                }
                for (; y != height; ++y) {  // Top border of eImage
                    for (imIter ptr = eImageIm->row_begin(y), end = eImageIm->row_end(y); ptr != end;
                         ++ptr, ++bi) {
                        *bi = *ptr;
                    }
                }
                assert(std::distance(edgePixels.begin(), bi) == nEdge);

                background = afw::math::makeStatistics(edgePixels, afw::math::MEDIAN).getValue();
            }
            *eImage -= background;
        }
    }
}

#define INSTANTIATE_IMAGE(IMAGE) template class PsfImagePca<IMAGE>;

#define INSTANTIATE(TYPE)                       \
    INSTANTIATE_IMAGE(afw::image::Image<TYPE>); \
    INSTANTIATE_IMAGE(afw::image::MaskedImage<TYPE>);

INSTANTIATE(float);

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
