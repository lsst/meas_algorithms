// -*- LSST-C++ -*-
#ifndef LSST_MEAS_ALGORITHMS_SHAPELETPSF_H
#define LSST_MEAS_ALGORITHMS_SHAPELETPSF_H

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

/**
 * @file
 *
 * @brief A PSF class that describes the PSF in terms of its shapelet decomposition.
 *
 * @author Mike Jarvis
 */
#include <vector>

#include "lsst/afw/table/Source.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/algorithms/PsfCandidate.h"
#include "lsst/meas/algorithms/ShapeletKernel.h"
//#include "lsst/meas/algorithms/PSF.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletPsfImpl;

    class ShapeletPsf : public afw::detection::Psf
    {
    public:
        typedef afw::detection::Psf Base;

        typedef float PixelT;
        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::image::Exposure<PixelT> Exposure;
        typedef lsst::meas::algorithms::PsfCandidate<PixelT>::PtrList PsfCandidateList;
        typedef Exposure::MaskedImageT::Image Image;
        typedef lsst::afw::math::SpatialCellSet SpatialCellSet;
        typedef lsst::afw::geom::Point2D Point;
        typedef lsst::afw::geom::Extent2I Extent;
        typedef lsst::afw::image::Color Color;
        typedef lsst::afw::math::Kernel Kernel;

        typedef boost::shared_ptr<ShapeletPsf> Ptr;
        typedef boost::shared_ptr<const ShapeletPsf> ConstPtr;

        /*!
         * @brief Construct the ShapeletPsf from a list of Sources
         *
         * Lupton and Jarvis discussed having the constructor take an optional Filter argument.
         * I think in this case, the filter info can be grabbed from the star candidates in psfCandidateList.
         * However, if this doesn't work for some reason, then we should add such an argument here.
         */
        ShapeletPsf(
            const Exposure& exposure,   ///< the exposure on which to measure the decomposition
            const PsfCandidateList& psfCandidateList,   ///< List of PSF candidates
            const Policy& policy        ///< a Policy constructed from ShapeletPsfDeterminer.ConfigClass
        );

        /*!
         * @brief Destructor needs to delete pImpl
         */
        ~ShapeletPsf();

        /*!
         * @brief Copy constructor does a shallow copy.
         *
         * The copy shares the details with the original, so all changes to 
         * either one affect the other.
         */
        ShapeletPsf(const ShapeletPsf& rhs);

        /*!
         * @brief Make a clone of this
         */
        virtual Base::Ptr clone() const { 
            return boost::make_shared<ShapeletPsf>(*this); 
        }

        /*!
         * @brief Get the cellSet of Psf candidates used for the interpolation.
         *
         * This isn't the same as the cellSet that was input in the constructor, because:
         *
         * 1) The candidates now have measured Shapelets
         * 2) Some candidates are marked bad if they were deemed to be outliers.
         */
        const SpatialCellSet& getCellSet() const;

// there is no implementation so I'm commenting it out for now -- RO 2011-02-14
//         /*!
//          * @brief Get the preferred size of the image.
//          *
//          * If you are flexible about the size of the image version of the
//          * kernel, then in routines like getKernel (in the base class Psf), you
//          * can omit the size argument, and just let the Psf class pick out a
//          * good size to use.  The way the base class implements this is to call
//          * this function to get a good size to use.
//          * FIXME: This doesn't seem to be the case right now.  Should I get rid
//          * of this method, or should be implement this as a virtual method
//          * in the base class?
//          *
//          * In our case, the size is chosen to be 10 sigma in each direction,
//          * where sigma is the shapelet scale used for measuring the Psf (which
//          * is designed to be optimal for the average star in the input cellSet).
//          */
//         Extent getDefaultExtent() const;

        /*!
         * @brief Get a general Kernel that varies across an image.
         *
         * [ Required virtual function from base class Psf ]
         *
         * The color term is currently ignored, but is provided for future 
         * implementation
         *
         * The width and height are optional.
         * If not specified, then the Kernel will choose an appropriate size 
         * automatically.
         */
        Kernel::ConstPtr doGetKernel(
            const Color& color    ///< color to interpolate to
        ) const;

        Kernel::Ptr doGetKernel(
            const Color& color    ///< color to interpolate to
        );

        /*! 
         * Use the base class implementation of doComputeImage.
         * It might be slightly inefficient because it does color 
         * interpolation first.  Then computes the image of the 
         * kernel at a particular point.
         *
         * Probably faster to interpolate on both color and position.
         * Then compute the image from the local kernel.
         *
         * But there are weird things with the normalize parameter
         * that I think are poor design, so I don't particularly want to 
         * duplicate them here.  Specifically, the normalizePeak
         * argument optionally sets the center value of the image to 1.
         *
         * First, it is erroneously (or perhaps presumptuously) called the 
         * peak, whereas not all kernels will have the peak at the center
         * (e.g. comatic PSF).
         *
         * Second, this is a different convention than was proposed in 
         * Kernel where the bool argument optionally sets the _sum_ to 1.
         *
         * I think any normalization operation should be moved to being 
         * a method of the Image class, and let the caller subsequently call
         * that method rather than use the confusing boolean argument here.
         * i.e.
         * image = computeImage(color, pos, size);
         * image.normalizeCenterToUnity();
         *
         * Much clearer than:
         * image = computeImage(color, pos, size, true) 
         *   // true here means normalize center pixel to unity
         */
#if 0
        virtual Image::Ptr doComputeImage(
            const Color& color,       ///< color to interpolate to
            const Point& ccdXY,       ///< position to interpolate to
            const Extent& size,       ///< width/height of Kernel image
            bool normalizePeak        ///< normalize the image (so center is 1)?
        ) const;
#endif


    private :
        ShapeletPsfImpl* pImpl;

        // Not implemented.
        ShapeletPsf& operator=(const ShapeletPsf& rhs);

    };

}}}

#endif
