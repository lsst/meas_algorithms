// -*- LSST-C++ -*-
#ifndef LSST_MEAS_ALGORITHMS_SHAPELETPSFCANDIDATE_H
#define LSST_MEAS_ALGORITHMS_SHAPELETPSFCANDIDATE_H
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
 * @brief A module for determining which objects are good PSF stars
 *
 * @author Mike Jarvis
 */

#include "lsst/afw/math/SpatialCell.h"
#include "lsst/afw/detection/Source.h"
#include "boost/shared_ptr.hpp"

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletPsfCandidate : 
        public lsst::afw::math::SpatialCellCandidate
    {
    public :
        typedef lsst::afw::math::SpatialCellCandidate base;
        typedef lsst::afw::detection::Source Source;

        typedef boost::shared_ptr<ShapeletPsfCandidate> Ptr;
        typedef boost::shared_ptr<const ShapeletPsfCandidate> ConstPtr;

        /*!
         * @brief Constructor takes position, size, and original source
         *
         * The object stores the position, and initial guess of the size
         * and takes whatever other information is required from 
         * the Source data.  Currently, only sky is used, but it 
         * seemed wise to have source in the constructor in case we
         * decide to store something else too.  For example, we
         * will eventually want to extract some kind of color 
         * information.
         *
         * I don't store the image, wcs, or weightImage in the 
         * candidate itself, since all candidates will have the 
         * same things for these.
         *
         * Once we start doing PCA-style interpolations that use
         * sources from different exposures, we will have to do
         * something different.  But that algorithm will require
         * a pretty significant re-thinking of the design anyway, 
         * so no point in worrying about that yet.
         *
         * (No destructor, copy constructor, or op=, since the defaults 
         * do the right thing.)
         *
         * @note FIXME: source should really be a Source::ConstPtr, 
         * but that is not currently defined in Source.
         */
        inline ShapeletPsfCandidate(
            double x,       ///< X position of candidate
            double y,       ///< Y position of candidate
            double size,    ///< Initial estimate of size
            Source::Ptr source   ///< Original source
        ) : 
            base(x,y), _size(size), _source(source), _rating(1.)
        {}

        /*!
         * @brief Set the shapelet decomposition
         */
        inline void setShapelet(Shapelet::ConstPtr shapelet)
        { 
            _shapelet = shapelet; 
            if (_shapelet->hasCovariance())
                _rating = _shapelet->getValues()(0) / sqrt((*_shapelet->getCovariance())(0,0));
        }

        /*!
         * @brief Get position
         */
        inline double getX() const { return base::getXCenter(); }
        inline double getY() const { return base::getYCenter(); }

        /*!
         * @brief Get size
         */
        inline double getSize() const { return _size; }

        /*!
         * @brief Get source
         */
        inline Source::Ptr getSource() const 
        { 
            assert(_source);
            return _source; 
        }

        /*!
         * @brief Get the shapelet decomposition
         */
        inline Shapelet::ConstPtr getShapelet() const 
        { 
            assert(_shapelet);
            return _shapelet; 
        }

        /*! 
         * @brief Check if shapelet decomposition is set
         */
        inline bool hasShapelet() const { return _shapelet; }

        /*!
         * @brief Define "goodness" of candidate for SpatialCell
         */
        inline double getCandidateRating() const { return _rating; }

        /*!
         * @brief Mark the candidate as BAD.
         *
         * SpatialCellCandidate::setStatus is a bit of a pain to use, since
         * there doesn't seem to be any way around explicitly specifying 
         * the full namespace and class specification of BAD.
         * So do it once here.
         */
        inline void setBad() 
        { 
            lsst::afw::math::SpatialCellCandidate::setStatus(
                lsst::afw::math::SpatialCellCandidate::BAD);
        }

    private :

        double _size;
        Source::Ptr _source;
        Shapelet::ConstPtr _shapelet;
        double _rating;
    };

}}}

#endif
