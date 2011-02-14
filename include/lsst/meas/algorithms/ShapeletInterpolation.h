// -*- LSST-C++ -*-
#ifndef LSST_MEAS_ALGORITHMS_SHAPELETINTERPOLATION_H
#define LSST_MEAS_ALGORITHMS_SHAPELETINTERPOLATION_H

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
 * @brief Defines the ShapeletInterpolation class
 *
 * @author Mike Jarvis
 */
#include "boost/shared_ptr.hpp"

#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "boost/shared_ptr.hpp"
#include "Eigen/Core"
#include <complex>

namespace lsst {
namespace meas {
namespace algorithms {

    class ShapeletInterpolationImpl;

    class ShapeletInterpolation 
    {
        /*!
         * @brief This class is an interpolator to get a Shapelet at an arbitrary position.
         *
         * In addition to providing an interface to an interpolation function
         * for Shapelets, it also measures the interpolating function.
         *
         * Thus, the constructor actually just loads some parameters from a
         * Policy file that are needed for the measurement and interpolation.
         *
         * The measurement and interpolation is then done with the 
         * calculate() function.
         * 
         * In practice, the ShapeletInterpolation object is made by 
         * ShapeletPsf, and is used by it to implement the interpolation.
         * I don't think there is any reason to make this object 
         * separately from a ShapeletPsf.
         */
    public:

        typedef boost::shared_ptr<ShapeletInterpolation> Ptr;
        typedef boost::shared_ptr<const ShapeletInterpolation> ConstPtr;

        typedef float PixelT;
        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::math::SpatialCellSet SpatialCellSet;
        typedef lsst::afw::image::Exposure<PixelT> Exposure;
        typedef lsst::afw::geom::PointD PointD;

        typedef lsst::meas::algorithms::Shapelet Shapelet;

        /*!
         * @brief Basic constructor just loads parameters from a policy file.
         *
         * (See the comment for the ShapeletPsf constructor for more details.)
         */
        ShapeletInterpolation(const Policy& policy);

        /*!
         * @brief Destructor needs to delete pImpl
         */
        ~ShapeletInterpolation() {};

        /*!
         * @brief Copy constructor does a shallow copy
         */
        ShapeletInterpolation(const ShapeletInterpolation& rhs);

        /*! 
         * @brief op= does a shallow copy
         */
        ShapeletInterpolation& operator=(const ShapeletInterpolation& rhs);

        /*!
         * @brief get the order of the shapelet
         */
        int getOrder() const;

        /*!
         * @brief get the order of the fit
         */
        int getFitOrder() const;

        /*!
         * @brief get the scale size of the shapelet
         */
        double getSigma() const;

        /*!
         * @brief the size of the shapelet vector
         */
        int getSize() const;

        /*!
         * @brief the number of fit coefficients
         */
        int getFitSize() const;

        /*!
         * @brief set a new value of sigma
         */
        void setSigma(double sigma);

        /*!
         * @brief Calculate the interpolation parameters from a SpatialCellSet.  
         *
         * The candidates must be ShapeletPsfCandidate's.
         */
        void calculate(
            SpatialCellSet::Ptr cellSet,    ///< the set of candidates
            const Exposure& exposure        ///< the exposure on which to measure the PSF
        );

        /*!
         * @brief Perform the interpolation
         */
        Shapelet::ConstPtr interpolate(const PointD& pos) const;
        Shapelet::ConstPtr interpolate(double x, double y) const;

        /*!
         * @brief Perform the interpolation for only one shapelet coefficient.
         */
        double interpolateSingleElement(const PointD& pos, int i) const;
        double interpolateSingleElement(double x, double y, int i) const;

    private :

        boost::shared_ptr<ShapeletInterpolationImpl> pImpl;
    };

}}}

#endif
