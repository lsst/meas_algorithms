// -*- LSST-C++ -*-
#if !defined(LSST_MEAS_ALGORITHMS_PSFCANDIDATE_H)
#define LSST_MEAS_ALGORITHMS_PSFCANDIDATE_H

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
 * @brief Class used by SpatialCell for spatial PSF fittig
 *
 * @ingroup algorithms
 */
#include <vector>

#include "boost/shared_ptr.hpp"

#include "lsst/afw.h"
#include "lsst/pex/policy.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/math/SpatialCell.h"
//#include "lsst/afw/cameraGeom.h"

namespace lsst {
namespace meas {
namespace algorithms {
    /** 
     * @brief Class stored in SpatialCells for spatial Psf fitting
     * 
     * PsfCandidate is a detection that may turn out to be a PSF. We'll
     * assign them to sets of SpatialCells; these sets will then be used to fit
     * a spatial model to the PSF.
     */
    template <typename ExposureT>
    class PsfCandidate : public lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT> {
        using lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>::_image;
    public:
        using lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>::getXCenter;
        using lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>::getYCenter;
        using lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>::getWidth;
        using lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>::getHeight;
    
        typedef boost::shared_ptr<PsfCandidate<ExposureT> > Ptr;
        typedef boost::shared_ptr<const PsfCandidate<ExposureT> > ConstPtr;
        typedef std::vector<Ptr > PtrList;

        typedef typename ExposureT::MaskedImageT MaskedImageT;
        
        /**
         * Construct a PsfCandidate from a specified source and image.
         *
         * The x/yCenter is set to source.getX/YAstrom()
         */
        PsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                     typename ExposureT::ConstPtr parentExposure ///< The image wherein lie the Sources
                    ) :
            lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>(source.getXAstrom(), source.getYAstrom()),
            _parentExposure(parentExposure),
            _distortion(),
            _detector(),
            _offsetImage(),
            _undistImage(),
            _undistOffsetImage(),
            _source(source),
            _haveImage(false),
            _haveUndistImage(false),
            _haveUndistOffsetImage(false),
            _amplitude(0.0), _var(1.0) {

            _stashDistortion();
        }
        
        /**
         * Construct a PsfCandidate from a specified source, image and xyCenter.
         */
        PsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                     typename ExposureT::ConstPtr parentExposure, ///< The image wherein lie the Sources
                     double xCenter,    ///< the desired x center
                     double yCenter     ///< the desired y center
                    ) :
            lsst::afw::math::SpatialCellImageCandidate<typename ExposureT::MaskedImageT>(xCenter, yCenter),
            _parentExposure(parentExposure),
            _distortion(),
            _detector(),
            _offsetImage(),
            _undistImage(),
            _undistOffsetImage(),
            _source(source),
            _haveImage(false),
            _haveUndistImage(false),
            _haveUndistOffsetImage(false),
            _amplitude(0.0), _var(1.0) {

            _stashDistortion();
        }
        
        /// Destructor
        virtual ~PsfCandidate() {};
        
        /**
         * Return Cell rating
         * 
         * @note Required method for use by SpatialCell
         */
        double getCandidateRating() const { return _source.getPsfFlux(); }
        
        /// Return the original Source
        lsst::afw::detection::Source const& getSource() const { return _source; }
        
        /// Return the best-fit amplitude
        double getAmplitude() const { return _amplitude; }
    
        /// Set the best-fit amplitude
        void setAmplitude(double amplitude) { _amplitude = amplitude; }

        /// Return the variance in use when fitting this object
        double getVar() const { return _var; }
    
        /// Set the variance to use when fitting this object
        void setVar(double var) { _var = var; }
    
        typename ExposureT::MaskedImageT::ConstPtr getImage() const;
        typename ExposureT::MaskedImageT::Ptr getOffsetImage(std::string const algorithm, unsigned int buffer) const;
        typename ExposureT::MaskedImageT::Ptr getUndistOffsetImage(std::string const algorithm,
                                                       unsigned int buffer) const;
        typename ExposureT::MaskedImageT::Ptr getUndistImage(int width, int height) const;

        
        /// Return the number of pixels being ignored around the candidate image's edge
        static int getBorderWidth() { return _border; }
    
        /// Set the number of pixels to ignore around the candidate image's edge
        static void setBorderWidth(int border) { _border = border; }
    private:
        typename ExposureT::ConstPtr _parentExposure; // the %image that the Sources are found in

        void _stashDistortion() {
            _haveDetector = _haveDistortion = false;
            if (_parentExposure->getDetector()) {
                _detector = _parentExposure->getDetector();
                _haveDetector = true;                
            }
            if (_haveDetector && _parentExposure->getDetector()->getDistortion()) { 
                _distortion = _parentExposure->getDetector()->getDistortion();
                _haveDistortion = true;
            }
        }
        
        typename ExposureT::MaskedImageT::Ptr offsetImage(typename ExposureT::MaskedImageT::Ptr img,
                                                          std::string const algorithm, unsigned int buffer);
        typename ExposureT::MaskedImageT::Ptr extractImage(unsigned int width, unsigned int height) const;
        typename ExposureT::MaskedImageT::Ptr mutable _offsetImage; // %image offset to put center on a pixel
        typename ExposureT::MaskedImageT::Ptr mutable _undistImage; // %image undistort
        typename ExposureT::MaskedImageT::Ptr mutable _undistOffsetImage; // %image undistorted and offset
        lsst::afw::detection::Source const _source; // the Source itself
        lsst::afw::cameraGeom::Distortion::Ptr _distortion;
        lsst::afw::cameraGeom::Detector::Ptr _detector;

        bool mutable _haveDetector;
        bool mutable _haveDistortion;
        bool mutable _haveImage;                    // do we have an Image to return?
        bool mutable _haveUndistImage;
        bool mutable _haveUndistOffsetImage;
        double _amplitude;                          // best-fit amplitude of current PSF model
        double _var;                                // variance to use when fitting this candidate
        static int _border;                         // width of border of ignored pixels around _image
        lsst::afw::geom::Point2D _xyCenter;
        static int _defaultWidth;
    };
    
    /**
     * Return a PsfCandidate of the right sort
     *
     * Cf. std::make_pair
     *
     * @note It may be desirable to check that ImagePtrT really is a shared_ptr<image>. 
     * The code is written this way to allow the compiler to deduce the argument types.
     */
    template <typename T>
    struct PsfCandidate_traits {
        typedef T Image;
    };
    
    template <typename T>
    struct PsfCandidate_traits<boost::shared_ptr<T> > {
        typedef T Image;
    };
    
    template <typename ImagePtrT>
    boost::shared_ptr<PsfCandidate<typename PsfCandidate_traits<ImagePtrT>::Image> >
    makePsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                     ImagePtrT image                             ///< The image wherein lies the object
                    )
    {
        typedef typename PsfCandidate_traits<ImagePtrT>::Image Image;
        return typename PsfCandidate<Image>::Ptr(new PsfCandidate<Image>(source, image));
    }
   
}}}

#endif
