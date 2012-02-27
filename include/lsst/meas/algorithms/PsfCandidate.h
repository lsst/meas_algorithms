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

#include "boost/make_shared.hpp"

#include "lsst/pex/policy.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/FootprintSet.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/math/SpatialCell.h"

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
    template <typename PixelT>
    class PsfCandidate : public lsst::afw::math::SpatialCellImageCandidate<PixelT> {
        using lsst::afw::math::SpatialCellImageCandidate<PixelT>::_image;
    public:
        using lsst::afw::math::SpatialCellImageCandidate<PixelT>::getXCenter;
        using lsst::afw::math::SpatialCellImageCandidate<PixelT>::getYCenter;
        using lsst::afw::math::SpatialCellImageCandidate<PixelT>::getWidth;
        using lsst::afw::math::SpatialCellImageCandidate<PixelT>::getHeight;
    
        typedef boost::shared_ptr<PsfCandidate<PixelT> > Ptr;
        typedef boost::shared_ptr<const PsfCandidate<PixelT> > ConstPtr;
        typedef std::vector<Ptr > PtrList;

        typedef lsst::afw::image::MaskedImage<PixelT,lsst::afw::image::MaskPixel,
                                              lsst::afw::image::VariancePixel> MaskedImageT;

        typedef lsst::afw::image::MaskPixel MaskPixel;

        /**
         * Construct a PsfCandidate from a specified source and image.
         *
         * The x/yCenter is set to source.getX/YAstrom()
         */
        PsfCandidate(PTR(afw::table::SourceRecord) const& source, ///< The detected Source
                     CONST_PTR(lsst::afw::image::Exposure<PixelT,lsst::afw::image::MaskPixel,
                               lsst::afw::image::VariancePixel>) parentExposure ///< The image wherein lie the Sources
        ) :
            lsst::afw::math::SpatialCellImageCandidate<PixelT>(source->getX(), source->getY()),
            _parentExposure(parentExposure),
            _offsetImage(),
            _undistImage(),
            _undistOffsetImage(),
            _source(source),
            _distortion(),
            _detector(),
            _haveImage(false),
            _haveUndistImage(false),
            _haveUndistOffsetImage(false),
            _amplitude(0.0), _var(1.0) {

            _stashDistortion();
        }
        
        /**
         * Construct a PsfCandidate from a specified source, image and xyCenter.
         */
        PsfCandidate(PTR(afw::table::SourceRecord) const& source, ///< The detected Source
                     CONST_PTR(lsst::afw::image::Exposure<PixelT,lsst::afw::image::MaskPixel,
                               lsst::afw::image::VariancePixel>) parentExposure, ///< The image wherein lie the Sources
                     double xCenter,    ///< the desired x center
                     double yCenter     ///< the desired y center
                    ) :
            lsst::afw::math::SpatialCellImageCandidate<PixelT>(xCenter, yCenter),
            _parentExposure(parentExposure),
            _offsetImage(),
            _undistImage(),
            _undistOffsetImage(),
            _source(source),
            _distortion(),
            _detector(),
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
        double getCandidateRating() const { return _source->getPsfFlux(); }
        
        /// Return the original Source
        PTR(afw::table::SourceRecord) getSource() const { return _source; }
        
        /// Return the best-fit amplitude
        double getAmplitude() const { return _amplitude; }
    
        /// Set the best-fit amplitude
        void setAmplitude(double amplitude) { _amplitude = amplitude; }

        /// Return the variance in use when fitting this object
        double getVar() const { return _var; }
    
        /// Set the variance to use when fitting this object
        void setVar(double var) { _var = var; }
    
        CONST_PTR(lsst::afw::image::MaskedImage<PixelT,
                  lsst::afw::image::MaskPixel,lsst::afw::image::VariancePixel>) getImage() const;
        CONST_PTR(lsst::afw::image::MaskedImage<PixelT,lsst::afw::image::MaskPixel,
                  lsst::afw::image::VariancePixel>) getImage(int width, int height) const;
        PTR(lsst::afw::image::MaskedImage<PixelT,
            lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel>) getOffsetImage(std::string const algorithm,
                                                             unsigned int buffer) const;
        PTR(lsst::afw::image::MaskedImage<PixelT,lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel>) getUndistOffsetImage(std::string const algorithm,
                                                                   unsigned int buffer,
                                                                   bool keepEdge=false) const;
        PTR(lsst::afw::image::MaskedImage<PixelT,lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel>) getUndistImage(int width, int height) const;
        PTR(lsst::afw::image::MaskedImage<PixelT,lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel>) getUndistImage() const;

        /// Return the number of pixels being ignored around the candidate image's edge
        static int getBorderWidth() { return _border; }
    
        /// Set the number of pixels to ignore around the candidate image's edge
        static void setBorderWidth(int border) { _border = border; }
    private:
        CONST_PTR(lsst::afw::image::Exposure<PixelT,lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel>) _parentExposure; // the %image that the Sources are found in

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
        
        PTR(afw::image::MaskedImage<PixelT,afw::image::MaskPixel,afw::image::VariancePixel>)
        offsetImage(
            PTR(afw::image::MaskedImage<PixelT,afw::image::MaskPixel,afw::image::VariancePixel>) img,
            std::string const algorithm,
            unsigned int buffer
        );
        
        PTR(afw::image::MaskedImage<PixelT,afw::image::MaskPixel,afw::image::VariancePixel>)
        extractImage(unsigned int width, unsigned int height) const;

        PTR(afw::image::MaskedImage<PixelT,afw::image::MaskPixel,afw::image::VariancePixel>) mutable _offsetImage; // %image offset to put center on a pixel
        PTR(afw::image::MaskedImage<PixelT,afw::image::MaskPixel,afw::image::VariancePixel>) mutable _undistImage; // %image undistort
        PTR(afw::image::MaskedImage<PixelT,afw::image::MaskPixel,afw::image::VariancePixel>) mutable _undistOffsetImage; // %image undistorted and offset
        PTR(afw::table::SourceRecord) _source; // the Source itself
        afw::cameraGeom::Distortion::ConstPtr _distortion;
        afw::cameraGeom::Detector::ConstPtr _detector;

        bool mutable _haveDetector;
        bool mutable _haveDistortion;
        bool mutable _haveImage;                    // do we have an Image to return?
        bool mutable _haveUndistImage;
        bool mutable _haveUndistOffsetImage;
        double _amplitude;                          // best-fit amplitude of current PSF model
        double _var;                                // variance to use when fitting this candidate
        static int _border;                         // width of border of ignored pixels around _image
        afw::geom::Point2D _xyCenter;
        static int _defaultWidth;
    };
    
    /**
     * Return a PsfCandidate of the right sort
     *
     * Cf. std::make_pair
     */
    template <typename PixelT>
    boost::shared_ptr<PsfCandidate<PixelT> >
    makePsfCandidate(PTR(afw::table::SourceRecord) const& source, ///< The detected Source
                     PTR(afw::image::Exposure<PixelT>) image    ///< The image wherein lies the object
                    )
    {
        return boost::make_shared< PsfCandidate<PixelT> >(source, image);
    }
   
}}}

#endif
