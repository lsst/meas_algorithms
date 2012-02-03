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

#include "lsst/afw.h"
#include "lsst/pex/policy.h"

#include "lsst/afw/detection/Psf.h"
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
    template <typename ImageT>
    class PsfCandidate : public lsst::afw::math::SpatialCellImageCandidate<ImageT> {
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;
    public:
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getXCenter;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getYCenter;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getWidth;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getHeight;
    
        typedef boost::shared_ptr<PsfCandidate<ImageT> > Ptr;
        typedef boost::shared_ptr<const PsfCandidate<ImageT> > ConstPtr;
        typedef std::vector<Ptr > PtrList;
        
        /**
         * Construct a PsfCandidate from a specified source and image.
         *
         * The x/yCenter is set to source.getX/YAstrom()
         */
        PsfCandidate(
            CONST_PTR(afw::table::SourceRecord) const & source, ///< The detected Source
            typename ImageT::ConstPtr parentImage ///< The image wherein lie the Sources
        ) :
            lsst::afw::math::SpatialCellImageCandidate<ImageT>(source->getX(), source->getY()),
            _parentImage(parentImage),
            _offsetImage(),
            _source(source),
            _haveImage(false), _amplitude(0.0), _var(1.0) { }
        
        /**
         * Construct a PsfCandidate from a specified source, image and xyCenter.
         */
        PsfCandidate(
            CONST_PTR(afw::table::SourceRecord) const& source, ///< The detected Source
            typename ImageT::ConstPtr parentImage, ///< The image wherein lie the Sources
            double xCenter,    ///< the desired x center
            double yCenter     ///< the desired y center
        ) :
            lsst::afw::math::SpatialCellImageCandidate<ImageT>(xCenter, yCenter),
            _parentImage(parentImage),
            _offsetImage(),
            _source(source),
            _haveImage(false), _amplitude(0.0), _var(1.0) { }
        
        /// Destructor
        virtual ~PsfCandidate() {};
        
        /**
         * Return Cell rating
         * 
         * @note Required method for use by SpatialCell
         */
        double getCandidateRating() const { return _source->getPsfFlux(); }
        
        /// Return the original Source
        CONST_PTR(afw::table::SourceRecord) getSource() const { return _source; }
        
        /// Return the best-fit amplitude
        double getAmplitude() const { return _amplitude; }
    
        /// Set the best-fit amplitude
        void setAmplitude(double amplitude) { _amplitude = amplitude; }

        /// Return the variance in use when fitting this object
        double getVar() const { return _var; }
    
        /// Set the variance to use when fitting this object
        void setVar(double var) { _var = var; }
    
        CONST_PTR(ImageT) getImage() const;
        PTR(ImageT) getOffsetImage(std::string const algorithm, unsigned int buffer) const;
    
        /// Return the number of pixels being ignored around the candidate image's edge
        static int getBorderWidth() { return _border; }
    
        /// Set the number of pixels to ignore around the candidate image's edge
        static void setBorderWidth(int border) { _border = border; }
    private:
        PTR(ImageT) extractImage(unsigned int width, unsigned int height) const;
        CONST_PTR(ImageT) _parentImage; // the %image that the Sources are found in
        PTR(ImageT) mutable _offsetImage; // %image offset to put center on a pixel
        CONST_PTR(afw::table::SourceRecord) _source; // the Source itself
        bool mutable _haveImage;                    // do we have an Image to return?
        double _amplitude;                          // best-fit amplitude of current PSF model
        double _var;                                // variance to use when fitting this candidate
        static int _border;                         // width of border of ignored pixels around _image
        lsst::afw::geom::Point2D _xyCenter;
    };
        
    template <typename ImageT>
    PTR(PsfCandidate<ImageT>)
    makePsfCandidate(
        CONST_PTR(afw::table::SourceRecord) const & source, ///< The detected Source
        PTR(ImageT)const & image                            ///< The image wherein lies the object
    ) {
        return boost::make_shared< PsfCandidate<ImageT> >(source, image);
    }
   
}}}

#endif
