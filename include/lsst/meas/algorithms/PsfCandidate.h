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
#include <memory>
#include <vector>

#include "lsst/geom/Point.h"
#include "lsst/afw/image/Exposure.h"
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
class PsfCandidate : public afw::math::SpatialCellImageCandidate {
public:
    typedef std::shared_ptr<PsfCandidate<PixelT>> Ptr;
    typedef std::shared_ptr<const PsfCandidate<PixelT>> ConstPtr;
    typedef std::vector<Ptr> PtrList;

    typedef afw::image::MaskedImage<PixelT> MaskedImageT;

    /**
     * Construct a PsfCandidate from a specified source and image.
     *
     * The x/yCenter is set to source.getX/YAstrom()
     */
    PsfCandidate(PTR(afw::table::SourceRecord) const & source,  ///< The detected Source
                 CONST_PTR(afw::image::Exposure<PixelT>)
                         parentExposure  ///< The image wherein lie the Sources
                 )
            : afw::math::SpatialCellImageCandidate(source->getX(), source->getY()),
              _parentExposure(parentExposure),
              _offsetImage(),
              _source(source),
              _image(nullptr),
              _amplitude(0.0),
              _var(1.0) {}

    /**
     * Construct a PsfCandidate from a specified source, image and xyCenter.
     */
    PsfCandidate(PTR(afw::table::SourceRecord) const & source,  ///< The detected Source
                 CONST_PTR(afw::image::Exposure<PixelT>)
                         parentExposure,  ///< The image wherein lie the Sources
                 double xCenter,          ///< the desired x center
                 double yCenter           ///< the desired y center
                 )
            : afw::math::SpatialCellImageCandidate(xCenter, yCenter),
              _parentExposure(parentExposure),
              _offsetImage(),
              _source(source),
              _image(nullptr),
              _amplitude(0.0),
              _var(1.0) {}

    /// Destructor
    virtual ~PsfCandidate(){};

    /**
     * Return Cell rating
     *
     * @note Required method for use by SpatialCell
     */
    double getCandidateRating() const { return _source->getPsfInstFlux(); }

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

    CONST_PTR(afw::image::MaskedImage<PixelT>) getMaskedImage() const;
    CONST_PTR(afw::image::MaskedImage<PixelT>) getMaskedImage(int width, int height) const;
    PTR(afw::image::MaskedImage<PixelT>)
    getOffsetImage(std::string const algorithm, unsigned int buffer) const;

    /// Return the number of pixels being ignored around the candidate image's edge
    static int getBorderWidth();

    /// Set the number of pixels to ignore around the candidate image's edge
    static void setBorderWidth(int border);

    /// Set threshold for rejecting pixels unconnected with the central footprint
    ///
    /// A non-positive threshold means that no threshold will be applied.
    static void setPixelThreshold(float threshold);

    /// Get threshold for rejecting pixels unconnected with the central footprint
    static float getPixelThreshold();

    /// Set whether blends are masked
    static void setMaskBlends(bool doMaskBlends);

    /// Get whether blends are masked
    static bool getMaskBlends();

private:
    CONST_PTR(afw::image::Exposure<PixelT>)
    _parentExposure;  // the %image that the Sources are found in

    PTR(afw::image::MaskedImage<PixelT>)
    offsetImage(PTR(afw::image::MaskedImage<PixelT>) img, std::string const algorithm, unsigned int buffer);

    PTR(afw::image::MaskedImage<PixelT>)
    extractImage(unsigned int width, unsigned int height) const;

    PTR(afw::image::MaskedImage<PixelT>) mutable _offsetImage;  // %image offset to put center on a pixel
    PTR(afw::table::SourceRecord) _source;                      // the Source itself

    mutable std::shared_ptr<afw::image::MaskedImage<PixelT>> _image;  // cutout image to return (cached)
    double _amplitude;   // best-fit amplitude of current PSF model
    double _var;         // variance to use when fitting this candidate
    static int _border;  // width of border of ignored pixels around _image
    geom::Point2D _xyCenter;
    static int _defaultWidth;
    static float _pixelThreshold;  ///< Threshold for masking pixels unconnected with central footprint
    static bool _doMaskBlends;     ///< Mask blends when extracting?
};

/**
 * Return a PsfCandidate of the right sort
 *
 * Cf. std::make_pair
 */
template <typename PixelT>
std::shared_ptr<PsfCandidate<PixelT>> makePsfCandidate(PTR(afw::table::SourceRecord)
                                                               const & source,  ///< The detected Source
                                                       PTR(afw::image::Exposure<PixelT>)
                                                               image  ///< The image wherein lies the object
                                                       ) {
    return std::make_shared<PsfCandidate<PixelT>>(source, image);
}

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif
