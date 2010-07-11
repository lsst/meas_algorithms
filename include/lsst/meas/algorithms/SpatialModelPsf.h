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
 
#if !defined(LSST_IP_DIFFIM_SPATIALMODELPSF_H)
#define LSST_IP_DIFFIM_SPATIALMODELPSF_H
/**
 * @file
 *
 * @brief Class used by SpatialCell for spatial Kernel fitting
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup afw
 */
#include "boost/shared_ptr.hpp"

#include "lsst/afw.h"
#include "lsst/pex/policy.h"

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

    typedef boost::shared_ptr<PsfCandidate> Ptr;
    
    /// Constructor
    PsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                 typename ImageT::ConstPtr parentImage ///< The image wherein lie the Sources
                ) :
        lsst::afw::math::SpatialCellImageCandidate<ImageT>(source.getXAstrom(), source.getYAstrom()),
        _parentImage(parentImage),
        _source(source),
        _haveImage(false), _amplitude(0.0), _var(1.0) { }
    
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

    /// Set the best-fit amplitude
    void setAmplitude(double amplitude) { _amplitude = amplitude; }
    /// Return the best-fit amplitude
    double getAmplitude() const { return _amplitude; }

    /// Set the variance to use when fitting this object
    void setVar(double var) { _var = var; }
    /// Return the variance in use when fitting this object
    double getVar() const { return _var; }

    typename ImageT::ConstPtr getImage() const;

    /// Return the number of pixels being ignored around the candidate image's edge
    static int getBorderWidth() {
        return _border;
    }

    /// Set the number of pixels to ignore around the candidate image's edge
    static void setBorderWidth(int border) {
        _border = border;
    }
private:
    typename ImageT::ConstPtr _parentImage; // the %image that the Sources are found in
    lsst::afw::detection::Source const _source; // the Source itself
    bool mutable _haveImage;                    // do we have an Image to return?
    double _amplitude;                          // best-fit amplitude of current PSF model
    double _var;                                // variance to use when fitting this candidate
    static int _border;                         // width of border of ignored pixels around _image
};
    
/**
 * Return a PsfCandidate of the right sort
 *
 * Cf. std::make_pair
 *
 * \note It may be desirable to check that ImagePtrT really is a shared_ptr<image>.  The code is written this
 * way to allow the compiler to deduce the argument types
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
    
template<typename PixelT>
std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, std::vector<double> >
createKernelFromPsfCandidates(lsst::afw::math::SpatialCellSet const& psfCells,
                              int const nEigenComponents,
                              int const spatialOrder,
                              int const ksize,
                              int const nStarPerCell=-1,
                              bool const constantWeight=true                              
                             );

template<typename PixelT>
int countPsfCandidates(lsst::afw::math::SpatialCellSet const& psfCells, int const nStarPerCell=-1);
    
template<typename PixelT>
std::pair<bool, double>
fitSpatialKernelFromPsfCandidates(lsst::afw::math::Kernel *kernel,
                                  lsst::afw::math::SpatialCellSet const& psfCells,
                                  int const nStarPerCell = -1,
                                  double const tolerance = 1e-5);
template<typename PixelT>
std::pair<bool, double>
fitSpatialKernelFromPsfCandidates(lsst::afw::math::Kernel *kernel,
                                  lsst::afw::math::SpatialCellSet const& psfCells,
                                  bool const doNonLinearFit,
                                  int const nStarPerCell = -1,
                                  double const tolerance = 1e-5);
    
template<typename ImageT>
double subtractPsf(lsst::meas::algorithms::PSF const& psf, ImageT *data, double x, double y);

template<typename Image>
std::pair<lsst::afw::math::Kernel::Ptr, double>
fitKernelToImage(lsst::afw::math::LinearCombinationKernel const& kernel,
                 Image const& image, lsst::afw::geom::Point2D const& pos);
    
}}}

#endif
