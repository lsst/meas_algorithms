// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Class used by SpatialCell for spatial Kernel fitting
 *
 * @author Andrew Becker, University of Washington
 *
 * @ingroup afw
 */

#ifndef LSST_IP_DIFFIM_SPATIALMODELPSF_H
#define LSST_IP_DIFFIM_SPATIALMODELPSF_H

#include "boost/shared_ptr.hpp"

#include "lsst/afw.h"
#include "lsst/pex/policy/Policy.h"
//#include "lsst/sdqa/SdqaRating.h"

#include "lsst/afw/math/SpatialCell.h"

namespace lsst {
namespace meas {
namespace algorithms {
    /** 
     * 
     * @brief Class stored in SpatialCells for spatial Psf fitting
     * 
     * PsfCandidate is a detection that may turn out to be a PSF. We'll
     * assign them to sets of SpatialCells; these sets will then be used to fit
     * a spatial model to the PSF.
     */    
    template <typename ImageT>
    class PsfCandidate : public lsst::afw::math::SpatialCellImageCandidate<ImageT> {
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getXCenter;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getYCenter;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getWidth;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getHeight;
        using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;
    public: 
        typedef boost::shared_ptr<PsfCandidate> Ptr;

        /**
         * Constructor
         */
        PsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                     typename ImageT::ConstPtr parentImage ///< The image wherein lie the Sources
                    ) :
            lsst::afw::math::SpatialCellImageCandidate<ImageT>(source.getXAstrom(), source.getYAstrom()),
            _parentImage(parentImage),
            _source(source),
            _haveImage(false) {
        }

        /**
         * Destructor
         */
        ~PsfCandidate() {};

        /**
         * Return Cell rating
         * 
         * @note Required method for use by SpatialCell
         */
        double getCandidateRating() const {
            return _source.getPsfMag();
        }
        
        /**
         * Return the original Source
         */
        lsst::afw::detection::Source const& getSource() const { return _source; }

        typename ImageT::ConstPtr getImage() const;
    private:
        typename ImageT::ConstPtr _parentImage; // the %image that the Sources are found in
        lsst::afw::detection::Source const _source; // the Source itself
        bool _haveImage;                            // do we have an Image to return?
    };

    /**
     * Return a PsfCandidate of the right sort
     *
     * Cf. std::make_pair
     */
    template <typename ImageT>
    typename PsfCandidate<ImageT>::Ptr
    makePsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                     typename ImageT::ConstPtr image ///< The image wherein lies the object
                    ) {
        
        return typename PsfCandidate<ImageT>::Ptr(new PsfCandidate<ImageT>(source, image));
    }

    std::pair<lsst::afw::math::LinearCombinationKernel::PtrT, std::vector<double> >
    createKernelFromPsfCandidates(lsst::afw::math::SpatialCellSet const& psfCells,
                                  int const nEigenComponents,
                                  int const spatialOrder,
                                  int const ksize,
                                  int const nStarPerCell=-1                                  
                                 );

}}}

#endif
