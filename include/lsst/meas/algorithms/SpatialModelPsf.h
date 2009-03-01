// -*- lsst-c++ -*-
/**
 * @file
 *
 * @brief Class used by SpatialModelCell for spatial Kernel fitting
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
    class PsfCandidate : public lsst::afw::math::SpatialCellCandidate {
    public: 
        typedef boost::shared_ptr<PsfCandidate> Ptr;

        /**
         * Constructor
         */
        PsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                     typename ImageT::ConstPtr image) : ///< The image wherein lies the object
            SpatialCellCandidate(source.getXAstrom(), source.getYAstrom()),
            _image(image),
            _flux(source.getPsfMag()) {
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
            return _flux;
        }
    private:
        typename ImageT::ConstPtr _image; // The image containing the object
        double _flux;                     // the object's flux
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
}}}

#endif


