// -*- LSST-C++ -*-
#ifndef LSST_COADD_ALGORITHMS_SIZEMAGNITUDESTARFINDER_H
#define LSST_COADD_ALGORITHMS_SIZEMAGNITUDESTARFINDER_H

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
#include "boost/shared_ptr.hpp"

#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/Match.h"  // because we can't forward-declare a typedef
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/meas/algorithms/PsfCandidate.h"
#include "lsst/pex/policy/Policy.h"

namespace lsst {
namespace meas {
namespace algorithms {

    class SizeMagnitudeStarSelectorImpl;

    class SizeMagnitudeStarSelector 
    {
    public:
        typedef float PixelT;
        typedef lsst::afw::image::Exposure<PixelT> Exposure;
        typedef lsst::afw::geom::PointD PointD;
        typedef lsst::pex::policy::Policy Policy;
        typedef lsst::afw::table::SourceRecord SourceRecord;
        typedef lsst::afw::table::SourceTable SourceTable;
        typedef lsst::afw::table::SourceCatalog SourceCatalog;
        typedef lsst::meas::algorithms::PsfCandidate<PixelT> PsfCandidateT;

        typedef PsfCandidateT::PtrList PsfCandidateList;

        /*!
         * @brief Create a size-magnitude star selector from a policy
         *
         * @param[in] policy: see policy/SizeMagnitudeStarSelectorDictionary.paf
         */
        SizeMagnitudeStarSelector(const Policy& policy);

        /*!
         * @brief Destructor
         */
        ~SizeMagnitudeStarSelector() {};

        /*!
         * @brief Find a set of stars from an input list of Sources
         *
         * This function uses an algorithm based on looking for a stellar locus
         * in the size-magnitude diagram that has constant size.  
         *
         * It basically construct a histogram of counts with respect to size,
         * starting with just the bright objects and then pushing down in magnitude
         * until the stellar peak starts to bleed into the galaxies, at 
         * which point it stops.
         */
        PsfCandidateList selectStars(
            const Exposure & exposure,   ///< the exposure on which the sources are found
            const SourceCatalog & sourceList, ///< The input list of sources to consider
            CONST_PTR(afw::table::ReferenceMatchVector) matches
            = CONST_PTR(afw::table::ReferenceMatchVector)() ///< astrometric matches (ignored by this algorithm)
        ) const;

    private :
    
        /*!
         * @brief Calculates a robust size measurement for a source.
         *
         * This measures the 2nd order shapelet decomposition of the object,
         * allowing the sigma to vary until the b11 term goes to zero.
         * When this happens, the sigma is basically the best-fit Gaussian sigma.
         * Also, this value of sigma gives the best S/N properties for the
         * rest of the shapelet vector.
         *
         * This size measurement has proven to be a good one to use for 
         * the size-magnitude star finder algorithm, since the size is very
         * stable with stellar magnitude.
         *
         * @note This may not be necessary.  Lupton says that there is already
         * a similar measurement being done in the LSST stack.  Should check
         * if the existing values are ok for star finder algorithm.
         * (I'm not sure which entry in Source they would be...)
         * Therefore, this function may be a good candidate for a having
         * its action be specifiable by a Policy parameter.
         */
        double calculateSourceSize(
            const SourceRecord & source, 
            const Exposure & exposure
        ) const;

        /*!
         * @brief This class is not intended to be copied.
         */
        SizeMagnitudeStarSelector(const SizeMagnitudeStarSelector& rhs);
        void operator=(const SizeMagnitudeStarSelector& rhs);

        boost::shared_ptr<SizeMagnitudeStarSelectorImpl> pImpl;
    };

}}}

#endif
