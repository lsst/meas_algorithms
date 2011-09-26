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
 
#if !defined(LSST_MEAS_ALGORITHMS_FLAGS_H)
#define LSST_MEAS_ALGORITHMS_FLAGS_H

namespace lsst { namespace meas { namespace algorithms {


/************************************************************************************************************/
/**
 * A class to provide a set of flags describing our processing
 *
 * This would be in MeasureSources, but that's inconvenient as it's templated
 */
struct Flags {
    enum {
        EDGE                            = 0x000001, ///< source is in region labelled EDGE
        SHAPE_SHIFT                     = 0x000002, ///< centroid shifted while estimating adaptive moments
        SHAPE_MAXITER                   = 0x000004, ///< too many iterations for adaptive moments
        SHAPE_UNWEIGHTED                = 0x000008, ///< "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_PSF            = 0x000010, ///< the PSF's "adaptive" moments are unweighted
        SHAPE_UNWEIGHTED_BAD            = 0x000020, ///< even the unweighted moments were bad
        PEAKCENTER                      = 0x000040, ///< given centre is position of peak pixel
        BINNED1                         = 0x000080, ///< source was found in 1x1 binned image
        INTERP                          = 0x000100, ///< source's footprint includes interpolated pixels
        INTERP_CENTER                   = 0x000200, ///< source's centre is close to interpolated pixels
        SATUR                           = 0x000400, ///< source's footprint includes saturated pixels
        SATUR_CENTER                    = 0x000800, ///< source's centre is close to saturated pixels
        DETECT_NEGATIVE                 = 0x001000, ///< source was detected as being significantly negative
        STAR                            = 0x002000, ///< source is thought to be point-like
        PSFSTAR                         = 0x004000, ///< source was used in PSF determination

        PHOTOM_NO_PSF                   = 0x008000, ///< NO Psf provided to photometry algorithm
        PHOTOM_NO_PEAK                  = 0x010000, ///< NO Peak provided to photometry algorithm
        PHOTOM_NO_SOURCE                = 0x020000, ///< NO source provided to photometry algorithm
        PHOTOM_NO_FOOTPRINT             = 0x040000, ///< NO FOOTPRINT provided to photometry algorithm

        SHAPELET_PHOTOM_NO_BASIS        = 0x080000, ///< ShapeletModelPhotometry configure without a basis
        SHAPELET_PHOTOM_BAD_MOMENTS     = 0x100000, ///< input moments are too large or not finite
        SHAPELET_PHOTOM_INVERSION_FAIL  = 0x200000, ///< ShapeletModelPhotometry failed
        SHAPELET_PHOTOM_INVERSION_UNSAFE= 0x400000, ///< ShapeletModelPhotometry should not be trusted
        SHAPELET_PHOTOM_GALAXY_FAIL     = 0x800000, ///< ShapeletModelPhotometry only fit a point source model
        ///ShapeletModelPhtoometry should be ignored in essentially all analyses
        SHAPELET_PHOTOM_BAD = PHOTOM_NO_PSF | PHOTOM_NO_SOURCE | PHOTOM_NO_FOOTPRINT | 
                SHAPELET_PHOTOM_NO_BASIS | SHAPELET_PHOTOM_BAD_MOMENTS | 
                SHAPELET_PHOTOM_INVERSION_FAIL | SHAPELET_PHOTOM_INVERSION_UNSAFE,
        /// Should this this object be ignored in essentially all analyses?
        BAD                       = EDGE|INTERP_CENTER|SATUR_CENTER

    };
};

}}} // namespace

#endif
