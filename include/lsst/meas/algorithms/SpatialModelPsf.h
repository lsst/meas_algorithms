// -*- LSST-C++ -*-
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
//#include "lsst/sdqa/SdqaRating.h"

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
    using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getXCenter;
    using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getYCenter;
    using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getWidth;
    using lsst::afw::math::SpatialCellImageCandidate<ImageT>::getHeight;
    using lsst::afw::math::SpatialCellImageCandidate<ImageT>::_image;
public: 
    typedef boost::shared_ptr<PsfCandidate> Ptr;
    
    /// Constructor
    PsfCandidate(lsst::afw::detection::Source const& source, ///< The detected Source
                 typename ImageT::ConstPtr parentImage ///< The image wherein lie the Sources
                ) :
        lsst::afw::math::SpatialCellImageCandidate<ImageT>(source.getXAstrom(), source.getYAstrom()),
        _parentImage(parentImage),
        _source(source),
        _haveImage(false), _amplitude(0.0) { }
    
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

    typename ImageT::ConstPtr getImage() const;
private:
    typename ImageT::ConstPtr _parentImage; // the %image that the Sources are found in
    lsst::afw::detection::Source const _source; // the Source itself
    bool mutable _haveImage;                    // do we have an Image to return?
    double _amplitude;                          // best-fit amplitude of current PSF model
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
std::pair<bool, double>
fitSpatialKernelFromPsfCandidates(lsst::afw::math::Kernel *kernel,
                                  lsst::afw::math::SpatialCellSet const& psfCells,
                                  int const nStarPerCell = -1,
                                  double const tolerance = 1e-5);
    
template<typename ImageT>
double subtractPsf(lsst::meas::algorithms::PSF const& psf, ImageT *data, double x, double y);
}}}

#endif
