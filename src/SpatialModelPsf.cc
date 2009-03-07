/*!
 * \brief Implementation of code to determine spatial model of PSF
 *
 * \file
 *
 * \ingroup algorithms
 */
#include <typeinfo>
#include <cmath>
#include "lsst/afw/image/ImagePca.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"

namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

/************************************************************************************************************/
/*
 * PsfCandidate's members
 */
/**
 * Return the %image at the position of the Source
 */
template <typename ImageT>
typename ImageT::ConstPtr lsst::meas::algorithms::PsfCandidate<ImageT>::getImage() const {
    if (!_haveImage) {
        int const width = getWidth() == 0 ? 15 : getWidth();
        int const height = getHeight() == 0 ? 15 : getHeight();
    
        std::pair<int, double> xCen = afwImage::positionToIndex(getXCenter(), true); // true => return the std::pair
        std::pair<int, double> yCen = afwImage::positionToIndex(getYCenter(), true);

        afwImage::PointI center(xCen.first, yCen.first); // integral part
        afwImage::BBox bbox(center - afwImage::PointI(width/2, height/2), width, height);
        bbox.shift(-_parentImage->getX0(), -_parentImage->getY0());
        
        typename ImageT::Ptr patch(new ImageT(*_parentImage, bbox, false)); // a shallow copy

        std::string algorithmName = "lanczos5";
        _image = afwMath::offsetImage(algorithmName, *patch, -xCen.second, -yCen.second);
    }
    
    return _image;
}

/************************************************************************************************************/

namespace {
    typedef float PixelT;

    // A class to pass around to all our PsfCandidates which builds the PcaImageSet
    class SetPcaImageVisitor : public afwMath::CandidateVisitor {
        typedef afwImage::Image<PixelT> ImageT;
        typedef afwImage::MaskedImage<PixelT> MaskedImageT;
    public:
        SetPcaImageVisitor(afwImage::ImagePca<ImageT> *imagePca // Set of Images to initialise
                          ) :
            afwMath::CandidateVisitor(),
            _imagePca(imagePca) {}
        
        // Called by SpatialCellSet::visitCandidates for each Candidate
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            PsfCandidate<MaskedImageT> *imCandidate = dynamic_cast<PsfCandidate<MaskedImageT> *>(candidate);
            if (imCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to PsfCandidate");
            }
            _imagePca->addImage(imCandidate->getImage()->getImage(), imCandidate->getSource().getPsfMag());
        }
    private:
        afwImage::ImagePca<ImageT> *_imagePca; // the ImagePca we're building
    };
}

/**
 * Return a Kernel::PtrT and a list of eigenvalues resulting from analysing the provided SpatialCellSet
 *
 * The Kernel is a LinearCombinationKernel of the first nEigenComponents eigenImages
 */
std::pair<afwMath::LinearCombinationKernel::PtrT, std::vector<double> > createKernelFromPsfCandidates(
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        int const nEigenComponents,     ///< number of eigen components to keep; <= 0 => infty
        int const spatialOrder,         ///< Order of spatial variation (cf. lsst::afw::math::PolynomialFunction2
        int const ksize,                ///< Size of generated Kernel images
        int const nStarPerCell          ///< max no. of stars per cell; <= 0 => infty
                                                                                    ) {
    typedef afwImage::Image<PixelT> ImageT;

    lsst::meas::algorithms::PsfCandidate<ImageT>::setWidth(ksize);
    lsst::meas::algorithms::PsfCandidate<ImageT>::setHeight(ksize);

    afwImage::ImagePca<ImageT> imagePca; // Here's the set of images we'll analyze

    SetPcaImageVisitor pcaVisitor(&imagePca);
    psfCells.visitCandidates(&pcaVisitor, nStarPerCell);

    //
    // Do a PCA decomposition of those PSF candidates
    //
    imagePca.analyze();
    
    std::vector<ImageT::Ptr> eigenImages = imagePca.getEigenImages();
    std::vector<double> eigenValues = imagePca.getEigenValues();
    int const nEigen = static_cast<int>(eigenValues.size());
    
    int const ncomp = (nEigenComponents <= 0 || nEigen < nEigenComponents) ? nEigen : nEigenComponents;
    //
    // Now build our LinearCombinationKernel; build the lists of basis functions
    // and spatial variation, then assemble the Kernel
    //
    afwMath::KernelList<afwMath::Kernel>  kernelList;
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    
    for (int i = 0; i != ncomp; ++i) {
        kernelList.push_back(afwMath::Kernel::PtrT(
                new afwMath::FixedKernel(afwImage::Image<afwMath::Kernel::PixelT>(*eigenImages[i], true))));

        afwMath::Kernel::SpatialFunctionPtr spatialFunction(new afwMath::PolynomialFunction2<double>(spatialOrder));
        spatialFunction->setParameter(0, 1.0); // the constant term; all others are 0
        spatialFunctionList.push_back(spatialFunction);
    }

    afwMath::LinearCombinationKernel::PtrT psf(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));

    return std::make_pair(psf, eigenValues);
}
    
//
// Explicit instantiations
//
// \cond
    template class PsfCandidate<afwImage::MaskedImage<float> >;
// \endcond
}}}
