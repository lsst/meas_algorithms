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
 * Return the %image at the position of the Source, centered in a pixel
 */
template <typename ImageT>
typename ImageT::ConstPtr lsst::meas::algorithms::PsfCandidate<ImageT>::getImage() const {
    int const width = getWidth() == 0 ? 15 : getWidth();
    int const height = getHeight() == 0 ? 15 : getHeight();

    if (_haveImage && (width != _image->getWidth() || height != _image->getHeight())) {
        _haveImage = false;
    }

    if (!_haveImage) {
        std::pair<int, double> xCen = afwImage::positionToIndex(getXCenter(), true); // true => return the std::pair
        std::pair<int, double> yCen = afwImage::positionToIndex(getYCenter(), true);

        afwImage::PointI center(xCen.first, yCen.first); // integral part
        afwImage::BBox bbox(center - afwImage::PointI(width/2, height/2), width, height);
        bbox.shift(-_parentImage->getX0(), -_parentImage->getY0());
        
        try {
            typename ImageT::Ptr patch(new ImageT(*_parentImage, bbox, false)); // a shallow copy

            std::string algorithmName = "lanczos5";
            _image = afwMath::offsetImage(*patch, -xCen.second, -yCen.second, algorithmName);
            
            _haveImage = true;
        } catch(lsst::pex::exceptions::LengthErrorException &e) {
            LSST_EXCEPT_ADD(e, "Setting image for PSF candidate");
            throw e;
        }
    }
    
    return _image;
}

/************************************************************************************************************/

namespace {
    // A class to pass around to all our PsfCandidates which builds the PcaImageSet
    template<typename PixelT>
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

            try {
                _imagePca->addImage(imCandidate->getImage()->getImage(), imCandidate->getSource().getPsfFlux());
            } catch(lsst::pex::exceptions::LengthErrorException &e) {
                return;
            }
        }
    private:
        afwImage::ImagePca<ImageT> *_imagePca; // the ImagePca we're building
    };
}

/************************************************************************************************************/
/**
 * Return a Kernel::Ptr and a list of eigenvalues resulting from analysing the provided SpatialCellSet
 *
 * The Kernel is a LinearCombinationKernel of the first nEigenComponents eigenImages
 *
 * N.b. This is templated over the Pixel type of the science image
 */
template<typename PixelT>
std::pair<afwMath::LinearCombinationKernel::Ptr, std::vector<double> > createKernelFromPsfCandidates(
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        int const nEigenComponents,     ///< number of eigen components to keep; <= 0 => infty
        int const spatialOrder,         ///< Order of spatial variation (cf. lsst::afw::math::PolynomialFunction2
        int const ksize,                ///< Size of generated Kernel images
        int const nStarPerCell          ///< max no. of stars per cell; <= 0 => infty
                                                                                    ) {
    typedef typename afwImage::Image<PixelT> ImageT;
    typedef typename afwImage::MaskedImage<PixelT> MaskedImageT;
    //
    // Set the sizes for PsfCandidates made from either Images or MaskedImages
    //
    lsst::meas::algorithms::PsfCandidate<ImageT>::setWidth(ksize);
    lsst::meas::algorithms::PsfCandidate<ImageT>::setHeight(ksize);
    lsst::meas::algorithms::PsfCandidate<MaskedImageT>::setWidth(ksize);
    lsst::meas::algorithms::PsfCandidate<MaskedImageT>::setHeight(ksize);

    afwImage::ImagePca<ImageT> imagePca; // Here's the set of images we'll analyze

    SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
    psfCells.visitCandidates(&importStarVisitor, nStarPerCell);

    //
    // Do a PCA decomposition of those PSF candidates
    //
    imagePca.analyze();
    
    std::vector<typename ImageT::Ptr> eigenImages = imagePca.getEigenImages();
    std::vector<double> eigenValues = imagePca.getEigenValues();
    int const nEigen = static_cast<int>(eigenValues.size());
    
    int const ncomp = (nEigenComponents <= 0 || nEigen < nEigenComponents) ? nEigen : nEigenComponents;
    //
    // Now build our LinearCombinationKernel; build the lists of basis functions
    // and spatial variation, then assemble the Kernel
    //
    afwMath::KernelList  kernelList;
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    
    for (int i = 0; i != ncomp; ++i) {
        kernelList.push_back(afwMath::Kernel::Ptr(
                new afwMath::FixedKernel(afwImage::Image<afwMath::Kernel::Pixel>(*eigenImages[i], true))));

        afwMath::Kernel::SpatialFunctionPtr spatialFunction(new afwMath::PolynomialFunction2<double>(spatialOrder));
        spatialFunction->setParameter(0, 1.0); // the constant term; all others are 0
        spatialFunctionList.push_back(spatialFunction);
    }

    afwMath::LinearCombinationKernel::Ptr psf(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));

    return std::make_pair(psf, eigenValues);
}


/************************************************************************************************************/
namespace {
/**
 * Fit the model mImage to the data;  the model is assumed to have been shifted to have the same centroid
 *
 * Return (chi^2, amplitude) where amplitude*model is the best fit to the data
 */
template<typename ModelImageT, typename DataImageT>
std::pair<double, double>
fitKernel(ModelImageT const& mImage,    // The model image at this point
          DataImageT const& data        // the data to fit
         ) {
    assert(data.getDimensions() == mImage.getDimensions());

    double sumMM = 0.0, sumMD = 0.0, sumDD = 0.0; // sums of model*model/variance etc.
    for (int y = 0; y != data.getHeight(); ++y) {
        typename ModelImageT::x_iterator mptr = mImage.row_begin(y);
        for (typename DataImageT::x_iterator ptr = data.row_begin(y), end = data.row_end(y);
             ptr != end; ++ptr, ++mptr) {
            double const M = (*mptr)[0];       // value of model
            double const D = ptr.image();      // value of data
            double const var = ptr.variance(); // data's variance
            if (var != 0.0) {                  // assume variance == 0 => infinity XXX
                double const iVar = 1.0/var;
                sumMM += M*M*iVar;
                sumMD += M*D*iVar;
                sumDD += D*D*iVar;
            }
        }
    }
        
    if (sumMM == 0.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "sum(data*data)/var == 0");
    }

    double const amp = sumMD/sumMM;     // estimate of amplitude of model at this point            
    double chi2 = sumDD - 2*amp*sumMD + amp*amp*sumMM; // chi^2 from this object

#if 1
    static volatile bool show = false;   // Display the centre of the image; set from gdb
        
    if (show) {
        int y = data.getHeight()/2;
        int x = data.getWidth()/2;
        int hsize = 2;
        printf("\ndata  ");
        for (int ii = -hsize; ii <= hsize; ++ii) {
            for (int jj = -hsize; jj <= hsize; ++jj) {
                printf("%7.1f ", data.at(x + ii, y - jj).image());
            }
            printf("  model  ");
            for (int jj = -hsize; jj <= hsize; ++jj) {
                printf("%7.1f ", amp*(*(mImage.at(x + ii, y - jj)))[0]);
            }
            printf("\n      ");
        }
    }
#endif
        
    return std::make_pair(chi2, amp);
}
}

/************************************************************************************************************/
/*
 * Fit for the spatial variation of the PSF parameters over the field
 */
namespace {
#if !defined(DOXYGEN)
#   include "Minuit/FCNBase.h"
#   include "Minuit/FunctionMinimum.h"
#   include "Minuit/MnMigrad.h"
#   include "Minuit/MnMinos.h"
#   include "Minuit/MnPrint.h"
#endif
    
    /// A class to pass around to all our PsfCandidates to evaluate the PSF fit's X^2 
    template<typename PixelT>
    class evalChi2Visitor : public afwMath::CandidateVisitor {
        typedef afwImage::Image<PixelT> Image;
        typedef afwImage::MaskedImage<PixelT> MaskedImage;

        typedef afwImage::Image<afwMath::Kernel::Pixel> KImage;
    public:
        evalChi2Visitor(afwMath::Kernel const& kernel
                        ) :
            afwMath::CandidateVisitor(),
            _chi2(0.0), _kernel(kernel),
            _kImage(KImage::Ptr(new KImage(kernel.getDimensions()))) {
        }
        
        void reset() {
            _chi2 = 0.0;
        }

        // Called by SpatialCellSet::visitCandidates for each Candidate
        void processCandidate(afwMath::SpatialCellCandidate *candidate) {
            PsfCandidate<MaskedImage> *imCandidate = dynamic_cast<PsfCandidate<MaskedImage> *>(candidate);
            if (imCandidate == NULL) {
                throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                                  "Failed to cast SpatialCellCandidate to PsfCandidate");
            }

            _kernel.computeImage(*_kImage, false,
                                 imCandidate->getSource().getXAstrom(), imCandidate->getSource().getYAstrom());
            typename MaskedImage::ConstPtr data;
            try {
                data = imCandidate->getImage();
            } catch(lsst::pex::exceptions::LengthErrorException &e) {
                return;
            }

            try {
                std::pair<double, double> result = fitKernel(*_kImage, *data);

                double dchi2 = result.first; // chi^2 from this object
#if 0
                double const amp = result.second; // estimate of amplitude of model at this point
#endif
                
                imCandidate->setChi2(dchi2);
                
                _chi2 += dchi2;
            } catch(lsst::pex::exceptions::RangeErrorException &e) {
                LSST_EXCEPT_ADD(e, (boost::format("Object at (%.2f, %.2f)") %
                                    imCandidate->getSource().getXAstrom() %
                                    imCandidate->getSource().getYAstrom()).str());
                throw e;
            }
        }

        // Return the computed chi^2
        double getValue() const { return _chi2; }

    private:
        double mutable _chi2;            // the desired chi^2
        afwMath::Kernel const& _kernel;  // the kernel
        typename KImage::Ptr mutable _kImage; // The Kernel at this point; a scratch copy
    };

    /************************************************************************************************************/
    /**
     * Fit a Kernel's spatial variability from a set of stars
     *
     * N.b. This is templated over the Pixel type of the science image
     */
    // Set the Kernel's spatial parameters from a vector of length(nComponents*nSpatialParams)
    void setSpatialParameters(afwMath::Kernel *kernel,
                              std::vector<double> const& coeffs,
                              int const nComponents,
                              int const nSpatialParams) {
        std::vector<std::vector<double> > kCoeffs; // coefficients rearranged for Kernel
        kCoeffs.reserve(nComponents);
        for (int i = 0; i != nComponents; ++i) {
            kCoeffs.push_back(std::vector<double>(nSpatialParams));
            std::copy(coeffs.begin() + i*nSpatialParams, coeffs.begin() + (i + 1)*nSpatialParams, kCoeffs[i].begin());
        }

        kernel->setSpatialParameters(kCoeffs);
    }
    //
    // The object that minuit minimises
    //
    template<typename PixelT>
    class MinimizeChi2 : public FCNBase {
    public:
        explicit MinimizeChi2(evalChi2Visitor<PixelT> const& chi2Visitor,
                              afwMath::Kernel *kernel,
                              afwMath::SpatialCellSet const& psfCells,
                              int nStarPerCell,
                              int nComponents,
                              int nSpatialParams
                             ) : _errorDef(1.0),
                                 _chi2Visitor(chi2Visitor),
                                 _kernel(kernel),
                                 _psfCells(psfCells),
                                 _nStarPerCell(nStarPerCell),
                                 _nComponents(nComponents),
                                 _nSpatialParams(nSpatialParams)
            {
        }
        /**
         * Error definition of the function. MINUIT defines Parameter errors as the
         * change in Parameter Value required to change the function Value by up. Normally,
         * for chisquared fits it is 1, and for negative log likelihood, its Value is 0.5.
         * If the user wants instead the 2-sigma errors for chisquared fits, it becomes 4,
         */
        double up() const {return _errorDef;}
        
        // Evaluate our cost function (in this case chi^2)
        double operator()(const std::vector<double>& coeffs) const {
            setSpatialParameters(_kernel, coeffs, _nComponents, _nSpatialParams);

            _psfCells.visitCandidates(&_chi2Visitor, _nStarPerCell);

            return _chi2Visitor.getValue();
        }
        
        void setErrorDef(double def) { _errorDef = def; }
    private:
        double _errorDef;               // how much cost function has changed at the +- 1 error points

        evalChi2Visitor<PixelT> const& _chi2Visitor;
        afwMath::Kernel *_kernel;
        afwMath::SpatialCellSet const& _psfCells;
        int _nStarPerCell;
        int _nComponents;
        int _nSpatialParams;
    };
}

/************************************************************************************************************/
    
template<typename PixelT>
std::pair<bool, double>
fitSpatialKernelFromPsfCandidates(
        afwMath::Kernel *kernel,                 ///< the Kernel to fit
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        int const nStarPerCell,                  ///< max no. of stars per cell; <= 0 => infty
        double const tolerance                   ///< Tolerance; how close chi^2 should be to true minimum
                                 ) {
    typedef typename afwImage::Image<PixelT> Image;

    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();
    //
    // visitor that evaluates the chi^2 of the current fit
    //
    evalChi2Visitor<PixelT> getChi2(*kernel);
    //
    // We have to unpack the Kernel coefficients into a linear array, coeffs
    //
    std::vector<double> coeffs;         // The coefficients we want to fit
    coeffs.assign(nComponents*nSpatialParams, 0.0);

    std::vector<double> stepSize;       // step sizes
    stepSize.assign(nComponents*nSpatialParams, 100);
    //
    // Translate that into minuit's language
    //
    MnUserParameters fitPar;
    std::vector<std::string> paramNames;
    paramNames.reserve(nComponents*nSpatialParams);
    
    for (int i = 0, c = 0; c != nComponents; ++c) {
        coeffs[i] = 1;                  // the constant part of each spatial order
        for (int s = 0; s != nSpatialParams; ++s, ++i) {
            paramNames.push_back((boost::format("C%d:%d") % c % s).str());
            fitPar.add(paramNames[i].c_str(), coeffs[i], stepSize[i]);
        }
    }
    fitPar.fix("C0:0");
    //
    // Create the minuit object that knows how to minimise our functor
    //
    MinimizeChi2<PixelT> minimizerFunc(getChi2, kernel, psfCells, nStarPerCell, nComponents, nSpatialParams);

    const double errorDef = 1.0;       // use +- 1sigma errors
    minimizerFunc.setErrorDef(errorDef);
    //
    // tell minuit about it
    //    
    MnMigrad migrad(minimizerFunc, fitPar);
    //
    // And let it loose
    //
    int maxFnCalls = 0;                 // i.e. unlimited
    FunctionMinimum min = migrad(maxFnCalls, tolerance/(1e-4*errorDef)); // minuit uses 0.1*1e-3*tolerance*errorDef

    float minChi2 = min.fval();
    bool const isValid = min.isValid() && std::isfinite(minChi2);
    
    if (isValid) {
        for (int i = 0; i != nComponents*nSpatialParams; ++i) {
            coeffs[i] = min.userState().value(i);
        }

        setSpatialParameters(kernel, coeffs, nComponents, nSpatialParams);
    }

#if 0                                   // Estimate errors;  we don't really need this
    MnMinos minos(minimizerFunc, min);
    for (int i = 0, c = 0; c != nComponents; ++c) {
        for (int s = 0; s != nSpatialParams; ++s, ++i) {
            char const *name = paramNames[i].c_str();
            printf("%s %g", name, min.userState().value(name));
            if (isValid && !fitPar.parameter(fitPar.index(name)).isFixed()) {
                printf(" (%g+%g)\n", minos(i).first, minos(i).second);
            }
            printf("\n");
        }
    }
#endif
    //
    // One time more through the Candidates setting their chi^2 values. We'll
    // do all the candidates this time, not just the first nStarPerCell
    //
    psfCells.visitCandidates(&getChi2, 0);
    
    return std::make_pair(isValid, minChi2);
}

/************************************************************************************************************/
/**
 * Subtract a PSF from an image at a given position
 */
template<typename MaskedImageT>
double subtractPsf(PSF const& psf,      ///< the PSF to subtract
                   MaskedImageT *data,  ///< Image to subtract PSF from
                   double x,            ///< column position
                   double y             ///< row position
                  ) {
    typedef afwImage::Image<afwMath::Kernel::Pixel> KernelImage;

    int const width = psf.getWidth();
    int const height = psf.getHeight();
    //
    // Shift the PSF to the proper centre
    //
    KernelImage::Ptr kImageD(new KernelImage(width, height));
    psf.getKernel()->computeImage(*kImageD, false, x, y);

    std::string algorithmName = "lanczos5";
    KernelImage::Ptr kImage = afwMath::offsetImage(*kImageD, x, y, algorithmName);
    //
    // Now find the proper sub-Image
    //
    afwImage::BBox bbox(afwImage::PointI(0, 0), width, height);
    bbox.shift(kImage->getX0() - width/2, kImage->getY0() - height/2);
    
    typename MaskedImageT::Ptr subData(new MaskedImageT(*data, bbox, false)); // a shallow copy
    //
    // Now we've got both; find the PSF's amplitude
    //
    try {
        std::pair<double, double> result = fitKernel(*kImage, *subData);
        double const chi2 = result.first; // chi^2 for fit
        double const amp = result.second; // estimate of amplitude of model at this point
        //
        // Convert kImage to the proper type so that I can subtract it.
        //
        typename MaskedImageT::Image::Ptr kImageF(new typename MaskedImageT::Image(*kImage, true)); // of data's type

        *kImageF *= amp;
        *subData->getImage() -= *kImageF;
        
        return chi2;
    } catch(lsst::pex::exceptions::RangeErrorException &e) {
        LSST_EXCEPT_ADD(e, (boost::format("Object at (%.2f, %.2f)") % x % y).str());
        throw e;
    }
}
    
/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
    typedef float Pixel;
    template class PsfCandidate<afwImage::MaskedImage<Pixel> >;

    template
    std::pair<lsst::afw::math::LinearCombinationKernel::Ptr, std::vector<double> >
    createKernelFromPsfCandidates<Pixel>(lsst::afw::math::SpatialCellSet const&,
                                          int const, int const, int const, int const);
    template
    std::pair<bool, double>
    fitSpatialKernelFromPsfCandidates<Pixel>(afwMath::Kernel *, afwMath::SpatialCellSet const&,
                                              int const, double const);

    template
    double subtractPsf(PSF const&, afwImage::MaskedImage<float> *, double, double);

/// \endcond
}}}
