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
 
/*!
 * @brief Implementation of code to determine spatial model of PSF
 *
 * @file
 *
 * @ingroup algorithms
 */
#if !defined(DOXYGEN)
#   include "Minuit2/FCNBase.h"
#   include "Minuit2/FunctionMinimum.h"
#   include "Minuit2/MnMigrad.h"
#   include "Minuit2/MnMinos.h"
#   include "Minuit2/MnPrint.h"
#endif

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/SVD"

#include "lsst/afw/image/ImagePca.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"
#include "lsst/meas/algorithms/PsfCandidate.h"

namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

// A class to pass around to all our PsfCandidates which builds the PcaImageSet
template<typename PixelT>
class SetPcaImageVisitor : public afwMath::CandidateVisitor {
    typedef afwImage::Image<PixelT> ImageT;
    typedef afwImage::MaskedImage<PixelT> MaskedImageT;
public:
    explicit SetPcaImageVisitor(
            afwImage::ImagePca<MaskedImageT> *imagePca, // Set of Images to initialise
            unsigned int const mask=0x0                    // Ignore pixels with any of these bits set
                               ) :
        afwMath::CandidateVisitor(),
        _imagePca(imagePca)
        {
            ;
        }
    
    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afwMath::SpatialCellCandidate *candidate) {
        PsfCandidate<MaskedImageT> *imCandidate = dynamic_cast<PsfCandidate<MaskedImageT> *>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }
        
        try {
            // Shift image to be centered in a pixel
            double const dx = afwImage::positionToIndex(imCandidate->getXCenter(), true).second;
            double const dy = afwImage::positionToIndex(imCandidate->getYCenter(), true).second;

            typename MaskedImageT::Ptr im =
                afwMath::offsetImage(*imCandidate->getImage(), -dx, -dy, "lanczos5");
            _imagePca->addImage(im, imCandidate->getSource().getPsfFlux());
        } catch(lsst::pex::exceptions::LengthErrorException &) {
            return;
        }
    }
private:
    afwImage::ImagePca<MaskedImageT> *_imagePca; // the ImagePca we're building
};

/************************************************************************************************************/
/// A class to pass around to all our PsfCandidates to count our candidates
template<typename PixelT>
class countVisitor : public afwMath::CandidateVisitor {
    typedef afwImage::MaskedImage<PixelT> MaskedImage;
public:
    explicit countVisitor() : afwMath::CandidateVisitor(), _n(0) {}
    
    void reset() {
        _n = 0;
    }
    
    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afwMath::SpatialCellCandidate *candidate) {
        PsfCandidate<MaskedImage> *imCandidate = dynamic_cast<PsfCandidate<MaskedImage> *>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }
        
        try {
            imCandidate->getImage();
        } catch(lsst::pex::exceptions::LengthErrorException &) {
            return;
        }
            
        ++_n;
    }
    
    // Return the number
    double getN() const { return _n; }
    
private:
    int mutable _n;                       // the desired number
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
        lsst::afw::geom::Extent2I const& dims, ///< Dimensions of image
        int const nEigenComponents,     ///< number of eigen components to keep; <= 0 => infty
        int const spatialOrder,         ///< Order of spatial variation (cf. afw::math::PolynomialFunction2)
        int const ksize,                ///< Size of generated Kernel images
        int const nStarPerCell,         ///< max no. of stars per cell; <= 0 => infty
        bool const constantWeight       ///< should each star have equal weight in the fit?
                                                                                                    )
{
    typedef typename afwImage::Image<PixelT> ImageT;
    typedef typename afwImage::MaskedImage<PixelT> MaskedImageT;
    //
    // Set the sizes for PsfCandidates made from either Images or MaskedImages
    //
    lsst::meas::algorithms::PsfCandidate<ImageT>::setWidth(ksize);
    lsst::meas::algorithms::PsfCandidate<ImageT>::setHeight(ksize);
    lsst::meas::algorithms::PsfCandidate<MaskedImageT>::setWidth(ksize);
    lsst::meas::algorithms::PsfCandidate<MaskedImageT>::setHeight(ksize);

    afwImage::ImagePca<MaskedImageT> imagePca(constantWeight); // Here's the set of images we'll analyze

    SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
    psfCells.visitCandidates(&importStarVisitor, nStarPerCell);

    //
    // Do a PCA decomposition of those PSF candidates.
    //
    // We have "gappy" data;  in other words we don't want to include any pixels with INTRP set
    //
    int niter = 10;                     // number of iterations of updateBadPixels
    double deltaLim = 10.0;             // acceptable value of delta, the max change due to updateBadPixels
    lsst::afw::image::MaskPixel const INTRP = afwImage::Mask<>::getPlaneBitMask("INTRP");
    lsst::afw::image::MaskPixel const CR = afwImage::Mask<>::getPlaneBitMask("CR");
    
    for (int i = 0; i != niter; ++i) {
        int const ncomp = (i == 0) ? 0 :
            ((nEigenComponents == 0) ? imagePca.getEigenImages().size() : nEigenComponents);
        double delta = imagePca.updateBadPixels(INTRP | CR, ncomp);
        if (i > 0 && delta < deltaLim) {
            break;
        }
        
        imagePca.analyze();
    }
    
    std::vector<typename MaskedImageT::Ptr> eigenImages = imagePca.getEigenImages();
    std::vector<double> eigenValues = imagePca.getEigenValues();
    int const nEigen = static_cast<int>(eigenValues.size());
    
    int const ncomp = (nEigenComponents <= 0 || nEigen < nEigenComponents) ? nEigen : nEigenComponents;
    //
    // Set the background level of the components to 0.0 to avoid coupling variable background
    // levels to the form of the Psf.  More precisely, we calculate the median of an outer "annulus"
    // of width bkg_width
    //
    int bkg_border = 2;
    for (int k = 0; k != ncomp; ++k) {
        ImageT const& im = *eigenImages[k]->getImage();
        
        if (bkg_border > im.getWidth()) {
            bkg_border = im.getWidth();
        }
        if (bkg_border > im.getHeight()) {
            bkg_border = im.getHeight();
        }

        double sum = 0;
        int n = 0;
        // Bottom and Top borders
        for (int i = 0; i != bkg_border; ++i) {
            typename ImageT::const_x_iterator
                ptrB = im.row_begin(i), ptrT = im.row_begin(im.getHeight() - i - 1);
            for (int j = 0; j != im.getWidth(); ++j, ++ptrB, ++ptrT) {
                sum += *ptrB + *ptrT;
                n += 2;
            }
        }
        for (int i = bkg_border; i < im.getHeight() - bkg_border; ++i) {
            // Left and Right borders
            typename ImageT::const_x_iterator
                ptrL = im.row_begin(i), ptrR = im.row_begin(im.getWidth() - i - 1);
            for (int j = 0; j != bkg_border; ++j, ++ptrL) {
                sum += *ptrL + *ptrR;
                n += 2;
            }
        }
        sum /= 2*(bkg_border*im.getWidth() + bkg_border*(im.getHeight() - 2*bkg_border));

        *eigenImages[k] -= sum;
    }
    //
    // Now build our LinearCombinationKernel; build the lists of basis functions
    // and spatial variation, then assemble the Kernel
    //
    afwMath::KernelList  kernelList;
    std::vector<afwMath::Kernel::SpatialFunctionPtr> spatialFunctionList;
    afwGeom::Box2D const range(afwGeom::Point2D(0.0, 0.0), afwGeom::Extent2D(dims.getX(), dims.getY()));

    for (int i = 0; i != ncomp; ++i) {
        kernelList.push_back(afwMath::Kernel::Ptr(new afwMath::FixedKernel(
                                      afwImage::Image<afwMath::Kernel::Pixel>(*eigenImages[i]->getImage(),true)
                                                                          )));

        afwMath::Kernel::SpatialFunctionPtr
            spatialFunction(new afwMath::Chebyshev1Function2<double>(spatialOrder, range));
        spatialFunction->setParameter(0, 1.0); // the constant term; all others are 0
        spatialFunctionList.push_back(spatialFunction);
    }

    afwMath::LinearCombinationKernel::Ptr
        psf(new afwMath::LinearCombinationKernel(kernelList, spatialFunctionList));

    return std::make_pair(psf, eigenValues);
}

/************************************************************************************************************/
/**
 * Count the number of candidates in use
 */
template<typename PixelT>
int countPsfCandidates(afwMath::SpatialCellSet const& psfCells,
                       int const nStarPerCell)
{
    countVisitor<PixelT> counter;
    psfCells.visitCandidates(&counter, nStarPerCell);

    return counter.getN();    
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
          DataImageT const& data,       // the data to fit
          double lambda = 0.0,          // floor for variance is lambda*data
          int const id=-1               // ID for this object; useful in debugging
         ) {
    assert(data.getDimensions() == mImage.getDimensions());
    assert(id == id);
    int const DETECTED = afwImage::Mask<>::getPlaneBitMask("DETECTED");

    double sumMM = 0.0, sumMD = 0.0, sumDD = 0.0; // sums of model*model/variance etc.
    int npix = 0;                                 // number of pixels used to evaluate chi^2
    for (int y = 0; y != data.getHeight(); ++y) {
        typename ModelImageT::x_iterator mptr = mImage.row_begin(y);
        for (typename DataImageT::x_iterator ptr = data.row_begin(y), end = data.row_end(y);
             ptr != end; ++ptr, ++mptr) {
            double const m = (*mptr)[0];       // value of model
            double const d = ptr.image();      // value of data
            double const var = ptr.variance() + lambda*d; // data's variance
            if (!(ptr.mask() & DETECTED)) {
                continue;
            }
            if (var != 0.0) {                  // assume variance == 0 => infinity XXX
                double const iVar = 1.0/var;
                npix++;
                sumMM += m*m*iVar;
                sumMD += m*d*iVar;
                sumDD += d*d*iVar;
            }
        }
    }
    
    if (npix == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "No good pixels");
    }
    if (sumMM == 0.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RangeErrorException, "sum(data*data)/var == 0");
    }

    double const amp = sumMD/sumMM;     // estimate of amplitude of model at this point            
    double const chi2 = (sumDD - 2*amp*sumMD + amp*amp*sumMM)/(npix - 1);

#if 0
    bool show = false;                  // Display the centre of the image; set from gdb
        
    if (show) {
        show = true;                    // you can jump here in gdb to set show if direct attempts fail
        int y = data.getHeight()/2;
        int x = data.getWidth()/2;
        int hsize = 2;
        printf("\ndata  ");
        for (int ii = -hsize; ii <= hsize; ++ii) {
            for (int jj = -hsize; jj <= hsize; ++jj) {
                printf("%7.1f ", data.at(x + jj, y - ii).image());
            }
            printf("  model  ");
            for (int jj = -hsize; jj <= hsize; ++jj) {
                printf("%7.1f ", amp*(*(mImage.at(x + jj, y - ii)))[0]);
            }
            printf("\n      ");
        }
        printf("%g  %.1f\n", amp, chi2);
    }
#endif
        
    return std::make_pair(chi2, amp);
}
}
    
/************************************************************************************************************/
/*
 * Fit for the spatial variation of the PSF parameters over the field
 */
/// A class to pass around to all our PsfCandidates to evaluate the PSF fit's X^2 
template<typename PixelT>
class evalChi2Visitor : public afwMath::CandidateVisitor {
    typedef afwImage::Image<PixelT> Image;
    typedef afwImage::MaskedImage<PixelT> MaskedImage;
    
    typedef afwImage::Image<afwMath::Kernel::Pixel> KImage;
public:
    explicit evalChi2Visitor(afwMath::Kernel const& kernel,
                             double lambda
                            ) :
        afwMath::CandidateVisitor(),
        _chi2(0.0), _kernel(kernel), _lambda(lambda),
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
        
        double const xcen = imCandidate->getSource().getXAstrom();
        double const ycen = imCandidate->getSource().getYAstrom();
        double const dx = afwImage::positionToIndex(xcen, true).second;
        double const dy = afwImage::positionToIndex(ycen, true).second;

        _kernel.computeImage(*_kImage, true, xcen, ycen);
        typename MaskedImage::Ptr data;
        try {
            data = afwMath::offsetImage(*imCandidate->getImage(), -dx, -dy);
        } catch(lsst::pex::exceptions::LengthErrorException &) {
            return;
        }
        
        try {
            std::pair<double, double> result = fitKernel(*_kImage, *data, _lambda,
                                                         imCandidate->getSource().getId());
            
            double dchi2 = result.first;      // chi^2 from this object
            double const amp = result.second; // estimate of amplitude of model at this point
            
            imCandidate->setChi2(dchi2);
            imCandidate->setAmplitude(amp);
            
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
    double _lambda;                  // floor for variance is _lambda*data
    typename KImage::Ptr mutable _kImage; // The Kernel at this point; a scratch copy
};
    
/********************************************************************************************************/
/**
 * Fit a Kernel's spatial variability from a set of stars
 */
// Set the Kernel's spatial parameters from a vector
void setSpatialParameters(afwMath::Kernel *kernel,
                          std::vector<double> const& coeffs
                         )
{
    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();

    assert (nComponents*nSpatialParams == static_cast<long>(coeffs.size()));

    std::vector<std::vector<double> > kCoeffs; // coefficients rearranged for Kernel
    kCoeffs.reserve(nComponents);
    for (int i = 0; i != nComponents; ++i) {
        kCoeffs.push_back(std::vector<double>(nSpatialParams));
        std::copy(coeffs.begin() + i*nSpatialParams,
                  coeffs.begin() + (i + 1)*nSpatialParams, kCoeffs[i].begin());
    }
    
    kernel->setSpatialParameters(kCoeffs);
}
    
/**
 * Fit a Kernel's spatial variability from a set of stars
 */
// Set the Kernel's spatial parameters from an Eigen::VectorXd
void setSpatialParameters(afwMath::Kernel *kernel,
                          Eigen::VectorXd const& vec
                         )
{
    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();

    assert (nComponents*nSpatialParams == vec.size());

    std::vector<std::vector<double> > kCoeffs; // coefficients rearranged for Kernel
    kCoeffs.reserve(nComponents);
    for (int i = 0; i != nComponents; ++i) {
        std::vector<double> spatialCoeffs(nSpatialParams);
        for (int j = 0; j != nSpatialParams; ++j) {
            spatialCoeffs[j] = vec[i*nSpatialParams + j];
        }
        kCoeffs.push_back(spatialCoeffs);
    }
    
    kernel->setSpatialParameters(kCoeffs);
}
    
//
// The object that minuit minimises
//
template<typename PixelT>
class MinimizeChi2 : public ROOT::Minuit2::FCNBase {
public:
    explicit MinimizeChi2(evalChi2Visitor<PixelT> & chi2Visitor,
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
                             _nSpatialParams(nSpatialParams) {}

/**
 * Error definition of the function. MINUIT defines Parameter errors as the
 * change in Parameter Value required to change the function Value by up. Normally,
 * for chisquared fits it is 1, and for negative log likelihood, its Value is 0.5.
 * If the user wants instead the 2-sigma errors for chisquared fits, it becomes 4,
 */
    double Up() const { return _errorDef; }
        
    // Evaluate our cost function (in this case chi^2)
    double operator()(const std::vector<double>& coeffs) const {
        setSpatialParameters(_kernel, coeffs);
        
        _psfCells.visitCandidates(&_chi2Visitor, _nStarPerCell);
        
        return _chi2Visitor.getValue();
    }
    
    void setErrorDef(double def) { _errorDef = def; }
private:
    double _errorDef;               // how much cost function has changed at the +- 1 error points
    
    evalChi2Visitor<PixelT>& _chi2Visitor;
    afwMath::Kernel *_kernel;
    afwMath::SpatialCellSet const& _psfCells;
    int _nStarPerCell;
    int _nComponents;
    int _nSpatialParams;
};
    
/************************************************************************************************************/
/**
 * Fit spatial kernel using full-nonlinear optimization to estimate candidate amplitudes
 */    
template<typename PixelT>
std::pair<bool, double>
fitSpatialKernelFromPsfCandidates(
        afwMath::Kernel *kernel,                 ///< the Kernel to fit
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        int const nStarPerCell,                  ///< max no. of stars per cell; <= 0 => infty
        double const tolerance,                  ///< Tolerance; how close chi^2 should be to true minimum
        double const lambda                      ///< floor for variance is lambda*data
                                 ) {
    typedef typename afwImage::Image<PixelT> Image;

    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();
    //
    // visitor that evaluates the chi^2 of the current fit
    //
    evalChi2Visitor<PixelT> getChi2(*kernel, lambda);
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
    ROOT::Minuit2::MnUserParameters fitPar;
    std::vector<std::string> paramNames;
    paramNames.reserve(nComponents*nSpatialParams);
    
    for (int i = 0, c = 0; c != nComponents; ++c) {
        coeffs[i] = 1;                  // the constant part of each spatial order
        for (int s = 0; s != nSpatialParams; ++s, ++i) {
            paramNames.push_back((boost::format("C%d:%d") % c % s).str());
            fitPar.Add(paramNames[i].c_str(), coeffs[i], stepSize[i]);
        }
    }
    fitPar.Fix("C0:0");
    //
    // Create the minuit object that knows how to minimise our functor
    //
    MinimizeChi2<PixelT> minimizerFunc(getChi2, kernel, psfCells, nStarPerCell, nComponents, nSpatialParams);

    double const errorDef = 1.0;       // use +- 1sigma errors
    minimizerFunc.setErrorDef(errorDef);
    //
    // tell minuit about it
    //    
    ROOT::Minuit2::MnMigrad migrad(minimizerFunc, fitPar);
    //
    // And let it loose
    //
    int maxFnCalls = 0;                 // i.e. unlimited
    ROOT::Minuit2::FunctionMinimum min =
        migrad(maxFnCalls, tolerance/(1e-4*errorDef)); // minuit uses 0.1*1e-3*tolerance*errorDef

    float minChi2 = min.Fval();
    bool const isValid = min.IsValid() && std::isfinite(minChi2);
    
    if (true || isValid) {              // calculate coeffs even in minuit is unhappy
        for (int i = 0; i != nComponents*nSpatialParams; ++i) {
            coeffs[i] = min.UserState().Value(i);
        }

        setSpatialParameters(kernel, coeffs);
    }

#if 0                                   // Estimate errors;  we don't really need this
    ROOT::Minuit2::MnMinos minos(minimizerFunc, min);
    for (int i = 0, c = 0; c != nComponents; ++c) {
        for (int s = 0; s != nSpatialParams; ++s, ++i) {
            char const *name = paramNames[i].c_str();
            printf("%s %g", name, min.UserState().Value(name));
            if (isValid && !fitPar.Parameter(fitPar.Index(name)).IsFixed()) {
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
    psfCells.visitAllCandidates(&getChi2, true);
    
    return std::make_pair(isValid, minChi2);
}
    
/************************************************************************************************************/
/**
 * Fit spatial kernel using approximate fluxes for candidates, and solving a linear system of equations
 */        
namespace {
/// A class to calculate the A and b matrices used to estimate the PSF's spatial structure
template<typename PixelT>
class FillABVisitor : public afwMath::CandidateVisitor {
    typedef afwImage::Image<PixelT> Image;
    typedef afwImage::MaskedImage<PixelT> MaskedImage;
    
    typedef afwImage::Image<afwMath::Kernel::Pixel> KImage;
public:
    explicit FillABVisitor(afwMath::LinearCombinationKernel const& kernel, // the Kernel we're fitting
                           double tau2=0.0                // floor to the per-candidate variance
                          ) :
        afwMath::CandidateVisitor(),
        _kernel(kernel),
        _tau2(tau2),
        _nSpatialParams(_kernel.getNSpatialParameters()),
        _nComponents(_kernel.getNKernelParameters()),
        _basisImgs(),
        _A(_nComponents*_nSpatialParams, _nComponents*_nSpatialParams),
        _b(_nComponents*_nSpatialParams),
        _basisDotBasis(_nComponents, _nComponents)
    {
        _basisImgs.resize(_nComponents);

        _A.setZero();
        _b.setZero();
        //
        // Get all the Kernel's components as Images
        //
        afwMath::KernelList const& kernels = _kernel.getKernelList(); // Kernel's components
        for (int i = 0; i != _nComponents; ++i) {
            _basisImgs[i] = typename KImage::Ptr(new KImage(kernels[i]->getDimensions()));
            kernels[i]->computeImage(*_basisImgs[i], false);
        }

        //
        // Calculate the inner products of the Kernel components once and for all
        //
        for (int i = 0; i != _nComponents; ++i) {
            for (int j = i; j != _nComponents; ++j) {
                _basisDotBasis(i, j) = _basisDotBasis(j, i) =
                    afwImage::innerProduct(*_basisImgs[i], *_basisImgs[j],
                                           PsfCandidate<MaskedImage>::getBorderWidth());
            }
        }
    }
    
    void reset() {}
    
    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afwMath::SpatialCellCandidate *candidate) {
        PsfCandidate<MaskedImage> *imCandidate = dynamic_cast<PsfCandidate<MaskedImage> *>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }

        typename MaskedImage::ConstPtr data;
        try {
            data = imCandidate->getImage();
        } catch(lsst::pex::exceptions::LengthErrorException &) {
            return;
        }
#if 0
        double const amp = imCandidate->getAmplitude();
#else
        /*
         * Estimate the amplitude based on the current basis functions.
         *
         * N.b. you have to be a little careful here.  Consider a PSF that is phi == (N0 + b*y*N1)/(1 + b*y)
         * where the amplitude of N0 and N1 is 1.0, so a star has profile I = A*(N0 + b*y*N1)/(1 + b*y)
         *
         * If we set the amplitude to be A = I(0)/phi(0) (i.e. the central value of the data and best-fit phi)
         * then the coefficient of N0 becomes 1/(1 + b*y) which makes the model non-linear in y.
         */
        double const xcen = imCandidate->getXCenter();
        double const ycen = imCandidate->getYCenter();

        double const dx = afwImage::positionToIndex(xcen, true).second;
        double const dy = afwImage::positionToIndex(ycen, true).second;
        double amp = 0.0;
        {
            std::pair<afwMath::Kernel::Ptr, std::pair<double, double> > ret =
                fitKernelToImage(_kernel, *data, afwGeom::PointD(xcen, ycen));
            
            afwMath::Kernel::Ptr bestFitKernel = ret.first;
            amp = ret.second.first;
            {
                afwImage::Image<afwMath::Kernel::Pixel> kImage(bestFitKernel->getDimensions());
                bestFitKernel->computeImage(kImage, false);
#define PRINT 0
#if PRINT
                std::cout << "Amp = " << imCandidate->getYCenter() << " " << amp << " ";
#endif
                amp = afwMath::makeStatistics(kImage, afwMath::MAX).getValue();
                

#if PRINT
                std::cout << amp << " data " <<
                    afwMath::makeStatistics(*data->getImage(), afwMath::SUM).getValue();
#endif
                {
                    std::vector<double> params = bestFitKernel->getKernelParameters();
                    std::vector<double> ksum = _kernel.getKernelSumList();

                    double sum = 0.0;
                    for (unsigned int i = 0; i != params.size(); ++i) {
                        sum += params[i]*ksum[i];
                    }
                    amp /= sum;
                }
                //amp *= params[0];
            }
#if PRINT
            std::cout << " correctedAmp " << amp/1e4 << std::endl;
#endif
        }
#endif
        
        double const var = imCandidate->getVar();
        double const ivar = 1/(var + _tau2); // Allow for floor on variance

        // Spatial params of all the components
        std::vector<std::vector<double> > params(_nComponents);
        for (int ic = 0; ic != _nComponents; ++ic) {
            params[ic] = _kernel.getSpatialFunction(ic)->getDFuncDParameters(xcen, ycen);
        }

        for (int i = 0, ic = 0; ic != _nComponents; ++ic) {
            typename KImage::Ptr tmp = afwMath::offsetImage(*_basisImgs[ic], dx, dy);

            double const basisDotData = afwImage::innerProduct(*tmp, *data->getImage(),
                                                               PsfCandidate<MaskedImage>::getBorderWidth());
            for (int is = 0; is != _nSpatialParams; ++is, ++i) {
                _b(i) += ivar*params[ic][is]*basisDotData/amp;
                
                for (int j = i, jc = ic; jc != _nComponents; ++jc) {
                    for (int js = (i == j) ? is : 0; js != _nSpatialParams; ++js, ++j) {
                        _A(i, j) += ivar*params[ic][is]*params[jc][js]*_basisDotBasis(ic, jc);
                        _A(j, i) = _A(i, j); // could do this after _A is fully calculated
                    }
                }
            }
        }
    }

    Eigen::MatrixXd const& getA() const { return _A; }
    Eigen::VectorXd const& getB() const { return _b; }
    
private:
    afwMath::LinearCombinationKernel const& _kernel;  // the kernel
    double _tau2;                    // variance floor added in quadrature to true candidate variance
    int const _nSpatialParams;       // number of spatial parameters
    int const _nComponents;          // number of basis functions
    std::vector<typename KImage::Ptr> _basisImgs; // basis function images from _kernel
    Eigen::MatrixXd _A;              // We'll solve the matrix equation A x = b for the Kernel's coefficients
    Eigen::VectorXd _b;
    Eigen::MatrixXd _basisDotBasis;  // the inner products of the  Kernel components
};


/// A class to set the best-fit PSF amplitude for an object
template<typename PixelT>
class setAmplitudeVisitor : public afwMath::CandidateVisitor {
    typedef afwImage::MaskedImage<PixelT> MaskedImage;
public:
    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afwMath::SpatialCellCandidate *candidate) {
        PsfCandidate<MaskedImage> *imCandidate = dynamic_cast<PsfCandidate<MaskedImage> *>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }
        imCandidate->setAmplitude(afwMath::makeStatistics(*imCandidate->getImage()->getImage(),
                                                          afwMath::MAX).getValue());
    }
};

}

template<typename PixelT>
std::pair<bool, double>
fitSpatialKernelFromPsfCandidates(
        afwMath::Kernel *kernel,                 ///< the Kernel to fit
        afwMath::SpatialCellSet const& psfCells, ///< A SpatialCellSet containing PsfCandidates
        bool const doNonLinearFit,               ///< Use the full-up nonlinear fitter
        int const nStarPerCell,                  ///< max no. of stars per cell; <= 0 => infty
        double const tolerance,                   ///< Tolerance; how close chi^2 should be to true minimum
        double const lambda                       ///< floor for variance is lambda*data
                                 )
{
    if (doNonLinearFit) {
        return fitSpatialKernelFromPsfCandidates<PixelT>(kernel, psfCells, nStarPerCell, tolerance);
    }

    double const tau = 0;               // softening for errors

    afwMath::LinearCombinationKernel const* lcKernel =
        dynamic_cast<afwMath::LinearCombinationKernel const*>(kernel);
    if (!lcKernel) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                          "Failed to cast Kernel to LinearCombinationKernel while building spatial PSF model");
    }
#if 0
    //
    // Set the initial amplitudes of all our candidates
    //
    setAmplitudeVisitor<PixelT> setAmplitude;
    psfCells.visitAllCandidates(&setAmplitude, true);
#endif
    //
    // visitor that fills out the A and b matrices (we'll solve A x = b for the coeffs, x)
    //
    FillABVisitor<PixelT> getAB(*lcKernel, tau);
    //
    // Actually visit all our candidates
    //
    psfCells.visitCandidates(&getAB, nStarPerCell, true);
    //
    // Extract A and b, and solve Ax = b
    //
    Eigen::MatrixXd const& A = getAB.getA();
    Eigen::VectorXd const& b = getAB.getB();
    assert(b.size() > 1);               // eigen has/had problems with 1x1 matrices; fix me if we fail here
    Eigen::VectorXd x(b.size());
    A.svd().solve(b, &x);
#if 0
    std::cout << "x " << x.transpose() << std::endl;

    if (x.cols() >= 6) {
        for (int i = 0; i != 6; ++i) {
            double xcen = 25; double ycen = 35 + 35*i;
            std::cout << "x, y " << xcen << " , " << ycen << " b "
                      << (x[3] + xcen*x[4] + ycen*x[5])/(x[0] + xcen*x[1] + ycen*x[2]) << std::endl;
        }
    }
#endif

    setSpatialParameters(kernel, x);
    //
    // One time more through the Candidates setting their chi^2 values. We'll
    // do all the candidates this time, not just the first nStarPerCell
    //
    // visitor that evaluates the chi^2 of the current fit
    //
    evalChi2Visitor<PixelT> getChi2(*kernel, lambda);

    psfCells.visitAllCandidates(&getChi2, true);
    
    return std::make_pair(true, getChi2.getValue());
}

/************************************************************************************************************/
/**
 * Subtract a PSF from an image at a given position
 */
template<typename MaskedImageT>
double subtractPsf(afwDetection::Psf const& psf,      ///< the PSF to subtract
                   MaskedImageT *data,  ///< Image to subtract PSF from
                   double x,            ///< column position
                   double y             ///< row position
                  )
{
    if (lsst::utils::isnan(x + y)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //
    // Get Psf candidate
    //
    afwDetection::Psf::Image::Ptr kImage = psf.computeImage(afwGeom::PointD(x, y));
    //
    // Now find the proper sub-Image
    //
    afwGeom::BoxI bbox = kImage->getBBox(afwImage::LOCAL);
    bbox.shift(kImage->getXY0() - data->getXY0());
    
    typename MaskedImageT::Ptr subData(new MaskedImageT(*data, bbox, afwImage::LOCAL, false)); // a shallow copy
    //
    // Now we've got both; find the PSF's amplitude
    //
    double lambda = 0.0;                // floor for variance is lambda*data
    try {
        std::pair<double, double> result = fitKernel(*kImage, *subData, lambda);
        double const chi2 = result.first; // chi^2 for fit
        double const amp = result.second; // estimate of amplitude of model at this point
        //
        // Convert kImage to the proper type so that I can subtract it.
        //
        typename MaskedImageT::Image::Ptr
            kImageF(new typename MaskedImageT::Image(*kImage, true)); // of data's type

        *kImageF *= amp;
        *subData->getImage() -= *kImageF;
        
        return chi2;
    } catch(lsst::pex::exceptions::RangeErrorException &e) {
        LSST_EXCEPT_ADD(e, (boost::format("Object at (%.2f, %.2f)") % x % y).str());
        throw e;
    }
}

/************************************************************************************************************/
/**
 * Fit a LinearCombinationKernel to an Image, allowing the coefficients of the components to vary
 *
 * @return std::pair(coefficients, std::pair(kernels, center amplitude))
 */
template<typename Image>
std::pair<std::vector<double>, std::pair<afwMath::KernelList, std::vector<double> > >
fitKernelParamsToImage(
        afwMath::LinearCombinationKernel const& kernel, ///< the Kernel to fit
        Image const& image,                             ///< the image to be fit
        afwGeom::Point2D const& pos                     ///< the position of the object
                )
{
    afwMath::KernelList kernels = kernel.getKernelList();         // the Kernels that kernel adds together
    int const nKernel = kernels.size();
    std::vector<typename afwImage::Image<afwMath::Kernel::Pixel>::Ptr>
        kernelImages(nKernel);          // images of each Kernel in kernels

    if (nKernel == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                          "Your kernel must have at least one component");
    }
    /*
     * Go through all the kernels, get a copy centered at the desired sub-pixel position, and then
     * extract a subImage from the parent image at the same place
     */ 
    int x0 = 0, y0 = 0;                 // the XY0() point of the shifted Kernel basis functions
    int ctrX = 0, ctrY = 0;

    afwImage::Image<afwMath::Kernel::Pixel> scr(kernel.getDimensions());
    for (int i = 0; i != nKernel; ++i) {
        assert (!kernels[i]->isSpatiallyVarying());
        
        if (i == 0) {
            ctrX = kernel.getCtrX();
            ctrY = kernel.getCtrY();
        }

        kernels[i]->computeImage(scr, false);
        kernelImages[i] = afwMath::offsetImage(scr, pos[0] - ctrX, pos[1] - ctrY);
        
        if (i == 0) {
            x0 = kernelImages[i]->getX0();
            y0 = kernelImages[i]->getY0();
        }
    }

    afwGeom::BoxI bbox(kernelImages[0]->getBBox(afwImage::PARENT));
    // allow for image's origin
    bbox.shift(afwGeom::ExtentI(-image.getX0(), -image.getY0()));
    // shallow copy
    Image const& subImage(Image(image, bbox, afwImage::LOCAL, false));  

    /*
     * Solve the linear problem  subImage = sum x_i K_i + epsilon; we solve this for x_i by constructing the
     * normal equations, A x = b
     */
    Eigen::MatrixXd A(nKernel, nKernel);
    Eigen::VectorXd b(nKernel);

    for (int i = 0; i != nKernel; ++i) {
        b(i) = afwImage::innerProduct(*kernelImages[i], *subImage.getImage());

        for (int j = i; j != nKernel; ++j) {
            A(i, j) = A(j, i) = afwImage::innerProduct(*kernelImages[i], *kernelImages[j]);
        }
    }
    Eigen::VectorXd x(nKernel);

    if (nKernel == 1) {
        x(0) = b(0)/A(0, 0);
    } else {
        A.svd().solve(b, &x);
    }

    std::vector<double> params(nKernel);
    afwMath::KernelList newKernels(nKernel);
    std::vector<double> amplitudes(nKernel);
    for (int i = 0; i != nKernel; ++i) {
        afwMath::Kernel::Ptr newKernel(new afwMath::FixedKernel(*kernelImages[i]));
        newKernel->setCtrX(x0 + static_cast<int>(newKernel->getWidth()/2));
        newKernel->setCtrY(y0 + static_cast<int>(newKernel->getHeight()/2));

        params[i] = x[i];
        newKernels[i] = newKernel;
        amplitudes[i] = (*kernelImages[i])(ctrX, ctrY);
    }

    return std::make_pair(params, std::make_pair(newKernels, amplitudes));
}


/************************************************************************************************************/
/**
 * Fit a LinearCombinationKernel to an Image, allowing the coefficients of the components to vary
 *
 * @return std::pair(best-fit kernel, std::pair(amp, chi^2))
 */
template<typename Image>
std::pair<afwMath::Kernel::Ptr, std::pair<double, double> >
fitKernelToImage(
        afwMath::LinearCombinationKernel const& kernel, ///< the Kernel to fit
        Image const& image,                             ///< the image to be fit
        afwGeom::Point2D const& pos                     ///< the position of the object
                )
{
    std::pair<std::vector<double>, std::pair<afwMath::KernelList, std::vector<double> > > const fit = 
        fitKernelParamsToImage(kernel, image, pos);
    std::vector<double> params = fit.first;
    afwMath::KernelList kernels = fit.second.first;
    std::vector<double> amplitudes = fit.second.second;
    int const nKernel = params.size();
    assert(kernels.size() == static_cast<unsigned int>(nKernel));
    assert(amplitudes.size() == static_cast<unsigned int>(nKernel));

    double amp = 0.0;
    for (int i = 0; i != nKernel; ++i) {
        amp += params[i] * amplitudes[i];
    }

    afwMath::Kernel::Ptr outputKernel(new afwMath::LinearCombinationKernel(kernels, params));
    double chisq = 0.0;
    outputKernel->setCtrX(kernels[0]->getCtrX());
    outputKernel->setCtrY(kernels[0]->getCtrY());

    return std::make_pair(outputKernel, std::make_pair(amp, chisq));
}


/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
    typedef float Pixel;
    template class PsfCandidate<afwImage::MaskedImage<Pixel> >;

    template
    std::pair<afwMath::LinearCombinationKernel::Ptr, std::vector<double> >
    createKernelFromPsfCandidates<Pixel>(afwMath::SpatialCellSet const&, afwGeom::Extent2I const&,
                                         int const, int const, int const, int const, bool const);
    template
    int countPsfCandidates<Pixel>(afwMath::SpatialCellSet const&, int const);

    template
    std::pair<bool, double>
    fitSpatialKernelFromPsfCandidates<Pixel>(afwMath::Kernel *, afwMath::SpatialCellSet const&,
                                             int const, double const, double const);
    template
    std::pair<bool, double>
    fitSpatialKernelFromPsfCandidates<Pixel>(afwMath::Kernel *, afwMath::SpatialCellSet const&, bool const,
                                             int const, double const, double const);

    template
    double subtractPsf(afwDetection::Psf const&, afwImage::MaskedImage<float> *, double, double);

    template
    std::pair<std::vector<double>, std::pair<afwMath::KernelList, std::vector<double> > >
    fitKernelParamsToImage(afwMath::LinearCombinationKernel const&,
                     afwImage::MaskedImage<Pixel> const&, afwGeom::Point2D const&);

    template
    std::pair<afwMath::Kernel::Ptr, std::pair<double, double> >
    fitKernelToImage(afwMath::LinearCombinationKernel const&,
                     afwImage::MaskedImage<Pixel> const&, afwGeom::Point2D const&);
/// \endcond
}}}
