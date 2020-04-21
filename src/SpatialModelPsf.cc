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
#include <numeric>

#if !defined(DOXYGEN)
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#endif

#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "Eigen/SVD"

#include "lsst/geom/Box.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/math/SpatialCell.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/meas/algorithms/ImagePca.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"
#include "lsst/meas/algorithms/PsfCandidate.h"

namespace lsst {
namespace meas {
namespace algorithms {

namespace {

int const WARP_BUFFER(1);                      // Buffer (border) around kernel to prevent warp issues
std::string const WARP_ALGORITHM("lanczos5");  // Name of warping algorithm to use

// A class to pass around to all our PsfCandidates which builds the PcaImageSet
template <typename PixelT>
class SetPcaImageVisitor : public afw::math::CandidateVisitor {
    typedef afw::image::Image<PixelT> ImageT;
    typedef afw::image::MaskedImage<PixelT> MaskedImageT;
    typedef afw::image::Exposure<PixelT> ExposureT;

public:
    explicit SetPcaImageVisitor(PsfImagePca<MaskedImageT>* imagePca,  // Set of Images to initialise
                                unsigned int const mask = 0x0  // Ignore pixels with any of these bits set
                                )
            : afw::math::CandidateVisitor(), _imagePca(imagePca) {
        ;
    }

    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afw::math::SpatialCellCandidate* candidate) {
        PsfCandidate<PixelT>* imCandidate = dynamic_cast<PsfCandidate<PixelT>*>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicError,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }

        try {
            std::shared_ptr<MaskedImageT> im = imCandidate->getOffsetImage(WARP_ALGORITHM, WARP_BUFFER);

            // static int count = 0;
            // im->writeFits(str(boost::format("cand%03d.fits") % count));
            // count += 1;

            afw::math::StatisticsControl sctrl;
            sctrl.setNanSafe(false);

            if (!std::isfinite(
                        afw::math::makeStatistics(*im->getImage(), afw::math::MAX, sctrl).getValue())) {
                throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                                  str(boost::format("Image at %d, %d contains NaN") %
                                      imCandidate->getXCenter() % imCandidate->getYCenter()));
            }
            if (!std::isfinite(
                        afw::math::makeStatistics(*im->getVariance(), afw::math::MAX, sctrl).getValue())) {
                throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeError,
                                  str(boost::format("Variance of Image at %d, %d contains NaN") %
                                      imCandidate->getXCenter() % imCandidate->getYCenter()));
            }

            _imagePca->addImage(im, imCandidate->getSource()->getPsfInstFlux());
        } catch (lsst::pex::exceptions::LengthError&) {
            return;
        }
    }

private:
    PsfImagePca<MaskedImageT>* _imagePca;  // the ImagePca we're building
};

/************************************************************************************************************/
/// A class to pass around to all our PsfCandidates to count our candidates
template <typename PixelT>
class countVisitor : public afw::math::CandidateVisitor {
    typedef afw::image::MaskedImage<PixelT> MaskedImage;
    typedef afw::image::Exposure<PixelT> Exposure;

public:
    explicit countVisitor() : afw::math::CandidateVisitor(), _n(0) {}

    void reset() { _n = 0; }

    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afw::math::SpatialCellCandidate* candidate) {
        PsfCandidate<PixelT>* imCandidate = dynamic_cast<PsfCandidate<PixelT>*>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicError,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }

        try {
            imCandidate->getMaskedImage();
        } catch (lsst::pex::exceptions::LengthError&) {
            return;
        }

        ++_n;
    }

    // Return the number
    double getN() const { return _n; }

private:
    int mutable _n;  // the desired number
};

/// Offset a kernel so that its sub-pixel position corresponds to that of some target image
///
/// We place the kernel in an oversized image (dimensions expanded by WARP_BUFFER*2) and resample that,
/// so that edge effects from resampling are minimised.
template <typename ImageT>
std::vector<std::shared_ptr<ImageT>> offsetKernel(
        afw::math::LinearCombinationKernel const& kernel,  ///< the Kernel to offset
        float dx, float dy                                 ///< Offset to apply
        ) {
    afw::math::KernelList kernels = kernel.getKernelList();      // The Kernels that kernel adds together
    unsigned int const nKernel = kernels.size();                 // Number of kernel components
    std::vector<std::shared_ptr<ImageT>> kernelImages(nKernel);  // Images of each Kernel in kernels
    if (nKernel == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "Kernel has no components");
    }

    ImageT scratch(kernel.getDimensions());  // Buffered scratch space
    for (unsigned int i = 0; i != nKernel; ++i) {
        kernels[i]->computeImage(scratch, false);
        kernelImages[i] = afw::math::offsetImage(scratch, dx, dy, WARP_ALGORITHM, WARP_BUFFER);
    }

    return kernelImages;
}

}  // Anonymous namespace

/************************************************************************************************************/
/**
 * Return a Kernel pointer and a list of eigenvalues resulting from analysing the provided SpatialCellSet
 *
 * The Kernel is a LinearCombinationKernel of the first nEigenComponents eigenImages
 *
 * N.b. This is templated over the Pixel type of the science image
 */
template <typename PixelT>
std::pair<std::shared_ptr<afw::math::LinearCombinationKernel>, std::vector<double>>
createKernelFromPsfCandidates(
        afw::math::SpatialCellSet const& psfCells,  ///< A SpatialCellSet containing PsfCandidates
        lsst::geom::Extent2I const& dims,           ///< Dimensions of image
        lsst::geom::Point2I const& xy0,             ///< Origin of image
        int const nEigenComponents,                 ///< number of eigen components to keep; <= 0 => infty
        int const spatialOrder,     ///< Order of spatial variation (cf. afw::math::PolynomialFunction2)
        int const ksize,            ///< Size of generated Kernel images
        int const nStarPerCell,     ///< max no. of stars per cell; <= 0 => infty
        bool const constantWeight,  ///< should each star have equal weight in the fit?
        int const border            ///< Border size for background subtraction
        ) {
    typedef typename afw::image::Image<PixelT> ImageT;
    typedef typename afw::image::MaskedImage<PixelT> MaskedImageT;

    //
    // Set the sizes for PsfCandidates made from either Images or MaskedImages
    //
    // lsst::meas::algorithms::PsfCandidate<ImageT>::setWidth(ksize);
    // lsst::meas::algorithms::PsfCandidate<ImageT>::setHeight(ksize);
    // lsst::meas::algorithms::PsfCandidate<MaskedImageT>::setWidth(ksize);
    // lsst::meas::algorithms::PsfCandidate<MaskedImageT>::setHeight(ksize);
    lsst::meas::algorithms::PsfCandidate<PixelT>::setWidth(ksize);
    lsst::meas::algorithms::PsfCandidate<PixelT>::setHeight(ksize);

    PsfImagePca<MaskedImageT> imagePca(constantWeight, border);  // Here's the set of images we'll analyze

    {
        SetPcaImageVisitor<PixelT> importStarVisitor(&imagePca);
        bool const ignoreExceptions = true;
        psfCells.visitCandidates(&importStarVisitor, nStarPerCell, ignoreExceptions);
    }

    //
    // Do a PCA decomposition of those PSF candidates.
    //
    // We have "gappy" data;  in other words we don't want to include any pixels with INTRP set
    //
    int niter = 10;          // number of iterations of updateBadPixels
    double deltaLim = 10.0;  // acceptable value of delta, the max change due to updateBadPixels
    lsst::afw::image::MaskPixel const BAD = afw::image::Mask<>::getPlaneBitMask("BAD");
    lsst::afw::image::MaskPixel const CR = afw::image::Mask<>::getPlaneBitMask("CR");
    lsst::afw::image::MaskPixel const INTRP = afw::image::Mask<>::getPlaneBitMask("INTRP");

    for (int i = 0; i != niter; ++i) {
        int const ncomp =
                (i == 0) ? 0
                         : ((nEigenComponents == 0) ? imagePca.getEigenImages().size() : nEigenComponents);
        double delta = imagePca.updateBadPixels(BAD | CR | INTRP, ncomp);
        if (i > 0 && delta < deltaLim) {
            break;
        }

        imagePca.analyze();
    }

    std::vector<std::shared_ptr<MaskedImageT>> eigenImages = imagePca.getEigenImages();
    std::vector<double> eigenValues = imagePca.getEigenValues();
    int const nEigen = static_cast<int>(eigenValues.size());

    int const ncomp = (nEigenComponents <= 0 || nEigen < nEigenComponents) ? nEigen : nEigenComponents;
    //
    // Set the background level of the components to 0.0 to avoid coupling variable background
    // levels to the form of the Psf.  More precisely, we calculate the mean of an outer "annulus"
    // of width bkg_border
    //
    for (int k = 0; k != ncomp; ++k) {
        ImageT const& im = *eigenImages[k]->getImage();

        int bkg_border = 2;
        if (bkg_border > im.getWidth()) {
            bkg_border = im.getWidth() / 2;
        }
        if (bkg_border > im.getHeight()) {
            bkg_border = im.getHeight() / 2;
        }

        double sum = 0;
        // Bottom and Top borders
        for (int i = 0; i != bkg_border; ++i) {
            typename ImageT::const_x_iterator ptrB = im.row_begin(i),
                                              ptrT = im.row_begin(im.getHeight() - i - 1);
            for (int j = 0; j != im.getWidth(); ++j, ++ptrB, ++ptrT) {
                sum += *ptrB + *ptrT;
            }
        }
        for (int i = bkg_border; i < im.getHeight() - bkg_border; ++i) {
            // Left and Right borders
            typename ImageT::const_x_iterator ptrL = im.row_begin(i),
                                              ptrR = im.row_begin(i) + im.getWidth() - bkg_border;
            for (int j = 0; j != bkg_border; ++j, ++ptrL, ++ptrR) {
                sum += *ptrL + *ptrR;
            }
        }
        sum /= 2 * (bkg_border * im.getWidth() + bkg_border * (im.getHeight() - 2 * bkg_border));

        *eigenImages[k] -= sum;
    }
    //
    // Now build our LinearCombinationKernel; build the lists of basis functions
    // and spatial variation, then assemble the Kernel
    //
    afw::math::KernelList kernelList;
    std::vector<afw::math::Kernel::SpatialFunctionPtr> spatialFunctionList;
    geom::Box2D const range = geom::Box2D(geom::Point2D(xy0), geom::Extent2D(dims));

    for (int i = 0; i != ncomp; ++i) {
        {
            // Enforce unit sum for kernel by construction
            // Zeroth component has unit sum
            // Other components have zero sum by normalising and then subtracting the zeroth component
            ImageT& image = *eigenImages[i]->getImage();
            double sum = std::accumulate(image.begin(true), image.end(true), 0.0);
            if (i == 0) {
                image /= sum;
            } else {
                for (typename ImageT::fast_iterator ptr0 = eigenImages[0]->getImage()->begin(true),
                                                    ptr1 = image.begin(true), end = image.end(true);
                     ptr1 != end; ++ptr0, ++ptr1) {
                    *ptr1 = *ptr1 / sum - *ptr0;
                }
            }
        }

        kernelList.push_back(std::shared_ptr<afw::math::Kernel>(new afw::math::FixedKernel(
                afw::image::Image<afw::math::Kernel::Pixel>(*eigenImages[i]->getImage(), true))));

        afw::math::Kernel::SpatialFunctionPtr
        //            spatialFunction(new afw::math::PolynomialFunction2<double>(spatialOrder));
        spatialFunction(new afw::math::Chebyshev1Function2<double>(spatialOrder, range));
        spatialFunction->setParameter(0, 1.0);  // the constant term; all others are 0
        spatialFunctionList.push_back(spatialFunction);
    }

    std::shared_ptr<afw::math::LinearCombinationKernel> psf(
            new afw::math::LinearCombinationKernel(kernelList, spatialFunctionList));

    return std::make_pair(psf, eigenValues);
}

/************************************************************************************************************/
/**
 * Count the number of candidates in use
 */
template <typename PixelT>
int countPsfCandidates(afw::math::SpatialCellSet const& psfCells, int const nStarPerCell) {
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
template <typename ModelImageT, typename DataImageT>
std::pair<double, double> fitKernel(ModelImageT const& mImage,  // The model image at this point
                                    DataImageT const& data,     // the data to fit
                                    double lambda = 0.0,        // floor for variance is lambda*data
                                    bool detected = true,       // only fit DETECTED pixels?
                                    int const id = -1           // ID for this object; useful in debugging
                                    ) {
    assert(data.getDimensions() == mImage.getDimensions());
    assert(id == id);
    int const DETECTED = afw::image::Mask<>::getPlaneBitMask("DETECTED");
    int const BAD = afw::image::Mask<>::getPlaneBitMask("CR") | afw::image::Mask<>::getPlaneBitMask("BAD");

    double sumMM = 0.0, sumMD = 0.0, sumDD = 0.0;  // sums of model*model/variance etc.
    int npix = 0;                                  // number of pixels used to evaluate chi^2
    for (int y = 0; y != data.getHeight(); ++y) {
        typename ModelImageT::x_iterator mptr = mImage.row_begin(y);
        for (typename DataImageT::x_iterator ptr = data.row_begin(y), end = data.row_end(y); ptr != end;
             ++ptr, ++mptr) {
            double const m = (*mptr)[0];                     // value of model
            double const d = ptr.image();                    // value of data
            double const var = ptr.variance() + lambda * d;  // data's variance
            if (detected && !(ptr.mask() & DETECTED)) {
                continue;
            }
            if (ptr.mask() & BAD) {
                continue;
            }
            if (var != 0.0) {  // assume variance == 0 => infinity XXX
                double const iVar = 1.0 / var;
                npix++;
                sumMM += m * m * iVar;
                sumMD += m * d * iVar;
                sumDD += d * d * iVar;
            }
        }
    }

    if (npix == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RangeError, "No good pixels");
    }
    if (sumMM == 0.0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RangeError, "sum(data*data)/var == 0");
    }

    double const amp = sumMD / sumMM;  // estimate of amplitude of model at this point
    double const chi2 = (sumDD - 2 * amp * sumMD + amp * amp * sumMM) / (npix - 1);

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
}  // namespace

/************************************************************************************************************/
/*
 * Fit for the spatial variation of the PSF parameters over the field
 */
/// A class to pass around to all our PsfCandidates to evaluate the PSF fit's X^2
template <typename PixelT>
class evalChi2Visitor : public afw::math::CandidateVisitor {
    typedef afw::image::Image<PixelT> Image;
    typedef afw::image::MaskedImage<PixelT> MaskedImage;
    typedef afw::image::Exposure<PixelT> Exposure;

    typedef afw::image::Image<afw::math::Kernel::Pixel> KImage;

public:
    explicit evalChi2Visitor(afw::math::Kernel const& kernel, double lambda)
            : afw::math::CandidateVisitor(),
              _chi2(0.0),
              _kernel(kernel),
              _lambda(lambda),
              _kImage(std::shared_ptr<KImage>(new KImage(kernel.getDimensions()))) {}

    void reset() { _chi2 = 0.0; }

    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afw::math::SpatialCellCandidate* candidate) {
        PsfCandidate<PixelT>* imCandidate = dynamic_cast<PsfCandidate<PixelT>*>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicError,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }

        double const xcen = imCandidate->getSource()->getX();
        double const ycen = imCandidate->getSource()->getY();

        _kernel.computeImage(*_kImage, true, xcen, ycen);
        std::shared_ptr<MaskedImage const> data;
        try {
            data = imCandidate->getOffsetImage(WARP_ALGORITHM, WARP_BUFFER);
        } catch (lsst::pex::exceptions::LengthError&) {
            return;
        }

        try {
            std::pair<double, double> result =
                    fitKernel(*_kImage, *data, _lambda, false, imCandidate->getSource()->getId());

            double dchi2 = result.first;       // chi^2 from this object
            double const amp = result.second;  // estimate of amplitude of model at this point

            imCandidate->setChi2(dchi2);
            imCandidate->setAmplitude(amp);

            _chi2 += dchi2;
        } catch (lsst::pex::exceptions::RangeError& e) {
            imCandidate->setStatus(afw::math::SpatialCellCandidate::BAD);
            imCandidate->setChi2(std::numeric_limits<double>::quiet_NaN());
            imCandidate->setAmplitude(std::numeric_limits<double>::quiet_NaN());
        }
    }

    // Return the computed chi^2
    double getValue() const { return _chi2; }

private:
    double mutable _chi2;                     // the desired chi^2
    afw::math::Kernel const& _kernel;         // the kernel
    double _lambda;                           // floor for variance is _lambda*data
    std::shared_ptr<KImage> mutable _kImage;  // The Kernel at this point; a scratch copy
};

/********************************************************************************************************/
/**
 * Fit a Kernel's spatial variability from a set of stars
 */
// Set the Kernel's spatial parameters from a vector
void setSpatialParameters(afw::math::Kernel* kernel, std::vector<double> const& coeffs) {
    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();

    assert(nComponents * nSpatialParams == static_cast<long>(coeffs.size()));

    std::vector<std::vector<double>> kCoeffs;  // coefficients rearranged for Kernel
    kCoeffs.reserve(nComponents);
    for (int i = 0; i != nComponents; ++i) {
        kCoeffs.push_back(std::vector<double>(nSpatialParams));
        std::copy(coeffs.begin() + i * nSpatialParams, coeffs.begin() + (i + 1) * nSpatialParams,
                  kCoeffs[i].begin());
    }

    kernel->setSpatialParameters(kCoeffs);
}

/**
 * Fit a Kernel's spatial variability from a set of stars
 */
// Set the Kernel's spatial parameters from an Eigen::VectorXd
void setSpatialParameters(afw::math::Kernel* kernel, Eigen::VectorXd const& vec) {
    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();

    assert(nComponents * nSpatialParams == vec.size());

    std::vector<std::vector<double>> kCoeffs;  // coefficients rearranged for Kernel
    kCoeffs.reserve(nComponents);
    for (int i = 0; i != nComponents; ++i) {
        std::vector<double> spatialCoeffs(nSpatialParams);
        for (int j = 0; j != nSpatialParams; ++j) {
            spatialCoeffs[j] = vec[i * nSpatialParams + j];
        }
        kCoeffs.push_back(spatialCoeffs);
    }

    kernel->setSpatialParameters(kCoeffs);
}

//
// The object that minuit minimises
//
template <typename PixelT>
class MinimizeChi2 : public ROOT::Minuit2::FCNBase {
public:
    explicit MinimizeChi2(evalChi2Visitor<PixelT>& chi2Visitor, afw::math::Kernel* kernel,
                          afw::math::SpatialCellSet const& psfCells, int nStarPerCell, int nComponents,
                          int nSpatialParams)
            : _errorDef(1.0),
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
    double _errorDef;  // how much cost function has changed at the +- 1 error points

    evalChi2Visitor<PixelT>& _chi2Visitor;
    afw::math::Kernel* _kernel;
    afw::math::SpatialCellSet const& _psfCells;
    int _nStarPerCell;
    int _nComponents;
    int _nSpatialParams;
};

/************************************************************************************************************/
/**
 * Fit spatial kernel using full-nonlinear optimization to estimate candidate amplitudes
 */
template <typename PixelT>
std::pair<bool, double> fitSpatialKernelFromPsfCandidates(
        afw::math::Kernel* kernel,                  ///< the Kernel to fit
        afw::math::SpatialCellSet const& psfCells,  ///< A SpatialCellSet containing PsfCandidates
        int const nStarPerCell,                     ///< max no. of stars per cell; <= 0 => infty
        double const tolerance,                     ///< Tolerance; how close chi^2 should be to true minimum
        double const lambda                         ///< floor for variance is lambda*data
        ) {
    int const nComponents = kernel->getNKernelParameters();
    int const nSpatialParams = kernel->getNSpatialParameters();
    //
    // visitor that evaluates the chi^2 of the current fit
    //
    evalChi2Visitor<PixelT> getChi2(*kernel, lambda);
    //
    // We have to unpack the Kernel coefficients into a linear array, coeffs
    //
    std::vector<double> coeffs;  // The coefficients we want to fit
    coeffs.assign(nComponents * nSpatialParams, 0.0);

    std::vector<double> stepSize;  // step sizes
    stepSize.assign(nComponents * nSpatialParams, 100);
    //
    // Translate that into minuit's language
    //
    ROOT::Minuit2::MnUserParameters fitPar;
    std::vector<std::string> paramNames;
    paramNames.reserve(nComponents * nSpatialParams);

    for (int i = 0, c = 0; c != nComponents; ++c) {
        coeffs[i] = 1;  // the constant part of each spatial order
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

    double const errorDef = 1.0;  // use +- 1sigma errors
    minimizerFunc.setErrorDef(errorDef);
    //
    // tell minuit about it
    //
    ROOT::Minuit2::MnMigrad migrad(minimizerFunc, fitPar);
    //
    // And let it loose
    //
    int maxFnCalls = 0;  // i.e. unlimited
    ROOT::Minuit2::FunctionMinimum min =
            migrad(maxFnCalls, tolerance / (1e-4 * errorDef));  // minuit uses 0.1*1e-3*tolerance*errorDef

    float minChi2 = min.Fval();
    bool const isValid = min.IsValid() && std::isfinite(minChi2);

    if (true || isValid) {  // calculate coeffs even in minuit is unhappy
        for (int i = 0; i != nComponents * nSpatialParams; ++i) {
            coeffs[i] = min.UserState().Value(i);
        }

        setSpatialParameters(kernel, coeffs);
    }

#if 0  // Estimate errors;  we don't really need this
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
///
/// Given a set of kernels, and postage stamps of stars, we want to generate a PSF:
///
/// PSF(u,v;x,y) = K_0(u,v) + sum_(i>0) a_i F_i(x,y) K_i(u,v)
///
/// where K_0 = k_0(u,v) / sum_(u,v) k_0(u,v)
/// and K_(i>0) = k_i(u,v) / sum_(u,v) k_i(u,v) - K_0(u,v)
///
/// The K_i, i > 0 have zero sum, while K_0 has unit sum, so the PSF sum will
/// always be unity by construction.  This is basically the Alard
/// (2000A&AS..144..363A) technique.
///
/// The kernels provided to us here (through the 'kernel' in the constructor)
/// are the K_i; the conversion from k_i to K_i is done elsewhere (e.g.,
/// createKernelFromPsfCandidates).
///
/// Then the problem may be expressed as the matrix problem, Ax = b, where:
///
/// A_(i,j) = sum_(x,y,u,v) F_i(x,y) K_i(u,v) F_j(x,y) K_j(u,v)
///
/// b_i = sum_(x,y,u,v) [D(u,v) - K_0(u,v)] F_i(x,y) K_i(u,v)
///
/// Here, the sum_(x,y,u,v) means over all the postage stamps (u,v) of all the
/// stars (different x,y positions).  We take x,y for the star as the x,y for
/// the entire postage stamp, which greatly reduces the computation burden of
/// recalculating the polynomial for every pixel.  This assumes that the spatial
/// variation is smooth and gradual.
///
/// Note that because the 0th component has no spatial variation (in a formal
/// sense; its spatial variation is accomplished through its being subtracted
/// from the other components), the 'A' matrix and 'b' vector have have
/// dimensions (nComponents-1)*nSpatialParams, rather than
/// nComponents*nSpatialParams.  This affects the bounds of some of the
/// iterations, below.
///
template <typename PixelT>
class FillABVisitor : public afw::math::CandidateVisitor {
    typedef afw::image::Image<PixelT> Image;
    typedef afw::image::MaskedImage<PixelT> MaskedImage;
    typedef afw::image::Exposure<PixelT> Exposure;

    typedef afw::image::Image<afw::math::Kernel::Pixel> KImage;

public:
    explicit FillABVisitor(afw::math::LinearCombinationKernel const& kernel,  // the Kernel we're fitting
                           double tau2 = 0.0  // floor to the per-candidate variance
                           )
            : afw::math::CandidateVisitor(),
              _kernel(kernel),
              _tau2(tau2),
              _nSpatialParams(_kernel.getNSpatialParameters()),
              _nComponents(_kernel.getNKernelParameters()),
              _basisImgs(),
              _A((_nComponents - 1) * _nSpatialParams, (_nComponents - 1) * _nSpatialParams),
              _b((_nComponents - 1) * _nSpatialParams),
              _basisDotBasis(_nComponents, _nComponents) {
        _basisImgs.resize(_nComponents);

        _A.setZero();
        _b.setZero();
        //
        // Get all the Kernel's components as Images
        //
        afw::math::KernelList const& kernels = _kernel.getKernelList();  // Kernel's components
        for (int i = 0; i != _nComponents; ++i) {
            _basisImgs[i] = std::shared_ptr<KImage>(new KImage(kernels[i]->getDimensions()));
            kernels[i]->computeImage(*_basisImgs[i], false);
        }

        //
        // Calculate the inner products of the Kernel components once and for all
        //
        for (int i = 1; i != _nComponents; ++i) {  // Don't need 0th component
            for (int j = i; j != _nComponents; ++j) {
                _basisDotBasis(i, j) = _basisDotBasis(j, i) = afw::image::innerProduct(
                        *_basisImgs[i], *_basisImgs[j], PsfCandidate<PixelT>::getBorderWidth());
            }
        }
    }

    void reset() {}

    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afw::math::SpatialCellCandidate* candidate) {
        PsfCandidate<PixelT>* imCandidate = dynamic_cast<PsfCandidate<PixelT>*>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicError,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }

        CONST_PTR(MaskedImage) data;
        try {
            data = imCandidate->getMaskedImage(_kernel.getWidth(), _kernel.getHeight());
        } catch (lsst::pex::exceptions::LengthError&) {
            return;
        }
        double const xcen = imCandidate->getXCenter();
        double const ycen = imCandidate->getYCenter();
        double const dx = afw::image::positionToIndex(xcen, true).second;
        double const dy = afw::image::positionToIndex(ycen, true).second;

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
        std::pair<std::shared_ptr<afw::math::Kernel>, std::pair<double, double>> ret =
                fitKernelToImage(_kernel, *data, geom::Point2D(xcen, ycen));
        double const amp = ret.second.first;
#endif

        double const var = imCandidate->getVar();
        double const ivar = 1 / (var + _tau2);  // Allow for floor on variance

        // Spatial params of all the components
        std::vector<std::vector<double>> params(_nComponents);
        for (int ic = 1; ic != _nComponents; ++ic) {  // Don't need params[0]
            params[ic] = _kernel.getSpatialFunction(ic)->getDFuncDParameters(xcen, ycen);
        }

        std::vector<std::shared_ptr<KImage>> basisImages = offsetKernel<KImage>(_kernel, dx, dy);

        // Prepare values for basis dot data
        // Scale data and subtract 0th component as part of unit kernel sum construction
        std::shared_ptr<Image> dataImage(new Image(*data->getImage(), true));
        typename KImage::fast_iterator bPtr = basisImages[0]->begin(true);
        for (typename Image::fast_iterator dPtr = dataImage->begin(true), end = dataImage->end(true);
             dPtr != end; ++dPtr, ++bPtr) {
            *dPtr = *dPtr / amp - *bPtr;
        }

        for (int i = 0, ic = 1; ic != _nComponents; ++ic) {  // Don't need 0th component now
            double const basisDotData = afw::image::innerProduct(*basisImages[ic], *dataImage,
                                                                 PsfCandidate<PixelT>::getBorderWidth());
            for (int is = 0; is != _nSpatialParams; ++is, ++i) {
                _b(i) += ivar * params[ic][is] * basisDotData;

                for (int j = i, jc = ic; jc != _nComponents; ++jc) {
                    for (int js = (i == j) ? is : 0; js != _nSpatialParams; ++js, ++j) {
                        _A(i, j) += ivar * params[ic][is] * params[jc][js] * _basisDotBasis(ic, jc);
                        _A(j, i) = _A(i, j);  // could do this after _A is fully calculated
                    }
                }
            }
        }
    }

    Eigen::MatrixXd const& getA() const { return _A; }
    Eigen::VectorXd const& getB() const { return _b; }

private:
    afw::math::LinearCombinationKernel const& _kernel;  // the kernel
    double _tau2;               // variance floor added in quadrature to true candidate variance
    int const _nSpatialParams;  // number of spatial parameters
    int const _nComponents;     // number of basis functions
    std::vector<std::shared_ptr<KImage>> _basisImgs;  // basis function images from _kernel
    Eigen::MatrixXd _A;  // We'll solve the matrix equation A x = b for the Kernel's coefficients
    Eigen::VectorXd _b;
    Eigen::MatrixXd _basisDotBasis;  // the inner products of the  Kernel components
};

/// A class to set the best-fit PSF amplitude for an object
template <typename PixelT>
class setAmplitudeVisitor : public afw::math::CandidateVisitor {
    typedef afw::image::MaskedImage<PixelT> MaskedImage;
    typedef afw::image::Exposure<PixelT> Exposure;

public:
    // Called by SpatialCellSet::visitCandidates for each Candidate
    void processCandidate(afw::math::SpatialCellCandidate* candidate) {
        PsfCandidate<PixelT>* imCandidate = dynamic_cast<PsfCandidate<PixelT>*>(candidate);
        if (imCandidate == NULL) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LogicError,
                              "Failed to cast SpatialCellCandidate to PsfCandidate");
        }
        imCandidate->setAmplitude(
                afw::math::makeStatistics(*imCandidate->getMaskedImage()->getImage(), afw::math::MAX)
                        .getValue());
    }
};

}  // namespace

template <typename PixelT>
std::pair<bool, double> fitSpatialKernelFromPsfCandidates(
        afw::math::Kernel* kernel,                  ///< the Kernel to fit
        afw::math::SpatialCellSet const& psfCells,  ///< A SpatialCellSet containing PsfCandidates
        bool const doNonLinearFit,                  ///< Use the full-up nonlinear fitter
        int const nStarPerCell,                     ///< max no. of stars per cell; <= 0 => infty
        double const tolerance,                     ///< Tolerance; how close chi^2 should be to true minimum
        double const lambda                         ///< floor for variance is lambda*data
        ) {
    if (doNonLinearFit) {
        return fitSpatialKernelFromPsfCandidates<PixelT>(kernel, psfCells, nStarPerCell, tolerance);
    }

    double const tau = 0;  // softening for errors

    afw::math::LinearCombinationKernel const* lcKernel =
            dynamic_cast<afw::math::LinearCombinationKernel const*>(kernel);
    if (!lcKernel) {
        throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicError,
                "Failed to cast Kernel to LinearCombinationKernel while building spatial PSF model");
    }
#if 1
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
    Eigen::VectorXd x0(b.size());  // Solution to matrix problem

    switch (b.size()) {
        case 0:  // One candidate, no spatial variability
            break;
        case 1:  // eigen can't/won't handle 1x1 matrices
            x0(0) = b(0) / A(0, 0);
            break;
        default:
            x0 = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
            break;
    }
#if 0
    std::cout << "A " << A << std::endl;
    std::cout << "b " << b.transpose() << std::endl;
    std::cout << "x " << x.transpose() << std::endl;

    afw::image::Image<double> img(b.size(), b.size());
    for (int j = 0; j < b.size(); ++j) {
        for (int i = 0; i < b.size(); ++i) {
            img(i, j) = A(i, j);
        }
    }
    img.writeFits("a.fits");

    if (x.cols() >= 6) {
        for (int i = 0; i != 6; ++i) {
            double xcen = 25; double ycen = 35 + 35*i;
            std::cout << "x, y " << xcen << " , " << ycen << " b "
                      << (x[3] + xcen*x[4] + ycen*x[5])/(x[0] + xcen*x[1] + ycen*x[2]) << std::endl;
        }
    }
#endif

    // Generate kernel parameters (including 0th component) from matrix solution
    Eigen::VectorXd x(kernel->getNKernelParameters() * kernel->getNSpatialParameters());  // Kernel parameters
    x(0) = 1.0;
    std::fill(x.data() + 1, x.data() + kernel->getNSpatialParameters(), 0.0);
    std::copy(x0.data(), x0.data() + x0.size(), x.data() + kernel->getNSpatialParameters());

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
template <typename MaskedImageT>
double subtractPsf(afw::detection::Psf const& psf,  ///< the PSF to subtract
                   MaskedImageT* data,              ///< Image to subtract PSF from
                   double x,                        ///< column position
                   double y,                        ///< row position
                   double psfFlux                   ///< object's PSF flux (if not NaN)
                   ) {
    if (std::isnan(x + y)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //
    // Get Psf candidate
    //
    std::shared_ptr<afw::detection::Psf::Image> kImage = psf.computeImage(geom::PointD(x, y));

    //
    // Now find the proper sub-Image
    //
    geom::BoxI bbox = kImage->getBBox();

    std::shared_ptr<MaskedImageT> subData(
            new MaskedImageT(*data, bbox, afw::image::PARENT, false));  // shallow copy
    //
    // Now we've got both; find the PSF's amplitude
    //
    double lambda = 0.0;  // floor for variance is lambda*data
    try {
        double chi2;  // chi^2 for fit
        double amp;   // estimate of amplitude of model at this point

        if (std::isnan(psfFlux)) {
            std::pair<double, double> result = fitKernel(*kImage, *subData, lambda, true);
            chi2 = result.first;  // chi^2 for fit
            amp = result.second;  // estimate of amplitude of model at this point
        } else {
            chi2 = std::numeric_limits<double>::quiet_NaN();
            amp = psfFlux / afw::math::makeStatistics(*kImage, afw::math::SUM).getValue();
        }
        //
        // Convert kImage to the proper type so that I can subtract it.
        //
        std::shared_ptr<typename MaskedImageT::Image> kImageF(
                new typename MaskedImageT::Image(*kImage, true));  // of data's type

        *kImageF *= amp;
        *subData->getImage() -= *kImageF;

        return chi2;
    } catch (lsst::pex::exceptions::RangeError& e) {
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
template <typename Image>
std::pair<std::vector<double>, afw::math::KernelList> fitKernelParamsToImage(
        afw::math::LinearCombinationKernel const& kernel,  ///< the Kernel to fit
        Image const& image,                                ///< the image to be fit
        geom::Point2D const& pos                           ///< the position of the object
        ) {
    typedef afw::image::Image<afw::math::Kernel::Pixel> KernelT;

    afw::math::KernelList kernels = kernel.getKernelList();  // the Kernels that kernel adds together
    int const nKernel = kernels.size();

    if (nKernel == 0) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthError, "Your kernel must have at least one component");
    }

    /*
     * Go through all the kernels, get a copy centered at the desired sub-pixel position, and then
     * extract a subImage from the parent image at the same place
     */
    std::vector<std::shared_ptr<KernelT>> kernelImages = offsetKernel<KernelT>(kernel, pos[0], pos[1]);
    geom::BoxI bbox(kernelImages[0]->getBBox());
    Image const& subImage(Image(image, bbox, afw::image::PARENT, false));  // shallow copy

    /*
     * Solve the linear problem  subImage = sum x_i K_i + epsilon; we solve this for x_i by constructing the
     * normal equations, A x = b
     */
    Eigen::MatrixXd A(nKernel, nKernel);
    Eigen::VectorXd b(nKernel);

    for (int i = 0; i != nKernel; ++i) {
        b(i) = afw::image::innerProduct(*kernelImages[i], *subImage.getImage());

        for (int j = i; j != nKernel; ++j) {
            A(i, j) = A(j, i) = afw::image::innerProduct(*kernelImages[i], *kernelImages[j]);
        }
    }
    Eigen::VectorXd x(nKernel);

    if (nKernel == 1) {
        x(0) = b(0) / A(0, 0);
    } else {
        x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    }

    // the XY0() point of the shifted Kernel basis functions
    geom::Point2I const xy0 = kernelImages[0]->getXY0();

    afw::math::KernelList newKernels(nKernel);
    std::vector<double> params(nKernel);
    for (int i = 0; i != nKernel; ++i) {
        std::shared_ptr<afw::math::Kernel> newKernel(new afw::math::FixedKernel(*kernelImages[i]));
        newKernel->setCtr(xy0 + newKernel->getDimensions() / 2);

        params[i] = x[i];
        newKernels[i] = newKernel;
    }

    return std::make_pair(params, newKernels);
}

/************************************************************************************************************/
/**
 * Fit a LinearCombinationKernel to an Image, allowing the coefficients of the components to vary
 *
 * @return std::pair(best-fit kernel, std::pair(amp, chi^2))
 */
template <typename Image>
std::pair<std::shared_ptr<afw::math::Kernel>, std::pair<double, double>> fitKernelToImage(
        afw::math::LinearCombinationKernel const& kernel,  ///< the Kernel to fit
        Image const& image,                                ///< the image to be fit
        geom::Point2D const& pos                           ///< the position of the object
        ) {
    std::pair<std::vector<double>, afw::math::KernelList> const fit =
            fitKernelParamsToImage(kernel, image, pos);
    std::vector<double> params = fit.first;
    afw::math::KernelList kernels = fit.second;
    int const nKernel = params.size();
    assert(kernels.size() == static_cast<unsigned int>(nKernel));

    double amp = 0.0;
    for (int i = 0; i != nKernel; ++i) {
        std::shared_ptr<afw::math::Kernel> base = kernels[i];
        std::shared_ptr<afw::math::FixedKernel> k = std::static_pointer_cast<afw::math::FixedKernel>(base);
        amp += params[i] * k->getSum();
    }

    std::shared_ptr<afw::math::Kernel> outputKernel(new afw::math::LinearCombinationKernel(kernels, params));
    double chisq = 0.0;
    outputKernel->setCtr(kernels[0]->getCtr());

    return std::make_pair(outputKernel, std::make_pair(amp, chisq));
}

/************************************************************************************************************/
//
// Explicit instantiations
//
/// \cond
typedef float Pixel;

template std::pair<std::shared_ptr<afw::math::LinearCombinationKernel>, std::vector<double>>
createKernelFromPsfCandidates<Pixel>(afw::math::SpatialCellSet const&, geom::Extent2I const&,
                                     geom::Point2I const&, int const, int const, int const, int const,
                                     bool const, int const);
template int countPsfCandidates<Pixel>(afw::math::SpatialCellSet const&, int const);

template std::pair<bool, double> fitSpatialKernelFromPsfCandidates<Pixel>(afw::math::Kernel*,
                                                                          afw::math::SpatialCellSet const&,
                                                                          int const, double const,
                                                                          double const);
template std::pair<bool, double> fitSpatialKernelFromPsfCandidates<Pixel>(afw::math::Kernel*,
                                                                          afw::math::SpatialCellSet const&,
                                                                          bool const, int const, double const,
                                                                          double const);

template double subtractPsf(afw::detection::Psf const&, afw::image::MaskedImage<float>*, double, double,
                            double);

template std::pair<std::vector<double>, afw::math::KernelList> fitKernelParamsToImage(
        afw::math::LinearCombinationKernel const&, afw::image::MaskedImage<Pixel> const&,
        geom::Point2D const&);

template std::pair<std::shared_ptr<afw::math::Kernel>, std::pair<double, double>> fitKernelToImage(
        afw::math::LinearCombinationKernel const&, afw::image::MaskedImage<Pixel> const&,
        geom::Point2D const&);
/// \endcond

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
