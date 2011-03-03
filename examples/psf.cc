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
 
#include "lsst/pex/policy.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/afw/math/Random.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/meas/algorithms/SpatialModelPsf.h"

namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace algorithms = lsst::meas::algorithms;

// A test case for SpatialModelPsf
int main() {
    int const width = 100;
    int const height = 301;
    afwImage::MaskedImage<float>::Ptr mi(new afwImage::MaskedImage<float>(width, height));
    *mi->getImage() = 0;
    float const sd = 3;                 // standard deviation of image
    *mi->getVariance() = sd*sd;
    mi->getMask()->addMaskPlane("DETECTED");
    
    double const FWHM = 5;
    int const ksize = 25;                         // size of desired kernel
    afwMath::Random rand;                         // make these tests repeatable by setting seed
    
    afwMath::randomGaussianImage(mi->getImage().get(), rand); // N(0, 1)
    *mi->getImage() *= sd;                                    // N(0, sd^2)
    
    std::pair<int, int> xy[] = {
        std::pair<int, int>(20, 20),
        std::pair<int, int>(60, 20),
        std::pair<int, int>(30, 35),
        std::pair<int, int>(50, 50),
        std::pair<int, int>(50, 130),
        std::pair<int, int>(70, 80),
        std::pair<int, int>(60, 210),
        std::pair<int, int>(20, 210)
    };
    
    for (int i = 0; i != sizeof(xy)/sizeof(xy[0]); ++i) {
        int x = xy[i].first, y = xy[i].second;
        afwDetection::Source source;
        
        double const flux = 10000 - 0*x - 10*y;
        
        double const sigma = 3 + 0.005*(y - mi->getHeight()/2);
        afwDetection::Psf::Ptr psf = afwDetection::createPsf("DoubleGaussian", ksize, ksize, sigma, 1, 0.1);
        afwImage::Image<float> im(*psf->computeImage(), true);
        im *= flux;
        afwImage::Image<float> smi(*mi->getImage(),
                                   afwImage::BBox(afwImage::PointI(x - ksize/2, y - ksize/2), ksize, ksize));
        
        float const dx = rand.uniform() - 0.5;
        float const dy = rand.uniform() - 0.5;
        {
            afwImage::Image<float>::Ptr oim = afwMath::offsetImage(im, dx, dy);
            
            smi += *oim;
        }
    }

#if 0
    mi->writeFits("foo.fits");
#endif

    afwDetection::Psf::Ptr psf = afwDetection::createPsf("DoubleGaussian", ksize, ksize,
                                                         FWHM/(2*sqrt(2*log(2))), 1, 0.1);

    afwMath::SpatialCellSet cellSet(afwImage::BBox(afwImage::PointI(0, 0), width, height), 100);
    afwDetection::FootprintSet<float> fs(*mi, afwDetection::Threshold(100), "DETECTED");
    afwDetection::FootprintSet<float>::FootprintList objects = *fs.getFootprints();
    //
    // Prepare to measure
    //
    lsst::pex::policy::Policy moPolicy;
    moPolicy.add("astrometry.SDSS.use", 1);
    moPolicy.add("source.astrom", "SDSS");
    moPolicy.add("photometry.NAIVE.radius", 3.0); // use NAIVE (== crude aperture)  photometry
    moPolicy.add("source.psfFlux", "NAIVE"); // Use the NAIVE flux in Source.getPsfFlux(); PSF would probably be better
    
    afwImage::Exposure<float>::Ptr exposure = afwImage::makeExposure(*mi);
    exposure->setPsf(psf);
    algorithms::MeasureSources<afwImage::Exposure<float> >::Ptr measureSources =
        algorithms::makeMeasureSources(exposure, moPolicy);
    
    afwDetection::SourceSet sourceList;
    for (unsigned int i = 0; i != objects.size(); ++i) {
        afwDetection::Source::Ptr source = afwDetection::Source::Ptr(new afwDetection::Source);
        sourceList.push_back(source);
        
        source->setId(i);
        source->setFlagForDetection(source->getFlagForDetection() | algorithms::Flags::BINNED1);
        source->setFootprint(objects[i]);

        measureSources->apply(source);

        algorithms::PsfCandidate<afwImage::MaskedImage<float> >::Ptr candidate = algorithms::makePsfCandidate(*source, mi);
        cellSet.insertCandidate(candidate);
    }

    // Convert our cellSet to a LinearCombinationKernel

    int const nEigenComponents = 2;
    int const spatialOrder  =    1;
    int const kernelSize =      31;
    int const nStarPerCell =     4;
    int const nIterForPsf =      5;

    algorithms::PsfCandidate<afwImage::MaskedImage<float> >::setWidth(kernelSize);
    algorithms::PsfCandidate<afwImage::MaskedImage<float> >::setHeight(kernelSize);

    for (int iter = 0; iter != nIterForPsf; ++iter) {
        algorithms::createKernelFromPsfCandidates<float>(cellSet, nEigenComponents, spatialOrder,
                                                         kernelSize, nStarPerCell);
    }
}
