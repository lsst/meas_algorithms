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
 
#include <iostream>
#include "lsst/detection/Measure.h"
#include "lsst/pex/exceptions/Exception.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/daf/persistence/BoostStorage.h"
#include "lsst/daf/persistence/DbStorage.h"
#include "lsst/daf/persistence/DbTsvStorage.h"
#include "lsst/daf/persistence/Formatter.h"
#include "lsst/daf/persistence/LogicalLocation.h"
#include "lsst/daf/persistence/Persistence.h"

namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace lsstDetection = lsst::detection;
namespace dafPersist = lsst::daf::persistence;
namespace dafBase = lsst::daf::base;
namespace pexe = lsst::pex::exceptions;

typedef float ImagePixel;
typedef unsigned int MaskPixel;

int main(int argc, char **argv) {

    lsst::pex::logging::Trace::setDestination(std::cout);
    lsst::pex::logging::Trace::setVerbosity(".", 1);
    
    Exposure<ImagePixel, maskPixelType> inpExposure;

    try {
        inpExposure.readFits(argv[1]);
    } catch (pexe::ExceptionStack &e) {
        std::cerr << "Error processing Exposure " << argv[1] << ": " << std::endl << e.what() << std::endl;
        exit(1);
    }
    
    // Unpack the MaskedImage from the Exposure
    
    MaskedImage<ImagePixel, maskPixelType> img = inpExposure.getMaskedImage();

    // Crudely estimate noise from mean of variance image - should do sigma clipping

    MaskedImage<ImagePixel, maskPixelType>::ImagePtrT varianceImagePtr = img.getVariance();
    //
    float noise = sqrt(mean_channel_value(varianceImagePtr->getIVw()));

    float thresh = atof(argv[2]);

    std::cout << "Using threshold: " << thresh*noise << std::endl;

    int const nPixMin = 5;   

    bool polarity = false;

    DetectionSet<ImagePixel, maskPixelType> ds(img,
                                               Threshold(thresh*noise, Threshold::VALUE, polarity),
                                               "FP", nPixMin);

    //  Write out the image for debug so that the DetectionSet can be seen in the mask plane   
    
    std::string outName(argv[1]);

    try {
        outName = outName + "_OP";
        img.writeFits(outName);
    } catch (pexe::ExceptionStack &e) {
        std::cerr << "Failed to write " << outName << ": " << e.what() << std::endl;
        exit(1);
    }
    
    
    // Now measure the footprints.  Background is explicitly zero, as expected for difference image
    
    Wcs imgWcs = inpExposure.getWcs();
    
    Measure<ImagePixel, maskPixelType> mimg(img, "FP");
    
    std::vector<Footprint::PtrType>& fpVec = ds.getFootprints();
    afwDetection::SourceVector outputDiaSources;

    for (unsigned int i = 0; i < fpVec.size(); i++) {
        Source::Ptr diaPtr(new afwDetection::Source);
        diaPtr->setId(i); // will need to include Exposure id here!
        mimg.measureSource(diaPtr, *fpVec[i], 0);
        // use imgWcs to put ra and dec into DiaSource
        Coord2D pixCoord(diaPtr->getColc(), diaPtr->getRowc());
        Coord2D skyCoord = imgWcs.colRowToRaDec(pixCoord);
        diaPtr->setRa(skyCoord[0]);
        diaPtr->setDec(skyCoord[1]);
        std::cout << boost::format("DiaSource: %ld %lf %lf -> %lf %lf\n") % diaPtr->getId() 
            % diaPtr->getColc() % diaPtr->getRowc() % diaPtr->getRa() % diaPtr->getDec();
        outputDiaSources.push_back(*diaPtr);
    }
    
    // Now, try persisting the DiaSourceVector
    
    // Define a blank Policy.
    lsst::pex::policy::Policy::Ptr policy(new lsst::pex::policy::Policy);
    
    // Get a unique id for this test.
    struct timeval tv;
    gettimeofday(&tv, 0);      
    long long testId = tv.tv_sec * 1000000LL + tv.tv_usec;

    std::ostringstream os;
    os << testId;
    std::string testIdString = os.str();

    dafBase::DataProperty::PtrType additionalData = dafBase::DataProperty::createPropertyNode("info");
    lsst::daf::base::DataProperty::PtrType child1(new lsst::daf::base::DataProperty("visitId", testId));
    lsst::daf::base::DataProperty::PtrType child2(new lsst::daf::base::DataProperty("sliceId", 0));
    additionalData->addProperty(child1);
    additionalData->addProperty(child2);

    Persistence::Ptr persist = Persistence::getPersistence(policy);
    Storage::List storageList;
    LogicalLocation pathLoc("outputDiaSources.boost." + testIdString);
    storageList.push_back(persist->getPersistStorage("BoostStorage", pathLoc));
    persist->persist(outputDiaSources, storageList, additionalData);

}
