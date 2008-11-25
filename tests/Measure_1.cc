// -*- lsst-c++ -*-
#include <iostream>
#include <lsst/detection/Measure.h>
#include <lsst/pex/exceptions/Exception.h>
#include <lsst/pex/logging/Trace.h>
#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/image/Wcs.h>
#include <lsst/afw/detection/Source.h>
#include "lsst/daf/persistence/BoostStorage.h"
#include "lsst/daf/persistence/DbStorage.h"
#include "lsst/daf/persistence/DbTsvStorage.h"
#include "lsst/daf/persistence/Formatter.h"
#include "lsst/daf/persistence/LogicalLocation.h"
#include "lsst/daf/persistence/Persistence.h"

using namespace lsst::afw::detection;
using namespace lsst::afw::image;
using namespace lsst::detection;
using namespace lsst::daf::persistence;
namespace pexe = lsst::pex::exceptions;

typedef float ImagePixelT;
typedef unsigned int MaskPixelT;

int main(int argc, char**argv) {

    lsst::pex::logging::Trace::setDestination(std::cout);
    lsst::pex::logging::Trace::setVerbosity(".", 1);
    
    Exposure<ImagePixelT, maskPixelType> inpExposure;

    try {
	 inpExposure.readFits(argv[1]);
    } catch (pexe::ExceptionStack &e) {
	 std::cerr << "Error processing Exposure " << argv[1] << ": " << std::endl << e.what() << std::endl;
	 exit(1);
    }

    // Unpack the MaskedImage from the Exposure

    MaskedImage<ImagePixelT, maskPixelType> img = inpExposure.getMaskedImage();

    // Crudely estimate noise from mean of variance image - should do sigma clipping

    MaskedImage<ImagePixelT, maskPixelType>::ImagePtrT varianceImagePtr = img.getVariance();
    //
    float noise = sqrt(mean_channel_value(varianceImagePtr->getIVw()));

    float thresh = atof(argv[2]);

    std::cout << "Using threshold: " << thresh*noise << std::endl;

    const int nPixMin = 5;   

    bool polarity = false;

    DetectionSet<ImagePixelT, maskPixelType> ds(img, Threshold(thresh*noise, Threshold::VALUE, polarity), "FP", nPixMin);

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

    Measure<ImagePixelT, maskPixelType> mimg(img, "FP");

    std::vector<Footprint::PtrType>& fpVec = ds.getFootprints();
    lsst::afw::detection::SourceVector outputDiaSources;

    for (unsigned int i=0; i<fpVec.size(); i++) {
	 Source::Ptr diaPtr(new lsst::afw::detection::Source);
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

    lsst::daf::base::DataProperty::PtrType additionalData = lsst::daf::base::DataProperty::createPropertyNode("info");
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
