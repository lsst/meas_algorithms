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
 
#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>

#include "lsst/base.h"
#include "lsst/utils/Utils.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Statistics.h"
#include "lsst/meas/algorithms/CR.h"

using namespace std;
namespace eups = lsst::utils::eups;
namespace pexLogging = lsst::pex::logging; 
namespace pexPolicy = lsst::pex::policy;
namespace afwDetection = lsst::afw::detection;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace measAlg = lsst::meas::algorithms;
namespace ex = lsst::pex::exceptions;

typedef afwImage::MaskedImage<float> MImageF;
typedef pexPolicy::Policy PPolicy;
typedef afwDetection::Footprint Footprint;

/****************************************************************************************************/

int main() {
    pexLogging::Trace::setVerbosity("meas.algorithms.CR", 2);

    afwImage::MaskedImage<float>::Ptr im;  // declare im outside scope of try block
    try {
        string afwdata = eups::productDir("afwdata");
        im = MImageF::Ptr(new MImageF(afwdata + "/CFHT/D4/cal-53535-i-797722_1"));
        im->getMask()->addMaskPlane("DETECTED");
    } catch(lsst::pex::exceptions::NotFoundException const& e) {
        cerr << e << endl;
        return 1;
    }

    double const fwhm = 5;              // pixels
    PTR(afwDetection::Psf) psf = afwDetection::createPsf("DoubleGaussian", 0, 0, fwhm/(2*sqrt(2*log(2))));

    PPolicy::Ptr policy;
    try {
        policy = PPolicy::Ptr(PPolicy::createPolicy(eups::productDir("meas_algorithms") +
                                                    "/policy/CrRejectDictionary.paf") );
    } catch(ex::NotFoundException const& e) {
        cerr << e << endl;
        return 1;
    }

    cout << boost::format("%dx%d\n") % im->getWidth() % im->getHeight();
    //
    // to work
    //
    double const background = afwMath::makeStatistics(*im->getImage(), afwMath::MEAN).getValue();
    std::vector<Footprint::Ptr> crs = measAlg::findCosmicRays(*im, *psf, background, *policy);
    
    cout << boost::format("Found %d CRs\n") % crs.size();

    return 0;
}
