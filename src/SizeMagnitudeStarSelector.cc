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
#include "lsst/pex/logging.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/shapelet/SizeMagnitudeStarSelectorAlgo.h"
#include "lsst/afw/math/SpatialCell.h"

namespace pexLogging = lsst::pex::logging;

namespace lsst { namespace meas { namespace algorithms {

// All of the functionality is imported from shapelet::SizeMagnitudeStarSelector
// Just repeat the constructor and destructor.
class SizeMagnitudeStarSelectorImpl : public shapelet::SizeMagnitudeStarSelectorAlgo
{
    typedef lsst::pex::policy::Policy Policy;
    typedef shapelet::SizeMagnitudeStarSelectorAlgo base;
    typedef shapelet::ConfigFile ConfigFile;

public :
    SizeMagnitudeStarSelectorImpl(ConfigFile & params, const Policy& policy) :
        base(params,""),
        _aperture(policy.getDouble("aperture"))
    {}

    ~SizeMagnitudeStarSelectorImpl() {}

    double getAperture() const { return _aperture; }

private :
    // This parameter is used by the surface layer, not SizeMagnitudeStarSelectorAlgo.
    double _aperture;
};

SizeMagnitudeStarSelector::SizeMagnitudeStarSelector(const Policy& policy)
{
    using shapelet::ConfigFile;

    // Convert Policy info into my ConfigFile format:
    ConfigFile params;
    params["minsize"] = policy.getDouble("minsize");
    params["maxsize"] = policy.getDouble("maxsize");
    params["logsize"] = policy.getBool("logsize");
    params["minmag"] = policy.getDouble("minmag");
    params["maxmag"] = policy.getDouble("maxmag");
    params["starfrac"] = policy.getDouble("starfrac");
    params["startn1"] = policy.getDouble("startn1");
    params["fitorder"] = policy.getInt("fitorder");
    params["fitsigclip"] = policy.getDouble("fitsigclip");
    params["starsperbin"] = policy.getInt("starsperbin");
    params["purityratio"] = policy.getDouble("purityratio");

    // The rest of these are just given defaults.  
    // In my experience, there is not much reason to make these
    // run-time modifiable.
    params["maxoutmag"] = params["maxmag"];
    params["ndivx"] = 1;
    params["ndivy"] = 1;
    params["magstep1"] = 0.25;
    params["magstep2"] = 0.1;
    params["reject1"] = params["fitsigclip"];
    params["reject2"] = params["fitsigclip"];
    params["maxratio1"] = 3.0 * double(params["purityratio"]);
    params["binsize1"] = 0.1;
    params["minbinsize"] = 0.01;
    params["miniter1"] = 3;
    params["miniter2"] = 2;
    params["startn2"] = 3 * int(params["startn1"]);
    params["okvalcount"] = 2;
    params["maxrms"] = 0.05;
    params["maxrefititer"] = 5;

    pImpl = boost::shared_ptr<SizeMagnitudeStarSelectorImpl>(new SizeMagnitudeStarSelectorImpl(params, policy));
}

double SizeMagnitudeStarSelector::calculateSourceSize(
    const SourceRecord & source, 
    const Exposure & exposure) const
{
    double sigma = sqrt(source.getIxx() + source.getIyy());
    if (!(sigma > 0.)) return -1.;
    Shapelet shape(4, sigma);
    double x = source.getX();
    double y = source.getY();
    PointD pos(x, y);
    if (shape.measureFromImage(
            source, pos, false, false, pImpl->getAperture(), exposure)) {
        return shape.getSigma();
    } else {
        return -1.;
    }
}

/*
 * Calculates a magnitude for a source.
 *
 * The star finder is written in terms of using magnitudes rather than
 * fluxes, whereas Source stores fluxes; So this just translates the flux into a magnitude.
 *
 * Note This function may also be a good candidate for a having
 * its action be specifiable by a Policy parameter.
 */
static double calculateSourceMagnitude(lsst::afw::table::SourceRecord const & source,
                                       SizeMagnitudeStarSelector::Exposure const& exposure
                                      )
{
    return exposure.getCalib()->getMagnitude(source.getApFlux()); 
}

SizeMagnitudeStarSelector::PsfCandidateList SizeMagnitudeStarSelector::selectStars(
    const Exposure& exposure,
    const SourceCatalog & sourceList) const
{
    pexLogging::Debug traceLog("meas.algorithms.SizeMagnitudeStarSelector"); // trace output goes here
    const unsigned int MIN_OBJ_TO_TRY = 30;

    typedef Exposure::MaskedImageT MaskedImage;
    std::vector<shapelet::PotentialStar*> maybeStars;

    // First get a list of potential stars
    const int nSources = sourceList.size();
    traceLog.debug<4>("%d candidate stars", nSources);
    for (int i=0; i<nSources; ++i) {
        double const x = sourceList[i].getX();
        double const y = sourceList[i].getY();
        double const size = calculateSourceSize(sourceList[i], exposure);
        double const mag = calculateSourceMagnitude(sourceList[i], exposure);

        shapelet::Position pos(x, y);

        // Range checking
        bool ok = true;
        if (!pImpl->isOkSize(size)) {
            ok = false;
        }
        if (ok && !pImpl->isOkMag(mag)) {
            ok = false;
        }
        traceLog.debug<5>("i, x, y, size, mag = %d %.1f, %.1f %g %g: %d", i, x, y, size, mag, ok);

        if (ok) {
            double logSize = pImpl->convertToLogSize(size);
            maybeStars.push_back(new shapelet::PotentialStar(pos, mag, logSize, i, ""));
        }
    }
    traceLog.debug<4>("Total potential stars = %d", maybeStars.size());
    if (maybeStars.size() < MIN_OBJ_TO_TRY) {
        // Too few objects for algorithm to have any chance of producing reasonable output.
         pex::logging::Log::getDefaultLog().log(pexLogging::Log::WARN,
                  str(boost::format("Only %d viable objects for star selection. "
                                    "This algorithm needs at least %d objects to try to find stars")
                      % maybeStars.size() % MIN_OBJ_TO_TRY));
                                                                
        return PsfCandidateList();
    }

    // Run the actual algorithm
    std::vector<shapelet::PotentialStar*> stars = pImpl->findStars(maybeStars);
    traceLog.debug<4>("Identified %d stars", stars.size());

    // Convert the results into a PsfCandidateList
    //MaskedImage::ConstPtr imagePtr = MaskedImage::ConstPtr(new MaskedImage(exposure.getMaskedImage(), false));
    CONST_PTR(Exposure) expPtr = CONST_PTR(Exposure)(new Exposure(exposure, false));
    PsfCandidateList psfCandidateList;
    const int nStars = stars.size();
    for (int k=0; k<nStars; ++k) {
        if (pImpl->isOkOutputMag(stars[k]->getMag())) {
            int i=stars[k]->getIndex();
            PTR(PsfCandidateT) psfCandidate(new PsfCandidateT(
                sourceList.get(i),
                expPtr,
                stars[k]->getPos().getX(),
                stars[k]->getPos().getY()));
            psfCandidateList.push_back(psfCandidate);
        }
    }

    return psfCandidateList;
}

}}} // namespace lsst::meas::algorithms
