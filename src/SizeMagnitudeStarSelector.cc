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
#include "lsst/meas/algorithms/SizeMagnitudeStarSelector.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/shapelet/SizeMagnitudeStarSelectorAlgo.h"
#include "lsst/afw/math/SpatialCell.h"

namespace lsst { namespace meas { namespace algorithms {

// All of the functionality is imported from shapelet::SizeMagnitudeStarSelector
// Just repeat the constructor and destructor.
class SizeMagnitudeStarSelectorImpl : public shapelet::SizeMagnitudeStarSelectorAlgo
{
    typedef lsst::pex::policy::Policy Policy;
    typedef shapelet::SizeMagnitudeStarSelectorAlgo base;
    typedef shapelet::ConfigFile ConfigFile;

public :
    SizeMagnitudeStarSelectorImpl(ConfigFile& params, const Policy& policy) :
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
    params["minsize"] = policy.getDouble("minSize");
    params["maxsize"] = policy.getDouble("maxSize");
    params["logsize"] = policy.getBool("isSizeLog");
    params["minmag"] = policy.getDouble("minMag");
    params["maxmag"] = policy.getDouble("maxMag");
    params["starfrac"] = policy.getDouble("starFrac");
    params["startn1"] = policy.getDouble("startN");
    params["fitorder"] = policy.getInt("fitOrder");
    params["fitsigclip"] = policy.getDouble("fitSigClip");
    params["starsperbin"] = policy.getInt("fitStars");
    params["purityratio"] = policy.getDouble("purity");

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
    const Source& source, 
    const Exposure& exposure) const
{
    double sigma = sqrt(source.getIxx() + source.getIyy());
    if (!(sigma > 0.)) return -1.;
    Shapelet shape(4, sigma);
    double x = getSourceX(source);
    double y = getSourceY(source);
    PointD pos(x, y);
    if (shape.measureFromImage(
            source, pos, false, false, pImpl->getAperture(), exposure)) {
        return shape.getSigma();
    } else {
        return -1.;
    }
}

double SizeMagnitudeStarSelector::calculateSourceMagnitude(const Source& source) const
{ return -2.5*log10(source.getPetroFlux()); }

double SizeMagnitudeStarSelector::getSourceX(const Source& source) const
{ return source.getXAstrom(); }
double SizeMagnitudeStarSelector::getSourceY(const Source& source) const
{ return source.getYAstrom(); }

SizeMagnitudeStarSelector::PsfCandidateList SizeMagnitudeStarSelector::selectStars(
    const Exposure& exposure,
    const SourceSet& sourceList) const
{
    const unsigned int MIN_OBJ_TO_TRY = 30;

    typedef Exposure::MaskedImageT MaskedImage;
    std::vector<shapelet::PotentialStar*> maybeStars;

    // First get a list of potential stars
    const int nSources = sourceList.size();
    for (int i=0; i<nSources; ++i) {
        //std::cout<<"Object "<<i<<"/"<<nSources<<std::endl;

        double x = getSourceX(*sourceList[i]);
        double y = getSourceY(*sourceList[i]);
        //std::cout<<"x,y = "<<x<<','<<y<<std::endl;
        double size = calculateSourceSize(*sourceList[i], exposure);
        //std::cout<<"size = "<<size<<std::endl;
        double mag = calculateSourceMagnitude(*sourceList[i]);
        //std::cout<<"mag = "<<mag<<std::endl;
        shapelet::Position pos(x, y);

        // Range checking
        if (!pImpl->isOkSize(size)) {
            continue;
        }
        if (!pImpl->isOkMag(mag)) {
            continue;
        }
        //std::cout<<"ok size, mag\n";

        double logSize = pImpl->convertToLogSize(size);
        maybeStars.push_back(
            new shapelet::PotentialStar(pos, mag, logSize, i, ""));
    }
    //std::cout<<"Total potential stars = "<<maybeStars.size()<<std::endl;
    if (maybeStars.size() < MIN_OBJ_TO_TRY) {
        // Too few objects for algorithm to have any chance of producing 
        // reasonable output.
        std::cerr<<"Only "<<maybeStars.size()<<" viable objects for star selection.\n";
        std::cerr<<"This algorithm needs at least "<<MIN_OBJ_TO_TRY<<" objects to try to find stars.\n";
        return PsfCandidateList();
    }

    // Run the actual algorithm
    std::vector<shapelet::PotentialStar*> stars = pImpl->findStars(maybeStars);
    //std::cout<<"Identified "<<stars.size()<<" stars\n";

    // Convert the results into a PsfCandidateList
    //MaskedImage::ConstPtr imagePtr = MaskedImage::ConstPtr(new MaskedImage(exposure.getMaskedImage(), false));
    Exposure::ConstPtr expPtr = Exposure::ConstPtr(new Exposure(exposure, false));
    PsfCandidateList psfCandidateList;
    const int nStars = stars.size();
    for (int k=0; k<nStars; ++k) {
        if (pImpl->isOkOutputMag(stars[k]->getMag())) {
            int i=stars[k]->getIndex();
            PsfCandidateT::Ptr psfCandidate(new PsfCandidateT(
                *(sourceList[i]),
                expPtr,
                stars[k]->getPos().getX(),
                stars[k]->getPos().getY()));
            psfCandidateList.push_back(psfCandidate);
        }
    }

    return psfCandidateList;
}

}}} // namespace lsst::meas::algorithms
