
#include "lsst/meas/algorithms/StarFinder.h"
#include "lsst/meas/algorithms/Shapelet.h"
#include "lsst/meas/algorithms/ShapeletPsfCandidate.h"
#include "lsst/meas/algorithms/shapelet/StarFinderAlgo.h"

namespace lsst {
namespace meas {
namespace algorithms {

    // All of the functionality is imported from shapelet::StarFinder
    // Just repeat the constructor and destructor.
    class StarFinderImpl : public shapelet::StarFinderAlgo
    {
        typedef lsst::pex::policy::Policy Policy;
        typedef shapelet::StarFinderAlgo base;
        typedef shapelet::ConfigFile ConfigFile;

    public :
        StarFinderImpl(ConfigFile& params, const Policy& policy) :
            base(params,""),
            _cellSize(policy.getInt("cellSize")),
            _aperture(policy.getDouble("aperture"))
        {}

        ~StarFinderImpl() {}

        int getCellSize() const { return _cellSize; }
        double getAperture() const { return _aperture; }

    private :
        // This parameter is used by the surface layer, not StarFinderAlgo.
        int _cellSize;
        double _aperture;
    };

    StarFinder::StarFinder(const Policy& policy)
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
        params["maxratio1"] = 3.*double(params["purityratio"]);
        params["binsize1"] = 0.1;
        params["minbinsize"] = 0.01;
        params["miniter1"] = 3;
        params["miniter2"] = 2;
        params["startn2"] = 3*int(params["startn1"]);
        params["okvalcount"] = 2;
        params["maxrms"] = 0.05;
        params["maxrefititer"] = 5;

        pImpl = new StarFinderImpl(params,policy);
    }

    StarFinder::~StarFinder()
    { delete pImpl; pImpl = 0; }

    double StarFinder::calculateSourceSize(
        const Source& source, 
        const MaskedImage& image, const Wcs& wcs) const
    {
        double sigma = sqrt(source.getIxx() + source.getIyy());
        Shapelet shape(4,sigma);
        double x = getSourceX(source);
        double y = getSourceY(source);
        PointD pos = lsst::afw::geom::makePointD(x,y);
        if (shape.measureFromImage(
                source,pos,false,false,pImpl->getAperture(),image,wcs)) {
            return shape.getSigma();
        } else {
            return -1.;
        }
    }

    double StarFinder::calculateSourceMagnitude(const Source& source) const
    { return -2.5*log10(source.getPetroFlux()); }

    double StarFinder::getSourceX(const Source& source) const
    { return source.getXAstrom(); }
    double StarFinder::getSourceY(const Source& source) const
    { return source.getYAstrom(); }

    lsst::afw::math::SpatialCellSet::Ptr StarFinder::findStars(
        const SourceSet& allObj, 
        const MaskedImage& image, const Wcs& wcs) const
    {
        using lsst::afw::image::PointI;
        using lsst::afw::image::BBox;
        using shapelet::PotentialStar;
        using shapelet::Position;

        std::vector<PotentialStar*> maybestars;

        // First get a list of potential stars
        const int nObj = allObj.size();

        for (int i=0; i<nObj; ++i) {

            double x = getSourceX(*allObj[i]);
            double y = getSourceY(*allObj[i]);
            double size = calculateSourceSize(*allObj[i],image,wcs);
            double mag = calculateSourceMagnitude(*allObj[i]);
            Position pos(x,y);

            // Range checking
            if (!pImpl->isOkSize(size)) {
                continue;
            }
            if (!pImpl->isOkMag(mag)) {
                continue;
            }

            double logSize = pImpl->convertToLogSize(size);
            maybestars.push_back(
                new PotentialStar(pos,mag,logSize,i,""));
        }

        // Run the actual algorithm
        std::vector<PotentialStar*> stars = pImpl->findStars(maybestars);

        // Convert the results back into a SpatialCellSet
        int cellSize = pImpl->getCellSize();
        SpatialCellSet::Ptr ret(new SpatialCellSet(
            // This way of making the BBox comes from Lupton's example.
            // Why doesn't Image have a getBBox() method?
            BBox(PointI(0, 0), (image.getImage())->getWidth(), 
                 (image.getImage())->getHeight()),
            cellSize,cellSize));
        const int nStars = stars.size();
        for (int k=0; k<nStars;++k) {
            if (pImpl->isOkOutputMag(stars[k]->getMag())) {
                int i=stars[k]->getIndex();
                ShapeletPsfCandidate::Ptr cand(
                    new ShapeletPsfCandidate(
                        stars[k]->getPos().getX(),
                        stars[k]->getPos().getY(),
                        stars[k]->getSize(),
                        allObj[i]));
                ret->insertCandidate(cand);
            }
        }

        return ret;
    }

}}}
