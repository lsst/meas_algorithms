#ifndef MeasAlgoShapeletStarFinderAlgo_H
#define MeasAlgoShapeletStarFinderAlgo_H

#include <vector>
#include <stdexcept>
#include "lsst/meas/algorithms/shapelet/ConfigFile.h"
#include "lsst/meas/algorithms/shapelet/Function2D.h"
#include "lsst/meas/algorithms/shapelet/PotentialStar.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class StarFinderException : 
        public std::runtime_error 
    {
    public:
        StarFinderException(const std::string& m) : std::runtime_error(m) {};
    };

    class StarFinderAlgo
    {
    public :
        StarFinderAlgo(const ConfigFile& params, std::string keyPrefix);

        std::vector<PotentialStar*> findStars(
            std::vector<PotentialStar*>& allObj);

        void findMinMax(
            const std::vector<PotentialStar*>& list,
            double *min, double *max, const Function2D& f);

        void rejectOutliers(
            std::vector<PotentialStar*>& list, 
            double nSigma, double minSigma, const Function2D& f);

        std::vector<PotentialStar*> getPeakList(
            const std::vector<PotentialStar*>& objList,
            double binSize, double minSize, double maxSize,
            int nStart, int minIter, double magStep, double maxSignifRatio,
            bool isFirstPass, const Function2D& f);

        void fitStellarSizes(
            Function2D *f, int order, double sigClip,
            const std::vector<PotentialStar*>& starList, double *outSigma);

        void roughlyFitBrightStars(
            const std::vector<PotentialStar*>& objList,
            Function2D *f,double *outSigma);

        void setParams(
            const ConfigFile& params, std::string keyPrefix,
            bool mustExist=false);

        bool isOkSize(const double size)
        { return size >= _minSize && size <= _maxSize; }

        bool isOkMag(const double mag)
        { return mag >= _minMag && mag <= _maxMag; }

        bool isOkOutputMag(const double mag)
        { return mag >= _minMag && mag <= _maxOutMag; }

        double convertToLogSize(const double size)
        { return _isSizeLog ? size : std::log(size); }

        double getMinSize() const { return _minSize; }
        double getMaxSize() const { return _maxSize; }
        double getMinMag() const { return _minMag; }
        double getMaxMag() const { return _maxMag; }

    private :

        double _minSize;       // The min and max size to consider
        double _maxSize;
        bool _isSizeLog;       // true if sizes are already log(size)
        double _minMag;        // The min and max magnitude to consider
        double _maxMag;
        double _maxOutMag;
        int _nDivX;         // Divide the image into nx by ny subsections
        int _nDivY;

        // Parameters when finding stars in each subdivision
        double _nStart1;       // # objects to start with on first pass 
        double _starFrac;      // Of these, how many are probably stars?
        double _magStep1;      // Step size in mag for each successive pass
        double _reject1;       // nsigma rejection.  (actually n x quartile size)
        double _maxRatio1;     // max ratio of count in valley / count in peak
        double _binSize1;      // binsize of histogram 
        int _okValCount;    // if valley count <= this, ok no matter what peak is
        int _minIter1;      // min number of mag steps to take
        double _maxRms;        // in initial linear fit, max rms to allow 

        // Parameters for last pass of whole image
        double _nStart2;       // # objects to start with
        double _magStep2;      // Step size in mag
        double _minBinSize;    // Min width of histogram bins
        double _reject2;       // n quartile rejection
        double _purityRatio;   // max ratio of count in valley / count in peak
        int _minIter2;      // min number of mag steps to take

        // Parameters for fitting results of subsections
        int _starsPerBin;   // How many stars per subsection?
        int _fitOrder;      // order of fit for size(x,y)
        double _fitSigClip;    // nsigma rejection in this fit
        int _maxRefitIter;  // max number of times to refit whole image

        bool _shouldOutputDesQa;  // Output warnings in DES QA format

    };

}}}}
#endif
