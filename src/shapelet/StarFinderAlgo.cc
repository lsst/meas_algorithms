
#include <algorithm>
#include <fstream>
#include <iostream>
#include <functional>
#include <vector>
#include <sstream>
#include <string>
#include <string>
#include <valarray>
#include "unistd.h"

#include "lsst/meas/algorithms/shapelet/StarFinderAlgo.h"
#include "lsst/meas/algorithms/shapelet/Legendre2D.h"
#include "lsst/meas/algorithms/shapelet/Bounds.h"
#include "lsst/meas/algorithms/shapelet/PotentialStar.h"
#include "lsst/meas/algorithms/shapelet/Histogram.h"
#include "lsst/meas/algorithms/shapelet/ConfigFile.h"
#include "lsst/meas/algorithms/shapelet/fspd.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

#define SFKeyAssign(var,key) \
    do { \
        if (mustexist || params.keyExists(keyPrefix + key)) { \
            var = params[ keyPrefix + key ];\
        } \
        xdbg << "SFAssign " key " = " << var; \
        xdbg << " from parameter "<<(keyPrefix + key)<<std::endl; \
    } while (false) 

    void StarFinderAlgo::setParams(
        const ConfigFile& params, std::string keyPrefix, bool mustexist)
    {
        SFKeyAssign(_minSize,"minsize");
        SFKeyAssign(_maxSize,"maxsize");
        SFKeyAssign(_isSizeLog,"logsize");

        SFKeyAssign(_minMag,"minmag");
        SFKeyAssign(_maxMag,"maxmag");
        SFKeyAssign(_maxOutMag,"maxoutmag");
        SFKeyAssign(_nDivX,"ndivx");
        SFKeyAssign(_nDivY,"ndivy");

        SFKeyAssign(_nStart1,"startn1");
        SFKeyAssign(_starFrac,"starfrac");
        SFKeyAssign(_magStep1,"magstep1");
        SFKeyAssign(_reject1,"reject1");
        SFKeyAssign(_maxRatio1,"maxratio1");
        SFKeyAssign(_binSize1,"binsize1");
        SFKeyAssign(_okValCount,"okvalcount");
        SFKeyAssign(_minIter1,"miniter1");
        SFKeyAssign(_maxRms,"maxrms");

        SFKeyAssign(_nStart2,"startn2");
        SFKeyAssign(_magStep2,"magstep2");
        SFKeyAssign(_minBinSize,"minbinsize");
        SFKeyAssign(_reject2,"reject2");
        SFKeyAssign(_purityRatio,"purityratio");
        SFKeyAssign(_minIter2,"miniter2");

        SFKeyAssign(_starsPerBin,"starsperbin");
        SFKeyAssign(_fitOrder,"fitorder");
        SFKeyAssign(_fitSigClip,"fitsigclip");
        SFKeyAssign(_maxRefitIter,"maxrefititer");
    }

    StarFinderAlgo::StarFinderAlgo(
        const ConfigFile& params, std::string keyPrefix)
    {
        std::string fspd((const char*)findstars_params_default,
                         findstars_params_default_len);
        std::istringstream is(fspd);
        ConfigFile defaultParams;
        defaultParams.setDelimiter("\t");
        defaultParams.setComment("#");
        defaultParams.read(is);
        setParams(defaultParams,"",true);

        setParams(params,keyPrefix);

        _shouldOutputDesQa = params.read("des_qa",false); 
    }

#undef SFKeyAssign

    std::vector<PotentialStar*> StarFinderAlgo::findStars(
        std::vector<PotentialStar*>& allObj)
    {
        dbg<<"starting FindStars, allobj.size() = "<<allObj.size()<<std::endl;

        // sort allObj by magnitude
        std::sort(allObj.begin(),allObj.end(),std::mem_fun(
                &PotentialStar::isBrighterThan));
        dbg<<"sorted allobj\n";

        // First stage:
        // Split area into sections
        // For each secion find as many stars as possible.
        // Use these stars to fit the variation in rsq across the image.

        // Make bounds of whole region called totalBounds
        Bounds totalBounds;
        const int nObj = allObj.size();
        for(int k=0;k<nObj;++k) totalBounds += allObj[k]->getPos();
        dbg<<"totalbounds = \n"<<totalBounds<<std::endl;

        // boundsar is the 3x3 (nDivX x ndivy) array of quadrant bounds
        // call it qBounds even though it won't be quadrants unless 2x2
        std::vector<Bounds> qBounds = totalBounds.divide(_nDivX,_nDivY);
        xdbg<<"made qbounds\n";

        // probStars will be our first pass list of probable stars
        std::vector<PotentialStar*> probStars;

        // For each section, find the stars and add to probStars
        const int nSection = qBounds.size();
        for(int i=0;i<nSection;++i) {
            dbg<<"i = "<<i<<": bounds = "<<qBounds[i]<<std::endl;

            // someObj are the objects in this section
            // Note that someObj is automatically sorted by magnitude, since
            // allObj was sorted.
            std::vector<PotentialStar*> someObj;
            for(int k=0;k<nObj;++k) {
                if (qBounds[i].includes(allObj[k]->getPos())) 
                    someObj.push_back(allObj[k]);
            }
            dbg<<"added "<<someObj.size()<<" obj\n";

            // Does a really quick and dirty fit to the bright stars
            // Basically it takes the 10 smallest of the 50 brightest objects,
            // finds the peakiest 5, then fits their sizes to a 1st order function.
            // It also gives us a rough value for the sigma
            Legendre2D fLinear(qBounds[i]);
            double sigma;
            roughlyFitBrightStars(someObj,&fLinear,&sigma);
            dbg<<"fit bright stars: sigma = "<<sigma<<std::endl;

            // Calculate the min and max values of the (adjusted) sizes
            double minSize,maxSize;
            findMinMax(someObj,&minSize,&maxSize,fLinear);
            dbg<<"min,max = "<<minSize<<','<<maxSize<<std::endl;

            // Find the objects clustered around the stellar peak.
            std::vector<PotentialStar*> qPeakList =
                getPeakList(
                    someObj,_binSize1,minSize,maxSize,
                    int(_nStart1*someObj.size()),_minIter1,_magStep1,_maxRatio1,
                    true,fLinear);
            const int nPeak = qPeakList.size();
            dbg<<"peaklist has "<<nPeak<<" stars\n";

            // Remove outliers using a median,percentile rejection scheme.
            // _binSize1/2. is the minimum value of "sigma".
            rejectOutliers(qPeakList,_reject1,_binSize1/2.,fLinear);
            dbg<<"rejected outliers, now have "<<nPeak<<" stars\n";

            // Use at most 10 (starsPerBin) stars per region to prevent one region 
            // from dominating the fit.  Use the 10 brightest stars to prevent being
            // position biased (as one would if it were based on size
            int nStarsExpected = int(_starFrac * someObj.size());
            if (nPeak < nStarsExpected) {
                if (nPeak < int(0.2 * nStarsExpected)) {
                    if (_shouldOutputDesQa) {
                        std::cout<<"STATUS3BEG Warning: Only "<<
                            qPeakList.size()<<" stars found in section "<<
                            i<<". STATUS3END"<<std::endl;
                    }
                    dbg<<"Warning: only "<<qPeakList.size()<<
                        " stars found in section "<<i<<
                        "  "<<qBounds[i]<<std::endl;
                }
                probStars.insert(probStars.end(),qPeakList.begin(),qPeakList.end());
            } else {
                std::sort(qPeakList.begin(),qPeakList.end(),
                          std::mem_fun(&PotentialStar::isBrighterThan));
                probStars.insert(probStars.end(),qPeakList.begin(),
                                 qPeakList.begin()+nStarsExpected);
            }
            dbg<<"added to probstars\n";
        }
        xdbg<<"done qbounds loop\n";
        int nStars = probStars.size();
        dbg<<"nstars = "<<nStars<<std::endl;

        // Now we have a first estimate of which objects are stars.
        // Fit a quadratic function to them to characterize the variation in size
        Legendre2D f(totalBounds);
        double sigma;
        fitStellarSizes(&f,_fitOrder,_fitSigClip,probStars,&sigma);

        // Second stage:
        // Use the fitted function for the size we just got to adjust
        // the measured sizes according to their position in the image.
        // Make the histogram using measured - predicted sizes (still log
        // size actually) which should then have the peak very close to 0.
        // Set the bin size for this new histogram to be the rms scatter of the
        // stellar peak from pass 1.
        // This time step by 0.25 mag (magStep2).

        // Find the values of minSize,maxSize for the whole thing with the new f
        double minSize,maxSize;
        findMinMax(allObj,&minSize,&maxSize,f);
        dbg<<"new minmax = "<<minSize<<','<<maxSize<<std::endl;

        // Find the objects clustered around the stellar peak.
        // Note that the binSize is a set fraction of sigma (0.5), so this
        // should make the stellar peak clearer.
        // Also, the f at the end means the functional fitted size will be
        // subtracted off before adding to the histogram.
        probStars = getPeakList(
            allObj,0.5*sigma,minSize,maxSize,
            int(_nStart2*allObj.size()),_minIter2,_magStep2,_purityRatio,false,f);
        nStars = probStars.size();
        dbg<<"probstars has "<<nStars<<" stars\n";

        // Remove outliers using a median,percentile rejection scheme.
        rejectOutliers(probStars,_reject2,sigma,f);
        dbg<<"rejected outliers, now have "<<probStars.size()<<" stars\n";
        // Worth doing twice, since first pass usually throws out a lot of junk.
        rejectOutliers(probStars,_reject2,sigma,f);  
        dbg<<"rejected outliers, now have "<<probStars.size()<<" stars\n";

        // If you get a bad fit the first time through, it can take a 
        // few passes to fix it.
        // And always do at least one refit.
        bool shouldRefit = true;
        for(int iter=0;shouldRefit && iter<_maxRefitIter;++iter) {
            dbg<<"starting refit\n"<<std::endl;
            shouldRefit = false;
            std::vector<std::vector<PotentialStar*> > starsList(nSection);

            // Add each star to the appropriate sublist
            for(int k=0;k<nStars;++k) {
                for(int i=0;i<nSection;++i) {
                    if(qBounds[i].includes(probStars[k]->getPos()))
                        starsList[i].push_back(probStars[k]);
                }
            }

            // Make fitList, the new list of the 10 brightest stars per section
            std::vector<PotentialStar*> fitList;
            for(int i=0;i<nSection;++i) {
                // if there are still < 10 stars, give a warning, and
                // just add all the stars to fitList
                if (int(starsList[i].size()) < _starsPerBin) {
                    if (_shouldOutputDesQa) {
                        std::cout<<"STATUS3BEG Warning: Only "<<
                            starsList[i].size()<<" stars in section "<<i<<
                            ". STATUS3END"<<std::endl;
                    }
                    dbg<<"Warning: only "<<starsList[i].size()<<
                        " stars in section ";
                    dbg<<i<<"  "<<qBounds[i]<<std::endl;
                    fitList.insert(fitList.end(),starsList[i].begin(),
                                   starsList[i].end());
                    shouldRefit = true;
                } else {
                    // sort the sublist by magnitude
                    std::sort(starsList[i].begin(),starsList[i].end(),
                              std::mem_fun(&PotentialStar::isBrighterThan));

                    // add the brightest 10 to fitList
                    fitList.insert(fitList.end(),starsList[i].begin(),
                                   starsList[i].begin()+_starsPerBin);
                }
            }
            // Do all the same stuff as before (refit f, get new min,max,
            // find the peakList again, and reject the outliers)
            nStars = fitList.size();
            fitStellarSizes(&f,_fitOrder,_fitSigClip,fitList,&sigma);
            findMinMax(allObj,&minSize,&maxSize,f);
            dbg<<"new minmax = "<<minSize<<','<<maxSize<<std::endl;
            probStars = getPeakList(
                allObj,0.5*sigma,minSize,maxSize,
                int(_nStart2*allObj.size()),_minIter2,_magStep2,
                _purityRatio,false,f);
            dbg<<"probstars has "<<probStars.size()<<" stars\n";
            rejectOutliers(probStars,_reject2,sigma,f);
            rejectOutliers(probStars,_reject2,sigma,f);
            dbg<<"rejected outliers, now have "<<probStars.size()<<" stars\n";
            if (int(probStars.size()) < nStars) {
                shouldRefit = true;
                dbg<<"fewer than "<<nStars<<" - so refit\n";
            }
            nStars = probStars.size();
        }

        dbg<<"done FindStars\n";

        for(int i=0;i<nStars;++i) {
            xdbg<<"stars["<<i<<"]: size = "<<probStars[i]->getSize();
            xdbg<<", size - f(pos) = ";
            xdbg<<probStars[i]->getSize()-f(probStars[i]->getPos())<<std::endl;
            xdbg<<probStars[i]->getLine()<<std::endl;
        }
        return probStars;
    }

    void StarFinderAlgo::findMinMax(
        const std::vector<PotentialStar*>& list, 
        double *min, double *max, const Function2D& f)
    {
        double min1,max1;
        min1 = max1 = list[0]->getSize()-f(list[0]->getPos());
        const int nStars = list.size();
        for(int k=1;k<nStars;++k) {
            double size = list[k]->getSize() - f(list[k]->getPos());
            if (size > max1) max1 = size;
            if (size < min1) min1 = size;
        }
        *min = min1;
        *max = max1;
    }

    void StarFinderAlgo::rejectOutliers(
        std::vector<PotentialStar*>& list, 
        double nSigma, double minSigma, const Function2D& f)
    {
        // This rejects outliers by finding the median and the quartile 
        // values, and calls sigma the average of the two quartile deviations.
        // "nSigma" is then how many of these "sigma" away from the median
        //  to consider something an outlier.

        const int nStars = list.size();
        if (nStars <= 4) return;
        // Find median, 1st and 3rd quartile stars;
        std::vector<PotentialStar*> modifList(nStars);
        for(int k=0;k<nStars;++k) {
            modifList[k] = new PotentialStar(*list[k]);
            double newSize = list[k]->getSize() - f(list[k]->getPos());
            modifList[k]->setSize(newSize);
        }
        std::sort(modifList.begin(),modifList.end(),
                  std::mem_fun(&PotentialStar::isSmallerThan));
        PotentialStar* mStar = modifList[modifList.size()/2];
        PotentialStar* q1Star = modifList[modifList.size()/4];
        PotentialStar* q3Star = modifList[modifList.size()*3/4];
        double median = mStar->getSize();
        double q1 = q1Star->getSize();
        double q3 = q3Star->getSize();
        double sigma = std::max(q3-median,median-q1);
        sigma = std::max(sigma,minSigma);
        xdbg<<"q1,m,q3 = "<<q1<<" , "<<median<<" , "<<q3<<std::endl;
        xdbg<<"sigma = "<<sigma<<", nsig * sigma = "<<nSigma*sigma<<std::endl;
        // Remove elements which std::abs(x->getSize()-median) > nSigma*sigma
        int j=0;
        for(int k=0; k<nStars;++k) {
            if (std::abs(modifList[k]->getSize()-median)>nSigma*sigma) {
                xdbg<<"k = "<<k<<", size = "<<modifList[k]->getSize();
                xdbg<<" might be an outlier (median = "<<median<<")\n";
            }
            list[j] = list[k];  // Note: update original list, not modified version.
            ++j;
        }
        if (j<nStars) list.erase(list.begin()+j,list.end());
        for(int k=0;k<nStars;++k) delete modifList[k];
    }

    std::vector<PotentialStar*> StarFinderAlgo::getPeakList(
        const std::vector<PotentialStar*>& objList,
        double binSize, double minSize, double maxSize,
        int nStart, int minIter, double magStep, double maxSignifRatio,
        bool isFirstPass, const Function2D& f)
    {
        if (binSize < _minBinSize) binSize = _minBinSize;

        // Make a histogram to find the stellar peak
        Histogram<PotentialStar*> hist(binSize,minSize,maxSize);

        // Start with the nStart brightest objects, then step up in
        // magnitude by magStep at a time until histogram peak
        // gets dilluted.
        int k;
        Assert(nStart<=int(objList.size()));
        for(k=0;k<nStart;++k) {
            hist.add(objList[k]->getSize()-f(objList[k]->getPos()),objList[k]);
        }
        xdbg<<"added "<<nStart<<" objects\n";
        if(XDEBUG && dbgout) hist.print(*dbgout,minSize,maxSize);

        // Get first estimates of the stellar and (first) galaxy peaks
        // along with the corresponding valleys.
        double peak1 = hist.findFirstPeakAfter(minSize);
        xdbg<<"peak1 = "<<peak1<<" ("<<hist[peak1]<<")\n";
        double valley1 = hist.findFirstValleyAfter(peak1);
        xdbg<<"valley1 = "<<valley1<<" ("<<hist[valley1]<<")\n";
        while (valley1 < 0.0) {
            // Since we subtract the fit, stars should be close to 0
            peak1 = hist.findFirstPeakAfter(valley1);
            xdbg<<"new peak1 = "<<peak1<<" ("<<hist[peak1]<<")\n";

            // Old code, did not do well with stupid no-galaxy images created
            // by DES pipeline.
            //valley1 = hist.findFirstValleyAfter(peak1);
            //xdbg<<"new valley1 = "<<valley1<<" ("<<hist[valley1]<<")\n";

            double newValley1 = hist.findFirstValleyAfter(peak1);
            if (!(newValley1 > valley1)) {
                std::ostringstream err;
                err<<"Couldn't find a stellar peak.  ";
                err<<"All the objects seem to be basically the same size.";
                throw StarFinderException(err.str());
            }
            valley1 = newValley1;
            xdbg<<"new valley1 = "<<valley1<<" ("<<hist[valley1]<<")\n";


        }
        double peak2 = hist.findFirstPeakAfter(valley1);
        xdbg<<"peak2 = "<<peak2<<" ("<<hist[peak2]<<")\n";
        // Sometimes the stellar peak is horned (eg. 212), in which case, find
        // the next peak
        double nBinsForHorn = isFirstPass ? 2.5 : 4.5;
        if (peak2-peak1 < nBinsForHorn*binSize && peak2+binSize < maxSize) {
            valley1 = hist.findFirstValleyAfter(peak2);
            peak2 = hist.findFirstPeakAfter(peak2+binSize);
            xdbg<<"horn pattern: next peak2 = "<<peak2<<" ("<<hist[peak2]<<")\n";
        }
        double valley2 = hist.findFirstValleyAfter(peak2);
        xdbg<<"valley2 = "<<valley2<<" ("<<hist[valley2]<<")\n";
        double valley0 = hist.findFirstValleyBefore(peak1);
        if (!isFirstPass && valley0 > 0.) valley0 = hist.findFirstValleyBefore(0.);
        xdbg<<"valley0 = "<<valley0<<" ("<<hist[valley0]<<")\n";

        // Get the stars in the first peak for the initial value of peakList
        // When we're all done it will be the list of objects at the stellar peak.
        std::vector<PotentialStar*> peakList=hist.getRefsInRange(valley0,valley1);

        // Define the highest mag so far in list
        double highMag = objList[nStart-1]->getMag();

        // Loop at least minIter times.  Done is set to false whenever
        // we've just found a new best peak.
        bool done=false;
        double prevValley=valley1;
        double prevPeak=0.;
        for(int count = 0; count < minIter || !done; ++count) {

            // Increment highMag
            highMag += magStep;
            xdbg<<"count = "<<count<<", highmag = "<<highMag<<std::endl;

            // Add new objects to histogram.
            const int nObj = objList.size();
            for(;k<nObj && objList[k]->getMag()<highMag;++k) {
                hist.add(objList[k]->getSize()-f(objList[k]->getPos()),objList[k]);
            }
            xdbg<<"hist has "<<k<<" objects\n";
            if(XDEBUG && dbgout) hist.print(*dbgout,valley0,valley2);

            xdbg<<"valley0 = "<<valley0<<std::endl;
            // Find the peak on the stellar side
            double peak = hist.findFirstPeakAfter(valley0,!isFirstPass);
            valley0 = hist.findFirstValleyBefore(peak); // for next pass
            xdbg<<"done findpeak: "<<peak<<std::endl;

            // Find the valley 
            // The !isFirstPass allows poisson noise fluctuations for the final pass
            double valley = hist.findFirstValleyAfter(peak,!isFirstPass);
            xdbg<<"done findvalley: "<<valley<<std::endl;

            // If valley ends up negative, find next peak and valley
            if (valley < 0.0) {
                double nextPeak = hist.findFirstPeakAfter(valley,!isFirstPass);
                // if next peak is closer to 0 than current one, use new values.
                while (std::abs(nextPeak) < std::abs(peak)) {
                    peak = nextPeak;
                    xdbg<<"negative valley - new peak: "<<peak<<std::endl;
                    valley = hist.findFirstValleyAfter(peak,!isFirstPass);
                    xdbg<<"new valley: "<<valley<<std::endl;
                    if (valley < 0.0) {
                        nextPeak = hist.findFirstPeakAfter(valley,!isFirstPass);
                    } else {
                        break;
                    }
                }
            }

            // For the first pass, be extra conservative and don't let the valley
            // position get much larger than a previous good value.
            // If the peak is also larger than prevValley, something weird is
            // going on, so bail on this pass.
            if (isFirstPass && valley > prevValley+binSize) {
                valley = prevValley;
                xdbg<<"moved valley to "<<valley<<std::endl;
                if (peak >= valley) { 
                    done = true;
                    continue; 
                }
            }

            // If new peak is larger than the previous valley, and the new peak 
            // isn't closer to 0, then reject it (probably added enough objects to
            // the valley to fill it in completely).
            if (peak > prevValley && std::abs(peak) > std::abs(prevPeak)) {
                xdbg<<"New peak passed previous valley.\n";
                done = true;
                continue;
            }

            // The refined counts allow the bins to be centered anywhere near
            // the respective peak and valley and finds how many objects
            // would fall in the optimally centered bin.
            int peakCount = hist.getRefinedPeakCount(&peak);
            xdbg<<"peakcount = "<<peakCount<<" centered at "<<peak<<std::endl;
            xdbg<<"before refined: valley = "<<valley<<std::endl;
            int valleyCount = hist.getRefinedValleyCount(&valley);
            xdbg<<"valleycount = "<<valleyCount<<" centered at "<<valley<<std::endl;

            // Check for the horned pattern and use the next valley if it is there.
            // But only if the new valley is at least as good as the first one.
            //double nBinsForHorn = isFirstPass ? 1.5 : 3.5;
            nBinsForHorn = 1.5;
            if ((valley-peak < nBinsForHorn*binSize) && 
                (hist.findFirstPeakAfter(valley) - valley) < nBinsForHorn*binSize) {
                double newValley = hist.findFirstValleyAfter(
                    hist.findFirstPeakAfter(valley));
                int newValleyCount = hist.getRefinedValleyCount(&newValley);
                if (((isFirstPass && newValleyCount <= valleyCount) ||
                     (!isFirstPass && valleyCount > 0 && newValleyCount == 0)) && 
                    newValley <= prevValley + nBinsForHorn*binSize) {

                    xdbg<<"horn pattern detected.  new valley = "<<
                        newValley<<std::endl;
                    xdbg<<"new valley count = "<<newValleyCount<<std::endl;
                    valley = newValley;
                    valleyCount = newValleyCount;
                }
            }

            // The significance of the peak is the ratio of the number in the
            // valley to the number in the peak
            // So _small_ significances are good.
            // Has to be this way, since valleyCount can be 0,
            // so peakCount/valleyCount can be undefined.
            double signif = (double)valleyCount/(double)peakCount;
            dbg<<"signif = "<<valleyCount<<"/"<<peakCount<<" = "<<signif<<std::endl;

            // If the valley count is <= okValCount then drop the significance
            // to 0 regardless of the peak count, since sometimes the peak is fairly
            // broad so peakCount isn't high enough to get a good ratio.
            if (isFirstPass && valleyCount <= _okValCount) {
                signif = 0.;
                dbg<<"reset signif to 0, since valleycount "<<valleyCount<<" <= ";
                dbg<<_okValCount<<std::endl;
            }

            // If we have a new good peak, get all the stars in it for peakList
            // Also make sure we repeat the loop by setting done = false
            if (signif <= maxSignifRatio) {
                dbg<<"signif = "<<signif<<" <= "<<maxSignifRatio<<std::endl;
                xdbg<<"good signif value\n";
                // Range is symmetrical about peakvalue with upper limit at valley
                double peakStart = std::min(2.*peak-valley,peak-binSize*0.5);
                std::vector<PotentialStar*> temp = 
                    hist.getRefsInRange(peakStart,valley);
                xdbg<<"get refs from "<<peakStart<<" to "<<valley<<std::endl;
                xdbg<<"got refs: "<<temp.size()<<std::endl;
                // If this list has fewer objects than a previous list,
                // and it doesn't just have a tighter range,
                // don't take it.  Surprisingly, this is actually possible, and it
                // means something weird happened
                if (temp.size() >= peakList.size() || 
                    (valley<prevValley && valley>prevPeak) || 
                    (peak > prevPeak) ) {

                    if (temp.size() >= peakList.size()) {
                        xdbg<<"tempsize = "<<temp.size()<<" >= "<<peakList.size();
                    } else if ((valley<prevValley && valley>prevPeak)) {
                        xdbg<<"valley = "<<valley<<" < "<<prevValley;
                    } else {
                        xdbg<<"peak = "<<peak<<" > "<<prevPeak;
                    }
                    xdbg<<" so keep adding...\n";
                    peakList = temp;
                    done = false;
                    prevValley = valley;
                    prevPeak = peak;
                } else {
                    xdbg<<"tempsize < "<<peakList.size();
                    xdbg<<" and peak,valley didn't contract, so stop adding now.\n";
                    done=true;
                }
            } else {
                dbg<<"signif = "<<signif<<" > "<<maxSignifRatio<<std::endl;
                done=true; // maybe.  still loop if count < minIter
            }

            // If we've already used all the stars, don't loop anymore
            if (k == nObj) {
                xdbg<<"no more to add\n";
                count = minIter; done = true; 
            }
        } // for count
        xdbg<<"done finding peaklist\n";
        return peakList;
    }

    void StarFinderAlgo::fitStellarSizes(
        Function2D *f, int order, double sigClip,
        const std::vector<PotentialStar*>& starList, double *outSigma)
    {
        // Make list of positions
        std::vector<Position> posList(starList.size());
        std::transform(starList.begin(),starList.end(),posList.begin(),
                       std::mem_fun(&PotentialStar::getPos));

        // Make list of sizes
        std::vector<double> sizelist(starList.size());
        std::transform(starList.begin(),starList.end(),sizelist.begin(),
                       std::mem_fun(&PotentialStar::getSize));

        // Use all stars in list (to start with anyway), so all true
        std::vector<bool> useList(starList.size(),true);

        if (int(starList.size()) <= (order+1)*(order+2)/2) {
            std::ostringstream err;
            err<<"Not enough stars to do fit.  ";
            err<<"Increase starsPerBin or decrease fitOrder.";
            throw StarFinderException(err.str());
        }

        // Do the fit.
        double chisq;
        int dof;
        xdbg<<"before outlier fit\n";
        f->outlierFit(order,sigClip,posList,sizelist,&useList,0,&chisq,&dof,0);
        xdbg<<"after outlier fit\n";

        xdbg<<"chisq,dof,sigma = "<<chisq<<','<<dof<<','<<
            sqrt(chisq/dof)<<std::endl;
        *outSigma = sqrt(chisq/dof);
    }

    void StarFinderAlgo::roughlyFitBrightStars(
        const std::vector<PotentialStar*>& objList,
        Function2D *f,double *outSigma)
    {
        // objList is already sorted by magnitude
        // make a new list with just the 50 brightest objects:
        std::vector<PotentialStar*> brightList(
            objList.begin(),
            objList.begin()+int(_nStart1*objList.size()));

        // now sort this list by size
        std::sort(brightList.begin(),brightList.end(),
                  std::mem_fun(&PotentialStar::isSmallerThan));

        xdbg<<"brightlist is:\n";
        const int twenty = std::min(20,int(brightList.size()));
        for(int i=0; i<twenty; ++i) {
            xdbg<<brightList[i]->getMag()<<" , "<<
                brightList[i]->getSize()<<std::endl;
        }
        if (brightList.size() > 20) {
            xdbg<<"...  (total of "<<brightList.size()<<" elements)\n";
        }
        // Of the smallest 20, find the 5 that form the tightest peak
        // Originally I hardcoded this number as 5, but now it is calculated
        // from _starFrac.
        int five = int(0.5*_starFrac*brightList.size());
        xdbg<<"'five' = "<<five<<std::endl;
        if (five < 5) {
            five = 5;
            xdbg<<"'five' => "<<five<<std::endl; 
        }
        if (3*five-1 >= int(brightList.size())) {
            std::ostringstream err;
            err<<"Too few objects in brightList.  Increase startn1.";
            throw StarFinderException(err.str());
        }
        int peakStart = 0;
        double peakWidth = brightList[five-1]->getSize()-brightList[0]->getSize();
        for(int k=1;k<=2*five;++k) {
            Assert(k+five-1<int(brightList.size()));
            double sizeRange = brightList[k+five-1]->getSize() -
                brightList[k]->getSize();
            if (sizeRange < peakWidth) {
                peakWidth = brightList[k+five-1]->getSize() -
                    brightList[k]->getSize();
                peakStart = k;
            }
        }
        xdbg<<"found peak starting at "<<peakStart<<
            " with width "<<peakWidth<<std::endl;

        int k1=peakStart, k2=peakStart+five;
        // Now expand list to be 2*binSize1 wide
        while (
            k1 > 0 && brightList[peakStart+2]->getSize() -
            brightList[k1-1]->getSize() < _binSize1) {
            --k1;
        }
        const int nBright = brightList.size();
        while (
            k2 < nBright && 
            brightList[k2]->getSize() - brightList[k1]->getSize() < 2.*_binSize1) {
            ++k2;
        }
        xdbg<<"expanded to "<<k1<<','<<k2<<std::endl;

        // Call these objects in the expanded peak stars
        std::vector<PotentialStar*> starList(
            brightList.begin()+k1, brightList.begin()+k2);

        // Fit f using a linear fit with 2 sigma clipping
        fitStellarSizes(f,0,2.,starList,outSigma);

        // if sigma is too big, drop either the first or last "star" and try again.
        while(*outSigma > _maxRms && starList.size() > 4) {
            double diff1 = starList.front()->getSize() -
                (*f)(starList.front()->getPos());
            double diff2 = starList.back()->getSize() -
                (*f)(starList.back()->getPos());

            if (std::abs(diff1) > std::abs(diff2)) {
                starList.erase(starList.begin());
                xdbg<<"Erased first star.  New size = "<<starList.size()<<std::endl;
            } else {
                starList.erase(starList.end()-1);
                xdbg<<"Erased last star.  New size = "<<starList.size()<<std::endl;
            }
            fitStellarSizes(f,0,2.,starList,outSigma);
        }
    }

}}}}
