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

#include <algorithm>

#include "lsst/meas/algorithms/shapelet/Histogram.h"
#include "lsst/meas/algorithms/shapelet/Form.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"
#include "lsst/meas/algorithms/shapelet/PotentialStar.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    template <typename T>
    Histogram<T>::Histogram(double binSize, double minValue, double maxValue) :
        _binSize(binSize)
    {
        // Add an extra bin to prevent rounding errors from giving you a
        // value outsize the allowed range.
        int nBins = int(ceil((maxValue-minValue)/_binSize))+1;
        _minValue = minValue - binSize/2.;
        _maxValue = maxValue + binSize/2.;
        _refs = std::vector<std::vector<T> >(nBins);
        _values = std::vector<std::vector<double> >(nBins);
        if (value(nBins-1) > _maxValue) _maxValue = value(nBins-1)+binSize/4.;

        dbg<<"made histogram:\n";
        dbg<<"minvalue="<<_minValue<<
            " index(minvalue)="<<index(_minValue)<<std::endl;
        dbg<<"maxvalue="<<_maxValue<<
            " index(maxvalue)="<<index(_maxValue)<<std::endl;
        dbg<<"mini=0  value(mini)="<<value(0)<<std::endl;
        dbg<<"maxi="<<_refs.size()-1;
        dbg<<" value(maxi)="<<value(_refs.size()-1)<<std::endl;
    }

    template <typename T>
    void Histogram<T>::add(double value,const T& ref)
    {
        _refs[index(value)].push_back(ref);
        _values[index(value)].push_back(value);
    }

    template <typename T>
    int Histogram<T>::getCount(int i) const
    {
        Assert(i<int(_refs.size()));
        Assert(_values[i].size() == _refs[i].size());
        return _refs[i].size();
    }

    template <typename T>
    double Histogram<T>::findPeak(double val1, double val2) const
    {
        if (val1<_minValue) val1 = _minValue;
        if (val2>_maxValue) val2 = _maxValue;
        int i1=index(val1), i2=index(val2);
        Assert(i1<int(_refs.size()));
        Assert(i2<int(_refs.size()));
        int iPeak=i1;
        int maxCount = getCount(i1);

        for(int i=i1+1;i<=i2;++i) {
            if (getCount(i) > maxCount) {
                maxCount = getCount(i);
                iPeak = i;
            }
        }
        return value(iPeak);
    }

    template <typename T>
    bool Histogram<T>::hasSinglePeak(double val1, double val2) const
    // Finds the first peak from the left... set to iPeak1
    // and the first peak from the right... set to iPeak2
    // Returns whether they are the same peak
    {
        if (val1<_minValue) val1 = _minValue;
        if (val2>_maxValue) val2 = _maxValue;
        int i1=index(val1), i2=index(val2);
        Assert(i1<int(_refs.size()));
        Assert(i2<int(_refs.size()));
        int iPeak1=i1,iPeak2=i2;
        while(iPeak1+1<=i2 && getCount(iPeak1+1)>=getCount(iPeak1)) ++iPeak1;
        while(iPeak2-1>=i1 && getCount(iPeak2-1)>=getCount(iPeak2)) --iPeak2;
        return iPeak1 >= iPeak2; // Usually ==, but if plateau, they swap
    }

    template <typename T>
    double Histogram<T>::findValley(double val1, double val2) const
    {
        if (val1<_minValue) val1 = _minValue;
        if (val2>_maxValue) val2 = _maxValue;
        int i1=index(val1), i2=index(val2);
        Assert(i1<int(_refs.size()));
        Assert(i2<int(_refs.size()));
        int iValley=i1;
        int minCount = getCount(i1);

        for(int i=i1+1;i<=i2;++i) {
            if (getCount(i) < minCount) {
                minCount = getCount(i);
                iValley = i;
            }
        }
        // Make sure you're not still sloping down.  If so, keep going.
        if (iValley==i2) {
            int nRefs = _refs.size();
            while(iValley+1 < nRefs && getCount(iValley+1) < getCount(iValley)) 
                ++iValley;
        }

        return value(iValley);
    }

    template <typename T>
    double Histogram<T>::findFirstValleyAfter(double val1,bool hasPoissonNoise) const
    {
        dbg<<"Start FindFirstValleyAfter "<<val1<<std::endl;
        if (val1<_minValue) val1 = _minValue;
        int i1=index(val1);
        Assert(i1<int(_refs.size()));
        // Keep going until it stops dropping
        // The sqrt bit allows for Poisson noise in the counts
        // Otherwise you can get caught early by ledges.
        // eg. 15 9 10 4 1 3 3 ... would stop at the 9 instead of at the 1
        int iValley=i1;
        double minCount = getCount(i1);
        int iCheck=i1+1;
        const int nRefs = _refs.size();
        while(iCheck < nRefs && getCount(iCheck) <= minCount +
              (hasPoissonNoise ? int(sqrt(float(getCount(iCheck)))) : 0)) {
            if (getCount(iCheck) < minCount) {
                minCount = getCount(iCheck);
                iValley = iCheck;
            }
            if (++iCheck == nRefs) break;
        }
        dbg<<"valley = "<<value(iValley)<<std::endl;
        return value(iValley);
    }

    template <typename T>
    double Histogram<T>::findFirstValleyBefore(double val1,bool hasPoissonNoise) const
    {
        if (val1>_maxValue) val1 = _maxValue;
        int i1=index(val1);
        Assert(i1<int(_refs.size()));
        // Keep going until it stops dropping
        // The sqrt bit allows for Poisson noise in the counts
        // Otherwise you can get caught early by ledges.
        // eg. 15 9 10 4 1 3 3 ... would stop at the 9 instead of at the 1
        int iValley=i1;
        double minCount = getCount(i1);
        if (i1 > 0) {
            int iCheck=i1-1;
            while(getCount(iCheck) <= minCount +
                  (hasPoissonNoise ? int(sqrt(float(getCount(iCheck)))) : 0)) {
                if (getCount(iCheck) < minCount) {
                    minCount = getCount(iCheck);
                    iValley = iCheck;
                }
                if (iCheck == 0) break;
                --iCheck;
            }
        }
        return value(iValley);
    }

    template <typename T>
    double Histogram<T>::findFirstValueAfter(double start) const
    {
        int i = (start<_minValue) ? index(_minValue) : index(start)+1;
        const int nValues = _values.size();
        while (i<nValues && _values[i].size() == 0) ++i;
        if (i == nValues) return _maxValue;
        else return *std::min_element(_values[i].begin(),_values[i].end());
    }

    template <typename T>
    double Histogram<T>::findFirstPeakAfter(double val1,bool hasPoissonNoise) const
    {
        if (val1<_minValue) val1 = _minValue;
        int i1=index(val1);
        const int nRefs = _refs.size();
        Assert(i1 < nRefs);
        // Keep going until it stops rising
        int iPeak=i1;
        double maxCount = getCount(i1);
        int iCheck=i1+1;
        while(iCheck < nRefs && maxCount <= getCount(iCheck) +
              (hasPoissonNoise ? int(sqrt(float(getCount(iCheck)))) : 0)) {
            if (getCount(iCheck) > maxCount) {
                maxCount = getCount(iCheck);
                iPeak = iCheck;
            }
            if (++iCheck == nRefs) break;
        }
        return value(iPeak);
    }

    template <typename T>
    double Histogram<T>::findFirstPeakBefore(double val1,bool hasPoissonNoise) const
    {
        if (val1>_maxValue) val1 = _maxValue;
        int i1=index(val1);
        Assert(i1 < int(_refs.size()));
        // Keep going until it stops rising
        int iPeak=i1;
        double maxCount = getCount(i1);
        if (i1 > 0) {
            int iCheck=i1-1;
            while(maxCount <= getCount(iCheck) +
                  (hasPoissonNoise ? int(sqrt(float(getCount(iCheck)))) : 0)) {
                if (getCount(iCheck) > maxCount) {
                    maxCount = getCount(iCheck);
                    iPeak = iCheck;
                }
                if (iCheck == 0) break;
                --iCheck;
            }
        }
        return value(iPeak);
    }

    template <typename T>
    int Histogram<T>::getTotalCountBefore(double val1) const
    {
        int count=0;
        int i1 = index(val1);
        Assert(i1 < int(_refs.size()));
        for(int i=0;i<i1;++i) count += getCount(i);
        return count;
    }

    template <typename T>
    int Histogram<T>::getTotalCountAfter(double val1) const
    {
        int count=0;
        int i1 = index(val1);
        const int nRefs = _refs.size();
        Assert(i1 < nRefs);
        for(int i=i1+1;i<nRefs;++i) count += getCount(i);
        return count;
    }

    template <typename T>
    int Histogram<T>::getRefinedPeakCount(double* peak) const
    // Takes the three bins around the peak and find the maximum number of
    // objects which fit into a bin width.  
    // It does this by sliding a window of width binSize along a sorted list
    // of values and counting how many objects fit into the window. 
    {
        std::vector<double> vals = 
            getValuesInRange(*peak-1.5*_binSize,
                             *peak+1.5*_binSize);
        Assert(vals.size()>=1);
        std::sort(vals.begin(),vals.end());
        int bestCount=0;
        double bestVal=*peak;
        // i is the object at the left side of the virtual bin
        // j is the first object after the end of the virtual bin
        const int nVals = vals.size();
        for(int i=0,j=0;j<nVals;++i) {
            // Determine j for this value of i
            while (j<nVals && vals[j] < vals[i]+_binSize) ++j;
            // If the corresponding count is better than bestCount, update it
            Assert(j>i);
            if (j-i > bestCount) {
                bestCount = j-i;
                bestVal = (vals[i]+vals[j-1])/2.;
            }
        }
        *peak = bestVal;
        return bestCount;
    }

    template <typename T>
    int Histogram<T>::getRefinedValleyCount(double* valley) const
    // Just like getRefinedPeakCount, but for a valley.
    // One slight difference is that for the best peak, you want to 
    // include as many objects as possible, so the "i" value is included
    // For the best valley, we imagine that the "i" value is just slightly
    // to the left of the start of the bin.  So the count is j-i-1 rather
    // than j-i.
    {
        std::vector<double> vals = 
            getValuesInRange(*valley-1.5*_binSize,
                             *valley+1.5*_binSize);
        const int nVals = vals.size();
        if (nVals == 0) return 0;
        std::sort(vals.begin(),vals.end());
        // start high since want to find the lowest count
        int bestCount=nVals;
        double bestVal=*valley;
        // i is the last object before the start of the virtual bin
        // j is the first object after the end of the virtual bin
        for(int i=0,j=0; i<nVals && vals[i]<*valley+0.5*_binSize; ++i) {
            // Determine j for this value of i
            while (j<nVals && vals[j] < vals[i]+_binSize) ++j;
            // If we can get rid of any i values, given this j, do so.
            double nextBinStart = j<nVals ? vals[j] : *valley+1.5*_binSize;
            while (i+1<nVals && vals[i+1]+_binSize < nextBinStart) ++i;
            // If the corresponding count is better than bestCount, update it
            Assert(j>i);
            if (j-i-1 < bestCount) {
                bestCount = j-i-1;
                bestVal = vals[i]+_binSize/2.;
            }
        }
        *valley = bestVal;
        return bestCount;
    }

    template <typename T>
    int Histogram<T>::operator[](double value) const
    {
        if (index(value) == int(_refs.size())) return 0;
        Assert(index(value) < int(_refs.size()));
        return getCount(index(value));
    }

    template <typename T>
    double Histogram<T>::findThresh(double minVal, double maxVal) const
    // This function finds the threshold point which maximizes the 
    // "Normalized Between Group Variance" (var).
    // For a given threshold point, we have two groups, the left side and 
    // the right side. 
    // Define fL and fR to be the fraction of points in each group.
    // xL and xR are the centroids of each group (the mean value)
    // vL and vR are the variances for each group.
    // Then the Between Group Variance is just:
    // bgv = fL fR (xR - xL)^2
    // Maximizing this works well for distributions which may have different
    // numbers of elements so long as the widths of the two distributions
    // are similar.
    // When the widths are not similar, it is better to normalize this by the
    // effective sigmas, or sqrt(variance).
    // var = fL fR (xR-xL)^2 / sqrt(vL vR)
    // Maximizing this seems to work well for separating stars from galaxies.
    //
    // To calculate this efficiently, we want to be able to keep track of all
    // the values as we go along the histogram to make it order N, rather than
    // N^2.
    //
    // nTot = Sumtot(N(i))
    // meanTot = Sumtot(N(i)*i)/Ntot
    // meanSqTot = Sumtot(N(i)*i*i)/Ntot
    //
    // (All sums below are from i=0..i1, "sumtot" above means from i=0..maxi)
    //
    // fL = Sum(N(i))/nTot
    // fR = 1-fL
    // xL = Sum(N(i)*i)/nTot/fL
    // xR = (meanTot - fL*xL)/fR
    // vL = Sum(N(i)*i*i)/nTot/fL - xL^2
    // vR = (meanSqTot - fL*(vL+xL^2))/fR
    {
        const double EPSILON = 1.e-6;

        dbg<<"starting FindThresh\n";
        if (dbgout) print(*dbgout);

        double nTot=0.,meanTot=0.,meanSqTot=0.;

        if (minVal < _minValue) minVal = _minValue;
        if (maxVal > _maxValue) maxVal = _maxValue;
        int i1 = index(minVal);
        int i2 = index(maxVal);
        Assert(i1 < int(_refs.size()));
        Assert(i2 < int(_refs.size()));

        // First calculate the "tot" values:
        for(int i=i1;i<=i2;++i) {
            double dcount = getCount(i);
            nTot += dcount;
            meanTot += dcount*i;
            meanSqTot += dcount*i*i;
        }
        meanTot /= nTot;
        meanSqTot /= nTot;
        dbg<<"ntot="<<nTot<<", meantot="<<meanTot<<
            ", meansqtot="<<meanSqTot<<std::endl;

        double sumN=0.,sumNi=0.,sumNii=0.,varBest=-1.;
        int iBest=-1;

        const double minVariance = 0.25;

        for(int i=i1;i<=i2;++i) {
            double dcount = getCount(i);
            sumN += dcount;
            sumNi += dcount*i;
            sumNii += dcount*i*i;
            double fL = sumN/nTot;
            double fR = 1.-fL;
            Assert(fL >=0. && fR >=0.);
            if (fL == 0. || fR == 0.) continue; // var = 0
            double xL = sumNi/nTot;  // actuall fL*xL now
            double xR = (meanTot-xL)/fR;
            xL /= fL;
            double vL = sumNii/nTot; // actually fL*(vL+xL^2)
            double vR = (meanSqTot-vL)/fR - xR*xR;
            vL = vL/fL - xL*xL;
            if (vL<minVariance) vL = minVariance;
            if (vR<minVariance) vR = minVariance;
            double var = fL*fR*std::pow(xR-xL,2)/sqrt(vL*vR);
            Assert(var >= 0.);
            if (var > varBest+EPSILON) {
                iBest = i;
                varBest = var;
            }
        }

        dbg<<"besti = "<<iBest<<", bestvar = "<<varBest<<std::endl;
        Assert(iBest >= 0);
        Assert(varBest > 0.);
        // +1 below forces thresh to be at least 1 bin away from first peak.
        // otherwise if peak is only 1 bin, can get iBest = that bin.
        dbg<<"returning thresh = "<<value(iBest)<<std::endl;
        return value(iBest);
    }

    template <typename T>
    std::vector<T> Histogram<T>::getRefsInRange(double min, double max) const
    {
        if (min < _minValue) min = _minValue;
        if (max > _maxValue) max = _maxValue;
        int i1=index(min);
        int i2=index(max);

        dbg<<"in getrefs: min,max = "<<min<<','<<max<<std::endl;

        std::vector<T> temp;
        const int nRefs1 = _refs[i1].size();
        for(int k=0; k<nRefs1; ++k) {
            if(_values[i1][k]>=min) {
                dbg<<"i1 - add ref @ "<<_values[i1][k]<<std::endl;
                temp.push_back(_refs[i1][k]);
            }
        }
        for(int i=i1+1;i<i2;++i) {
            if (_refs[i].size() > 0) 
                dbg<<"i - add all ("<<_refs[i].size()<<") refs near "<<
                    _values[i].front()<<std::endl;
            temp.insert(temp.end(),_refs[i].begin(),_refs[i].end());
        }
        const int nRefs2 = _refs[i2].size();
        for(int k=0; k<nRefs2; ++k) {
            if(_values[i2][k] <= max) {
                dbg<<"i2 - add ref @ "<<_values[i2][k]<<std::endl;
                temp.push_back(_refs[i2][k]);
            }
        }
        return temp;
    }

    template <typename T>
    std::vector<double> Histogram<T>::getValuesInRange(
        double min, double max) const
    {
        if (min < _minValue) min = _minValue;
        if (max > _maxValue) max = _maxValue;
        int i1=index(min);
        int i2=index(max);

        std::vector<double> temp;
        const int nVals1 = _values[i1].size();
        for(int k=0;k<nVals1;++k)
            if(_values[i1][k]>=min) temp.push_back(_values[i1][k]);
        for(int i=i1+1;i<i2;++i)
            temp.insert(temp.end(),_values[i].begin(),_values[i].end());
        const int nVals2 = _values[i2].size();
        for(int k=0;k<nVals2;++k)
            if(_values[i2][k]<=max) temp.push_back(_values[i2][k]);
        return temp;
    }

    template <typename T>
    int Histogram<T>::index(double value) const
    {
        Assert (value >= _minValue && value <= _maxValue);
        return int(floor((value - _minValue)/_binSize));
    }

    template <typename T>
    double Histogram<T>::value(int i) const
    {
        Assert(i<int(_refs.size()));
        return _binSize*(i+0.5)+_minValue;
    }

    template <typename T>
    void Histogram<T>::print(std::ostream& fout,double val1,double val2) const
    {
        if (val1 < _minValue) val1 = _minValue;
        if (val2 > _maxValue) val2 = _maxValue;
        int i1=index(val1);
        int i2=index(val2);
        Assert(i1 < int(_refs.size()));
        Assert(i2 < int(_refs.size()));
        Form sci2; sci2.sci().width(10).prec(2);

        for(int i=i1;i<=i2;++i) {
            //if ((i-i1)%5 == 0) 
            fout << sci2(value(i));
            //else fout << "           ";
            for(int j=0;j<getCount(i);++j) fout<<'*';
            fout<<std::endl;
        }
    }

    template class Histogram<PotentialStar*>;

}}}}
