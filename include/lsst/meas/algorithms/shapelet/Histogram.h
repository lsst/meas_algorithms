#ifndef MeasAlgoShapeletHistogram_H
#define MeasAlgoShapeletHistogram_H

#include <vector>
#include <iostream>

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    template <typename T>
    class Histogram {

    public:

        Histogram(double binSize, double minValue, double maxValue);
        ~Histogram() {}

        void add(double value,const T& ref);
        double findPeak(double minVal, double maxVal) const;
        bool hasSinglePeak(double minVal, double maxVal) const;
        double findValley(double minVal, double maxVal) const;
        double findFirstValueAfter(double start) const;
        double findFirstValleyAfter(double val1,bool hasPoissonNoise=false) const;
        double findFirstValleyBefore(double val1,bool hasPoissonNoise=false) const;
        double findFirstPeakAfter(double val1,bool hasPoissonNoise=false) const;
        double findFirstPeakBefore(double val1,bool hasPoissonNoise=false) const;
        int getTotalCountBefore(double val1) const;
        int getTotalCountAfter(double val1) const;
        int operator[](double value) const;
        double findThresh(double minVal, double maxVal) const;
        std::vector<T> getRefsInRange(double min, double max) const;
        std::vector<double> getValuesInRange(double min, double max) const;
        int getRefinedPeakCount(double* peak) const;
        int getRefinedValleyCount(double* valley) const;
        void print(std::ostream& fout,double val1=-1.e10,double val2=1.e10) const;

    private:

        int getCount(int i) const;
        int index(double value) const;
        double value(int i) const;

        double _binSize,_minValue,_maxValue;
        std::vector<std::vector<T> > _refs;
        std::vector<std::vector<double> > _values;
    };

}}}}

#endif
