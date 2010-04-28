#ifndef MeasAlgoShapeletEllipse_H
#define MeasAlgoShapeletEllipse_H

#include <complex>
#include <vector>
#include <ostream>
#include "lsst/meas/algorithms/shapelet/Pixel.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
//#include "lsst/meas/algorithms/shapelet/TimeVars.h"
#include "lsst/meas/algorithms/shapelet/MyMatrix.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class Ellipse 
    {

    public :

        Ellipse() :
            _cen(0.), _gamma(0.), _mu(0.),
            _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false), 
            _fPsf(1.0), _shouldDoTimings(false) {}

        Ellipse(double x, double y, double g1, double g2, double mu) :
            _cen(x,y), _gamma(g1,g2), _mu(mu), 
            _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false), 
            _fPsf(1.0), _shouldDoTimings(false) {}

        Ellipse(double vals[]) :
            _cen(vals[0],vals[1]), _gamma(vals[2],vals[3]), _mu(vals[4]),
            _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false),
            _fPsf(1.0), _shouldDoTimings(false) {}

        bool measure(const std::vector<PixelList>& pix, 
                     const std::vector<BVec>& psf,
                     int order, double sigma, bool shouldUseInteg, long& flag, 
                     DMatrix* cov=0, 
                     BVec* bret=0, DMatrix* bcov=0);
        bool measure(const std::vector<PixelList>& pix, 
                     int order, double sigma, bool shouldUseInteg, long& flag, 
                     DMatrix* cov=0, 
                     BVec* bret=0, DMatrix* bcov=0);

        void crudeMeasure(const PixelList& pix, double sigma);
        void crudeMeasure(
            const std::vector<PixelList>& pix, double sigma);

        void peakCentroid(const PixelList& pix, double maxr);

        void measureShapelet(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>& psf, BVec& bret, int order,
            DMatrix* bcov=0) const;
        void measureShapelet(
            const std::vector<PixelList>& pix, 
            BVec& bret, int order, DMatrix* bcov=0) const;

        void write(std::ostream& os) const
        { os << _cen<<" "<<_gamma<<" "<<_mu; }

        std::complex<double> getCen() const 
        { return _cen; }
        std::complex<double> getGamma() const 
        { return _gamma; }
        double getMu() const { return _mu; }

        void setCen(const std::complex<double>& cen)
        { _cen = cen; }
        void setGamma(const std::complex<double>& gamma)
        { _gamma = gamma; }
        void setMu(const double mu) { _mu = mu; }

        void fixCen() { _isFixedCen = true; }
        void fixGam() { _isFixedGamma = true; }
        void fixMu() { _isFixedMu = true; }

        void unfixCen() { _isFixedCen = false; }
        void unfixGam() { _isFixedGamma = false; }
        void unfixMu() { _isFixedMu = false; }

        bool isFixedCen() const { return _isFixedCen; }
        bool isFixedGamma() const { return _isFixedGamma; }
        bool isFixedMu() const { return _isFixedMu; }

        void setFP(double fp) { _fPsf = fp; }
        double getFP() const { return _fPsf; }

#if 0
        void doTimings() { _shouldDoTimings = true; }
        void resetTimes() { _times.reset(); }
        const EllipseTimes& getTimes() { return _times; }
#endif

    private :

        bool doMeasure(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>* psf, int order, double sigma,
            bool shouldUseInteg, long& flag, DMatrix* cov=0, 
            BVec* bret=0, DMatrix* bcov=0);
        void doMeasureShapelet(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>* psf, BVec& bret, int order,
            DMatrix* bcov=0) const;

        std::complex<double> _cen;
        std::complex<double> _gamma;
        double _mu; 

        bool _isFixedCen,_isFixedGamma,_isFixedMu;
        double _fPsf;

        bool _shouldDoTimings;

#if 0
        EllipseTimes _times;
#endif
    };

    inline std::ostream& operator<<(std::ostream& os, const Ellipse& s)
    { s.write(os); return os; }

}}}}

#endif
