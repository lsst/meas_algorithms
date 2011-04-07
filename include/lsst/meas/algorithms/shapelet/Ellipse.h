#ifndef MeasAlgoShapeletEllipse_H
#define MeasAlgoShapeletEllipse_H

#include <complex>
#include <vector>
#include <ostream>
#include "lsst/meas/algorithms/shapelet/Pixel.h"
#include "lsst/meas/algorithms/shapelet/BVec.h"
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
            _shouldDoTimings(false) {}

        Ellipse(std::complex<double> cen, std::complex<double> gamma,
                std::complex<double> mu) :
            _cen(cen), _gamma(gamma), _mu(mu), 
            _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false), 
            _shouldDoTimings(false) {}

        Ellipse(double vals[]) :
            _cen(vals[0],vals[1]), _gamma(vals[2],vals[3]), _mu(vals[4]),
            _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false),
            _shouldDoTimings(false) {}

        // Copy constructor and op= do not copy fixed-ness.  
        // They only copy the tranformation itself.
        Ellipse(const Ellipse& e2) :
            _cen(e2.getCen()), _gamma(e2.getGamma()),
            _mu(e2.getMu(),e2.getTheta()),
            _isFixedCen(false), _isFixedGamma(false), _isFixedMu(false),
            _shouldDoTimings(false) {}

        Ellipse& operator=(const Ellipse& e2)
        { 
            _cen = e2.getCen();
            _gamma = e2.getGamma();
            _mu = std::complex<double>(e2.getMu(),e2.getTheta());
            return *this;
        }

        // Find the ellipse transformation in which the galaxy is observed
        // to be "round", where round means:
        // b10 = 0, so galaxy is centroided, 
        //          iff _isFixedCen == false
        // b11 = 0, so shapelet is using optimal sigma, 
        //          iff _isFixedMu == false, and not deconvolving by psf.
        // b20 = 0, so galaxy is round, 
        //          iff _isFixedGamma == false
        // The first one below is for the galaxy as observed, and the 
        // second is for the galaxy before being convolved by the psf.
        bool measure(
            const std::vector<PixelList>& pix, 
            int order, int order2, int maxm,
            double sigma, long& flag, double thresh, DMatrix* cov=0);
        bool measure(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>& psf,
            int order, int order2, int maxm,
            double sigma, long& flag, double thresh, DMatrix* cov=0);

        // Given a measured shapelet vector, find the ellipse transformation
        // in which this shapelet vector would be observed to be round.
        bool findRoundFrame(
            const BVec& b, bool psf, int galOrder2, double thresh,
            long& flag, DMatrix* cov=0);

        // Do a really simple and fast measurement to get a good starting point.
        // Only does centroid and size.
        void crudeMeasure(const PixelList& pix, double sigma);
        void crudeMeasure(const std::vector<PixelList>& pix, double sigma);

        // Even simpler starting point, just doing the centroid.
        void peakCentroid(const PixelList& pix, double maxr);

        // Measure the shapelet decomposition in the frame described by 
        // the current Ellipse parameters.
        // The order parameter is allowed to be less than the order of bret,
        // in which case the rest of the vector is set to zeros.
        // order2 is the order for the intermediate steps in the calculation.
        // Again, the first is for the galaxy as observed (the "native" fit),
        // and the second is for the galaxy before being convolved by the psf.
        bool measureShapelet(
            const std::vector<PixelList>& pix, BVec& bret,
            int order, int order2, int maxm, DMatrix* bcov=0) const;
        bool measureShapelet(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>& psf, BVec& bret,
            int order, int order2, int maxm, DMatrix* bcov=0) const;

        // An alternative measurement that uses the formula:
        // <psi_pq | psi_st> = delta_ps delta_qt
        // Since this formulat is only really accurate in the limit of an
        // infinite field and infinitely small pixels, it is not really useful
        // for real data.  However, it is useful for testing purposes.
        bool altMeasureShapelet(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>& psf, BVec& bret, int order, int order2,
            double pixScale, DMatrix* bcov=0) const;
        bool altMeasureShapelet(
            const std::vector<PixelList>& pix, BVec& bret, int order, int order2,
            double pixScale, DMatrix* bcov=0) const;

        // Find the effective Ellipse transformation from combining the 
        // current transformation with a second one, either before or after.
        void preShiftBy(const std::complex<double>& cen,
                        const std::complex<double>& gamma,
                        const std::complex<double>& mu);
        void postShiftBy(const std::complex<double>& cen,
                         const std::complex<double>& gamma,
                         const std::complex<double>& mu);

        // Find the rotation-free transformation that distorts a 
        // circle to the same final shape as the current transformation does.
        // In otherwords, we move the rotation part from the end of the 
        // transformation to the beginning, since a rotation of a circle 
        // is a null operation.
        // After doing this method, getTheta() will return 0.
        void removeRotation();

        std::complex<double> getCen() const { return _cen; }
        std::complex<double> getGamma() const { return _gamma; }
        double getMu() const { return real(_mu); }
        double getTheta() const { return imag(_mu); }

        void setCen(const std::complex<double>& cen) { _cen = cen; }
        void setGamma(const std::complex<double>& gamma) { _gamma = gamma; }
        void setMu(double mu) { _mu = std::complex<double>(mu,getTheta()); }
        void setTheta(double theta) { _mu = std::complex<double>(getMu(),theta); }

        void fixCen() { _isFixedCen = true; }
        void fixGam() { _isFixedGamma = true; }
        void fixMu() { _isFixedMu = true; }

        void unfixCen() { _isFixedCen = false; }
        void unfixGam() { _isFixedGamma = false; }
        void unfixMu() { _isFixedMu = false; }

        bool isFixedCen() const { return _isFixedCen; }
        bool isFixedGamma() const { return _isFixedGamma; }
        bool isFixedMu() const { return _isFixedMu; }

        void doTimings() { _shouldDoTimings = true; }

        void write(std::ostream& os) const
        { os << _cen<<" "<<_gamma<<" "<<_mu; }

    private :

        bool doMeasure(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>* psf, int order, int order2, int maxm,
            double sigma, long& flag, double thresh, DMatrix* cov=0);

        bool doMeasureShapelet(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>* psf, BVec& bret,
            int order, int order2, int maxm, DMatrix* bcov=0) const;

        bool doAltMeasureShapelet(
            const std::vector<PixelList>& pix, 
            const std::vector<BVec>* psf, BVec& bret,
            int order, int order2, double pixScale, DMatrix* bcov=0) const;

        std::complex<double> _cen;
        std::complex<double> _gamma;
        std::complex<double> _mu; 

        bool _isFixedCen,_isFixedGamma,_isFixedMu;

        bool _shouldDoTimings;

    };

    inline std::ostream& operator<<(std::ostream& os, const Ellipse& s)
    { s.write(os); return os; }

    // The Miralda-Escude formula for adding two shears.
    // We transform to distortion first, then use his formula, then 
    // transform back to shear.
    std::complex<double> addShears(
        const std::complex<double> g1, const std::complex<double> g2);

    // Find the net Ellipse parameters from doing (c1,g1,m1) first,
    // then (c2,g2,m2).
    // The final transformation is returned as (c3,g3,m3).
    void CalculateShift(
        std::complex<double> c1, std::complex<double> g1, std::complex<double> m1,
        std::complex<double> c2, std::complex<double> g2, std::complex<double> m2,
        std::complex<double>& c3, std::complex<double>& g3, 
        std::complex<double>& m3);

}}}}

#endif
