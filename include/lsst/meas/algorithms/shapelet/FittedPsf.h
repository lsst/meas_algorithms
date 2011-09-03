#ifndef FitPsf_H
#define FitPsf_H

#include <vector>
#include <iostream>
#include <memory>

#include "lsst/meas/algorithms/shapelet/BVec.h"
#include "lsst/meas/algorithms/shapelet/Bounds.h"
#include "lsst/meas/algorithms/shapelet/ConfigFile.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class FittedPsfAtXY;

    class FittedPsf {

    public :

#if 0
        // Make from PsfCatalog
        FittedPsf(PsfCatalog& psfcat, const ConfigFile& params, PsfLog& log);
#endif

        // Setup FittedPsf, but don't assign to the values yet.
        // Should be followed by either read or calculate.
        FittedPsf(const ConfigFile& params);

        int getPsfOrder() const { return _psfOrder; }
        int getPsfSize() const { return (_psfOrder+1)*(_psfOrder+2)/2; }
        int getFitOrder() const { return _fitOrder; }
        int getFitSize() const { return _fitSize; }
        int getNpca() const { return _nPca; }
        double getSigma() const { return _sigma; }
        void setSigma(double sigma) { _sigma = sigma; }

        double getXMin() const { return _bounds.getXMin(); }
        double getXMax() const { return _bounds.getXMax(); }
        double getYMin() const { return _bounds.getYMin(); }
        double getYMax() const { return _bounds.getYMax(); }
        const Bounds& getBounds() const { return _bounds; }

#if 0
        void write() const;
        void writeAscii(std::string file) const;
        void writeFits(std::string file) const;
        void writeFitsOld(std::string file) const;

        void read();
        void read(std::string file);
        void readAscii(std::string file);
        void readFits(std::string file);
#endif

        void calculate(
            const std::vector<Position>& pos,
            const std::vector<BVec>& psf,
            const std::vector<double>& nu,
            std::vector<long>& flags);

        void interpolate(Position pos, BVec& b) const
        {
            Assert(_avePsf.get());
            //Assert(_mV.get());
            Assert(b.getOrder() == _psfOrder);
            Assert(static_cast<int>(b.size()) == _avePsf->size());
            b.setSigma(_sigma);
            interpolateVector(pos,TMV_vview(b.vec()));
        }

        double interpolateSingleElement(Position pos, int i) const;

        // This next construct with FittedPsfAtXY allows you to write:
        // b = psf(pos);
        // instead of:
        // psf.interpolate(pos,b);
        // Both do the same thing.  I just like the first notation better.
        friend class FittedPsfAtXY;
        inline FittedPsfAtXY operator()(Position pos) const; // below...

        BVec getMean() const;

    private :

        void interpolateVector(Position pos, DVectorView b) const;

        const ConfigFile& _params;

        int _psfOrder;
        double _sigma;
        int _fitOrder;
        int _fitSize;
        int _nPca;
        Bounds _bounds;
        std::auto_ptr<DVector> _avePsf;
#if USE_TMV
        std::auto_ptr<tmv::Matrix<double,tmv::RowMajor> > _mV;
#else
        std::auto_ptr<DMatrix> _mV_transpose;
#endif
        std::auto_ptr<DMatrix> _f;
    };

    class FittedPsfAtXY : public AssignableToBVec
    {

    public :

        FittedPsfAtXY(const FittedPsf& psf, Position pos) :
            _psf(psf), _pos(pos) 
        {}

        size_t size() const { return _psf._avePsf->size(); }
        int getOrder() const { return _psf.getPsfOrder(); }
        double getSigma() const { return _psf.getSigma(); }

        void assignTo(BVec& b) const
        { _psf.interpolate(_pos,b); }

    private :

        const FittedPsf& _psf;
        Position _pos;
    };

    inline FittedPsfAtXY FittedPsf::operator()(Position pos) const
    { return FittedPsfAtXY(*this,pos); }

}}}}

#endif
