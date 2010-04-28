#ifndef MeasAlgoShapeletPotentialStar_H
#define MeasAlgoShapeletPotentialStar_H

#include "lsst/meas/algorithms/shapelet/Bounds.h"
#include <string>

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class PotentialStar {

    public:
        PotentialStar(Position pos, double mag, double size, long index,
                      const std::string& line) :
            _pos(pos), _mag(mag), _size(size), _index(index), _line(line) 
        {}

        ~PotentialStar() {}

        const Position& getPos() const { return _pos; }

        double getMag() const { return _mag; }

        double getSize() const { return _size; }

        long getIndex() const { return _index; }

        const std::string& getLine() const { return _line; }

        void setSize(double newsize) { _size = newsize; }

        bool isBrighterThan(const PotentialStar* rhs) const
        { return _mag < rhs->_mag; }

        bool isSmallerThan(const PotentialStar* rhs) const
        { return _size < rhs->_size; }

    private:

        Position _pos;
        double _mag;
        double _size;
        long _index;
        std::string _line;

    };

}}}}

#endif

