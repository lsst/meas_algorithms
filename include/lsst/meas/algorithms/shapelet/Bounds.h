#ifndef MeasAlgoShapeletBounds_H
#define MeasAlgoShapeletBounds_H

#include <vector>
#include <complex>

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    class Position
    {

    public:

        Position() : _z(0.,0.) {}

        Position(const Position& rhs) : _z(rhs._z) {}

        Position(const std::complex<double>& rhs) : _z(rhs) {}

        ~Position() {}

        Position(double x, double y) : _z(x,y) {}

        Position& operator=(const Position& rhs) 
        { _z = rhs._z; return *this; }

        Position& operator=(const std::complex<double>& rhs) 
        { _z=rhs; return *this; }

        operator std::complex<double>() { return _z; }

        std::complex<double> operator-(const Position& p2) const 
        { return _z - p2._z; }

        std::complex<double> operator-(const std::complex<double>& z2) const 
        { return _z - z2; }

        Position& operator*=(double x) 
        { _z *= x; return *this; }

        Position& operator/=(double x) 
        { _z /= x; return *this; }

        double getX() const { return(_z.real()); }

        double getY() const { return(_z.imag()); }

        void read(std::istream& fin);

        void write(std::ostream& fout) const;

    private :

        std::complex<double> _z;

    }; // Position

    inline std::complex<double> operator-(
        const std::complex<double>& z1, const Position& p2)
    { return (p2-z1); }

    inline std::complex<double> operator/(
        const Position& p1, double x)
    { Position p2 = p1;  p2 /= x; return p2; }

    inline std::ostream& operator<<(std::ostream& os, const Position& pos)
    { pos.write(os); return os; }

    inline std::istream& operator>>(std::istream& os, Position& pos)
    { pos.read(os); return os; }

    class Bounds 
    {
        // Basically just a rectangle.  This is used to keep track of the bounds of
        // catalogs and fields.  You can set values, but generally you just keep
        // including positions of each galaxy or the bounds of each catalog
        // respectively using the += operators

    public:

        Bounds(double x1, double x2, double y1, double y2) :
            _isDefined(true),_xMin(x1),_xMax(x2),_yMin(y1),_yMax(y2) {}

        Bounds(const Position& pos) :
            _isDefined(true), _xMin(pos.getX()), _xMax(pos.getX()),
            _yMin(pos.getY()), _yMax(pos.getY()) {}

        Bounds(): _isDefined(false),_xMin(0.),_xMax(0.),_yMin(0.),_yMax(0.) {}

        ~Bounds() {}

        void setXMin(double x) { _xMin = x; _isDefined = true; }

        void setXMax(double x) { _xMax = x; _isDefined = true; }

        void setYMin(double y) { _yMin = y; _isDefined = true; }

        void setYMax(double y) { _yMax = y; _isDefined = true; }

        double getXMin() const { return _xMin; }

        double getXMax() const { return _xMax; }

        double getYMin() const { return _yMin; }

        double getYMax() const { return _yMax; }

        bool isDefined() const { return _isDefined; }

        Position getCenter() const;

        void operator+=(const Position& pos);

        void operator+=(const Bounds& rec);

        bool operator==(const Bounds& rhs) const;

        void snap(double d);

        void addBorder(double d);

        void addXBorder(double d);

        void addYBorder(double d);

        Bounds operator&(const Bounds& rhs) const; // Finds intersection

        bool includes(const Position& pos) const;

        bool includes(double x, double y) const;

        bool includes(const Bounds& b2) const;

        bool intersects(const Bounds& b2) const;

        double getArea() const;

        std::vector<Bounds> quarter() const;

        std::vector<Bounds> divide(int nx, int ny) const;

        void write(std::ostream& fout) const;

        void read(std::istream& fin);

        Position get00() const { return Position(_xMin,_yMin); }

        Position get01() const { return Position(_xMin,_yMax); }

        Position get10() const { return Position(_xMax,_yMin); }

        Position get11() const { return Position(_xMax,_yMax); }

        bool isWide() const { return (_xMax-_xMin > _yMax-_yMin); }

        bool isTall() const { return (_xMax-_xMin < _yMax-_yMin); }

    private:

        bool _isDefined;
        double _xMin,_xMax,_yMin,_yMax;

    };

    inline std::ostream& operator<<(std::ostream& fout, const Bounds& b)
    { b.write(fout); return fout;}

    inline std::istream& operator>>(std::istream& fin, Bounds& b)
    { b.read(fin); return fin;}


}}}}

#endif
