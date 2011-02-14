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

#include <vector>

#include "lsst/meas/algorithms/shapelet/Bounds.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    void Position::read(std::istream& fin)
    { double x,y; fin >> x >> y; _z = std::complex<double>(x,y); }

    void Position::write(std::ostream& fout) const
    { fout << getX() << " " << getY() <<" "; }

    void Bounds::operator+=(const Position& pos)
        // Expand the bounds to include the given position.
    {
        if (_isDefined) {
            if (pos.getX() < _xMin) _xMin = pos.getX();
            else if (pos.getX() > _xMax) _xMax = pos.getX();
            if (pos.getY() < _yMin) _yMin = pos.getY();
            else if (pos.getY() > _yMax) _yMax = pos.getY();
        } else {
            _xMin = _xMax = pos.getX();
            _yMin = _yMax = pos.getY();
            _isDefined = true;
        }
    }

    void Bounds::operator+=(const Bounds& rhs)
        // Expand the bounds to include the given bounds
    {
        if (!rhs.isDefined()) return;
        if (_isDefined) {
            if (rhs.getXMin() < _xMin) _xMin = rhs.getXMin();
            if (rhs.getXMax() > _xMax) _xMax = rhs.getXMax();
            if (rhs.getYMin() < _yMin) _yMin = rhs.getYMin();
            if (rhs.getYMax() > _yMax) _yMax = rhs.getYMax();
        } else {
            *this = rhs;
            _isDefined = true;
        }
    }

    bool Bounds::operator==(const Bounds& rhs) const
    {
        if (!_isDefined) return (!rhs._isDefined);
        else return (_xMin == rhs._xMin && _xMax == rhs._xMax &&
                     _yMin == rhs._yMin && _yMax == rhs._yMax);
    }


    Position Bounds::getCenter() const
    {
        return Position((_xMin + _xMax)/2.,(_yMin + _yMax)/2.);
    }

    Bounds Bounds::operator&(const Bounds& rhs) const
    {
        Bounds temp(_xMin<rhs._xMin ? rhs._xMin : _xMin,
                    _xMax>rhs._xMax ? rhs._xMax : _xMax,
                    _yMin<rhs._yMin ? rhs._yMin : _yMin,
                    _yMax>rhs._yMax ? rhs._yMax : _yMax);
        if (temp._xMin>temp._xMax || temp._yMin>temp._yMax) return Bounds();
        else return temp;
    }

    void Bounds::snap(double d)
    {
        if (_isDefined) {
            _xMax = d*ceil(_xMax/d);
            _xMin = d*floor(_xMin/d);
            _yMax = d*ceil(_yMax/d);
            _yMin = d*floor(_yMin/d);
        }
    }

    void Bounds::addXBorder(double d)
    {
        if (_isDefined) { _xMax += d; _xMin -= d; }
    }

    void Bounds::addYBorder(double d)
    {
        if (_isDefined) { _yMax += d; _yMin -= d; }
    }

    void Bounds::addBorder(double d)
    {
        if (_isDefined) { _xMax += d; _xMin -= d; _yMax += d; _yMin -= d; }
    }

    bool Bounds::includes(const Position& pos) const
    {
        return (_isDefined && pos.getX()<=_xMax && pos.getX()>=_xMin &&
                pos.getY()<=_yMax && pos.getY()>=_yMin);
    }

    bool Bounds::includes(double x, double y) const
    { return includes(Position(x,y)); }

    bool Bounds::includes(const Bounds& b2) const
    {
        if (!b2.isDefined()) return true;
        else return (
            _isDefined &&
            b2._xMin >= _xMin && b2._xMax <= _xMax &&
            b2._yMin >= _yMin && b2._yMax <= _yMax);
    }

    bool Bounds::intersects(const Bounds& b2) const
    {
        return (
            _isDefined && b2._isDefined &&
            !(b2._xMin >= _xMax) &&
            !(b2._xMax <= _xMin) &&
            !(b2._yMin >= _yMax) &&
            !(b2._yMax <= _yMin) );
    }

    double Bounds::getArea() const
    { return _isDefined ? (_xMax-_xMin)*(_yMax-_yMin) : 0.; }

    std::vector<Bounds> Bounds::quarter() const
    { return divide(2,2); }

    std::vector<Bounds> Bounds::divide(int nx, int ny) const
    {
        if (!_isDefined) return std::vector<Bounds>(nx*ny);
        std::vector<Bounds> temp;
        temp.reserve(nx*ny);
        std::vector<double> x(nx+1);
        std::vector<double> y(ny+1);
        x[0] = _xMin;  x[nx] = _xMax;
        y[0] = _yMin;  y[ny] = _yMax;
        double xstep = (_xMax-_xMin)/nx;
        double ystep = (_yMax-_yMin)/ny;
        for(int i=1;i<nx;++i) x[i] = x[0]+i*xstep;
        for(int j=1;j<ny;++j) y[j] = y[0]+j*ystep;
        for(int i=0;i<nx;++i) for(int j=0;j<ny;++j)
            temp.push_back(Bounds(x[i],x[i+1],y[j],y[j+1]));
        return temp;
    }

    void Bounds::write(std::ostream& fout) const
    { fout << _xMin << ' ' << _xMax << ' ' << _yMin << ' ' << _yMax << ' '; }

    void Bounds::read(std::istream& fin)
    { fin >> _xMin >> _xMax >> _yMin >> _yMax; _isDefined = true; }

}}}}
