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

#include <cassert>
#include <cmath>
#include <vector>

#include "lsst/meas/algorithms/shapelet/BinomFact.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    double fact(int i)
    {
        // return i!
        assert(i>=0);
        static std::vector<double> f(10);
        static bool first=true;

        if (first) {
            f[0] = f[1] = 1.;
            for(int j=2;j<10;j++) f[j] = f[j-1]*(double)j;
            first = false;
        }
        if (i>=(int)f.size()) {
#ifdef _OPENMP
#pragma omp critical (fact)
#endif
            {
                for(int j=f.size();j<=i;j++)
                    f.push_back(f[j-1]*(double)j);
            }
            assert(i==(int)f.size()-1);
        }
        assert(i<(int)f.size());
        return f[i];
    }

    double sqrtfact(int i)
    {
        // return sqrt(i!)
        static std::vector<double> f(10);
        static bool first=true;
        if (first) {
            f[0] = f[1] = 1.;
            for(int j=2;j<10;j++) f[j] = f[j-1]*sqrt((double)j);
            first = false;
        }
        if (i>=(int)f.size()) {
#ifdef _OPENMP
#pragma omp critical (sqrtfact)
#endif
            {
                for(int j=f.size();j<=i;j++)
                    f.push_back(f[j-1]*sqrt((double)j));
            }
        }
        assert(i<(int)f.size());
        return f[i];
    }

    double binom(int i,int j)
    {
        // return iCj, i!/(j!(i-j)!)
        static std::vector<std::vector<double> > f(10);
        static bool first=true;
        if (first) {
            f[0] = std::vector<double>(1,1.);
            f[1] = std::vector<double>(2,1.);
            for(int i1=2;i1<10;i1++) {
                f[i1] = std::vector<double>(i1+1);
                f[i1][0] = f[i1][i1] = 1.;
                for(int j1=1;j1<i1;j1++) f[i1][j1] = f[i1-1][j1-1] + f[i1-1][j1];
            }
            first = false;
        }
        if (i>=(int)f.size()) {
#ifdef _OPENMP
#pragma omp critical (binom)
#endif
            {
                for(int i1=f.size();i1<=i;i1++) {
                    f.push_back(std::vector<double>(i1+1,1.));
                    for(int j1=1;j1<i1;j1++) f[i1][j1] = f[i1-1][j1-1] + f[i1-1][j1];
                }
            }
            assert(i==(int)f.size()-1);
        }
        if (j<0 || j>i) return 0.;
        assert(i<(int)f.size());
        assert(j<(int)f[i].size());
        return f[i][j];
    }

    double sqrtn(int i)
    {
        // return sqrt(i)
        static std::vector<double> f(10);
        static bool first=true;
        if (first) {
            for(int j=0;j<10;j++) f[j] = std::sqrt((double)j);
            first = false;
        }
        if (i>=(int)f.size()) {
#ifdef _OPENMP
#pragma omp critical (sqrtn)
#endif
            {
                for(int j=f.size();j<=i;j++)
                    f.push_back(std::sqrt((double)j));
            }
        }
        assert(i<(int)f.size());
        return f[i];
    }



}}}}
