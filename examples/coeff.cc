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
 
// Compute a coefficient image for weighting of Sinc Aperture photometry
//
#include <iostream>
#include <sys/time.h>
#include "lsst/afw.h"
#include "lsst/meas/algorithms/detail/SincPhotometry.h"
#include "lsst/meas/algorithms/Photometry.h"

using namespace std;
namespace algorithms = lsst::meas::algorithms;
namespace image = lsst::afw::image;

class Timer {
  public:
    timeval start() {
        gettimeofday(&this->timer[0], NULL);
        return this->timer[0];
    }

    timeval stop() {
        gettimeofday(&this->timer[1], NULL);
        return this->timer[1];
    }

    int duration() const {
        int secs(this->timer[1].tv_sec - this->timer[0].tv_sec);
        int usecs(this->timer[1].tv_usec - this->timer[0].tv_usec);
        
        if(usecs < 0) {
            --secs;
            usecs += 1000000;
        }
        return static_cast<int>(secs * 1000 + usecs / 1000.0 + 0.5);
    }
    timeval timer[2];
};



void runAndPrint(int alg, double rad1, double rad2, double taper) {

    double posAngle = 0.0;
    double ellip = 0.0;
    
    Timer tm;
    tm.start();
    image::Image<double>::Ptr cimage;
    switch (alg) {
      case 1:
        cimage = algorithms::detail::calcImageRealSpace<double>(rad1, rad2, taper);
        break;
      case 2:
        cimage = algorithms::detail::calcImageKSpaceReal<double>(rad1, rad2);
        break;
      case 3:
        cimage = algorithms::detail::calcImageKSpaceCplx<double>(rad1, rad2, posAngle, ellip);
        break;
      default:
          lsst::afw::geom::ellipses::Axes const axes(rad2, rad2*(1.0-ellip), posAngle);
          cimage = algorithms::photometry::SincCoeffs<double>::calculate(axes, rad1/rad2);
        break;
    }
    tm.stop();

    printf("Computed coeff image with rad1=%.1f, rad2=%.1f, taper=%.2f in %.3f sec.\n",
           rad1, rad2, taper, 1.0*tm.duration()/1.0e3);
    char s[32];
    sprintf(s, "cimg%1d-%03.1f-%03.1f.fits", alg, rad1, rad2);
    cimage->writeFits(s);

}

int main(int argc, char *argv[]) {

    double rad1 = 4.0;
    double rad2 = 8.0;
    double taper = 0.1;
    double sigma = 1.5;
    if (argc == 4) {
        rad1 = atof(argv[1]);
        rad2 = atof(argv[2]);
        sigma = atof(argv[3]);
    }

    runAndPrint(0, rad1, rad2, taper);
    runAndPrint(1, rad1, rad2, taper);
    runAndPrint(2, rad1, rad2, taper);
    runAndPrint(3, rad1, rad2, taper);
    
}
