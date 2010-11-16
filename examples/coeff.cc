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



int main(int argc, char *argv[]) {

    double innerRadius = 4.0;
    double radius = 8.0;
    double taperwidth = 1.0;
    if (argc == 4) {
        innerRadius = atof(argv[1]);
        radius = atof(argv[2]);
        taperwidth = atof(argv[3]);
    }

    Timer tm1, tm2;

    tm1.start();
    image::Image<double>::Ptr cimage1 = algorithms::detail::getCoeffImage<double>(innerRadius, radius, taperwidth);
    tm1.stop();

    printf("Computed coeff image with innerRadius=%.1f, radius=%.1f, taper=%.2f in %.3f sec\n",
           innerRadius, radius, taperwidth, 1.0*tm1.duration()/1.0e3);
    char s[32];
    sprintf(s, "cimg-%03.1f-%03.1f-%04.2f.fits", innerRadius, radius, taperwidth);
    cimage1->writeFits(s);

    tm2.start();
    image::Image<double>::Ptr cimage2 = algorithms::detail::getCoeffImageFft(innerRadius, radius);
    tm2.stop();

    printf("Computed coeff image with innerRadius=%.1f, radius=%.1f, taper=%.2f in %.3f sec\n",
           innerRadius, radius, taperwidth, 1.0*tm2.duration()/1.0e3);
    sprintf(s, "cimg-%03.1f-%03.1f.fits", innerRadius, radius);
    cimage2->writeFits(s);

    
    
}
