#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

# -*- lsst-python -*-

import lsst.meas.algorithms as measAlgorithms

def main():

    # Make a big elliptical annular coeff image.
    r1 = 30.0
    r2 = 50.0
    pa = (3.1415926/180) * 30.0
    ell = 0.7
    img = measAlgorithms.getCoeffImage(r1, r2, pa, ell);
    img.writeFits("cellip-%.1f-%.1f-%.1f-%.1f.fits" % (r1, r2, pa, ell))
                          

if __name__ == '__main__':
    main()
    
