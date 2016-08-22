#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
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
"""Example use of SubtractBackgroundTask
"""
from __future__ import absolute_import, division, print_function
import os.path

import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.meas.algorithms import SubtractBackgroundTask


def loadData():
    """Load the data we need to run the example"""
    mypath = lsst.utils.getPackageDir('afwdata')
    imFile = os.path.join(mypath, "CFHT", "D4", "cal-53535-i-797722_small_1.fits")
    return afwImage.ExposureF(imFile)


def run():
    """Subtract background
    """
    # create the task
    backgroundConfig = SubtractBackgroundTask.ConfigClass()
    backgroundTask = SubtractBackgroundTask(config=backgroundConfig)

    # load the data
    exposure = loadData()

    # subtract an initial estimate of background level
    bgRes = backgroundTask.run(exposure=exposure)
    background = bgRes.background

    # compute mean and variance of the background
    backgroundImage = background.getImage()
    s = afwMath.makeStatistics(backgroundImage, afwMath.MEAN | afwMath.VARIANCE)
    bgmean = s.getValue(afwMath.MEAN)
    bgvar = s.getValue(afwMath.VARIANCE)
    print("background mean=%0.1f; variance=%0.1f" % (bgmean, bgvar))


if __name__ == "__main__":
    run()
