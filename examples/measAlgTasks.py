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

from __future__ import print_function
from builtins import range
import os
import sys
import numpy as np

import lsst.utils
import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as measAlg
from lsst.meas.algorithms.detection import SourceDetectionTask
from lsst.meas.base import SingleFrameMeasurementTask


def loadData():
    """Prepare the data we need to run the example"""

    # Load sample input from disk
    afwdataDir = lsst.utils.getPackageDir('afwdata')
    imFile = os.path.join(afwdataDir, "CFHT", "D4", "cal-53535-i-797722_small_1.fits")
    exposure = afwImage.ExposureF(imFile)
    psf = measAlg.SingleGaussianPsf(21, 21, 2)
    exposure.setPsf(psf)

    im = exposure.getMaskedImage().getImage()
    im -= float(np.median(im.getArray()))

    return exposure


def run(display=False):
    exposure = loadData()
    schema = afwTable.SourceTable.makeMinimalSchema()
    #
    # Create the detection task
    #
    config = SourceDetectionTask.ConfigClass()
    config.thresholdPolarity = "both"
    config.background.isNanSafe = True
    config.thresholdValue = 3
    detectionTask = SourceDetectionTask(config=config, schema=schema)
    #
    # And the measurement Task
    #
    config = SingleFrameMeasurementTask.ConfigClass()

    config.algorithms.names = ["base_SdssCentroid", "base_SdssShape", "base_CircularApertureFlux"]
    config.algorithms["base_CircularApertureFlux"].radii = [1, 2, 4, 8, 16] # pixels

    config.slots.instFlux = None
    config.slots.modelFlux = None
    config.slots.psfFlux = None

    algMetadata = dafBase.PropertyList()
    measureTask = SingleFrameMeasurementTask(schema, algMetadata=algMetadata, config=config)
    radii = algMetadata.get("base_CircularApertureFlux_radii")
    #
    # Create the output table
    #
    tab = afwTable.SourceTable.make(schema)
    #
    # Process the data
    #
    result = detectionTask.run(tab, exposure)

    sources = result.sources

    print("Found %d sources (%d +ve, %d -ve)" % (len(sources), result.fpSets.numPos, result.fpSets.numNeg))

    measureTask.run(sources, exposure)
    if display:                         # display on ds9 (see also --debug argparse option)
        frame = 1
        ds9.mtv(exposure, frame=frame)

        with ds9.Buffering():
            for s in sources:
                xy = s.getCentroid()
                ds9.dot('+', *xy, ctype=ds9.CYAN if s.get("flags.negative") else ds9.GREEN, frame=frame)
                ds9.dot(s.getShape(), *xy, ctype=ds9.RED, frame=frame)

                for i in range(s.get("flux.aperture.nProfile")):
                    ds9.dot('o', *xy, size=radii[i], ctype=ds9.YELLOW, frame=frame)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Demonstrate the use of Source{Detection,Measurement}Task")

    parser.add_argument('--debug', '-d', action="store_true", help="Load debug.py?", default=False)
    parser.add_argument('--ds9', action="store_true", help="Display sources on ds9", default=False)

    args = parser.parse_args()

    if args.debug:
        try:
            import debug
        except ImportError as e:
            print(e, file=sys.stderr)

    run(display=args.ds9)
