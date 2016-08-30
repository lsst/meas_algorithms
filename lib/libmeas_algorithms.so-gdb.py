import os.path
import sys
import lsst.utils
import gdb
#
# Adjust the load path to include lsst.gdb, bypassing the regular lsstimport mechanism as
# the version of python running within gdb may not be the same as we are using for lsst processing
#

measAlgorithmsDir = lsst.utils.getPackageDir('meas_algorithms')
printerDir = os.path.join(measAlgorithmsDir, "python", "lsst", "gdb")
if printerDir not in sys.path:
    sys.path.append(printerDir)

import meas.algorithms.printers

meas.algorithms.printers.register(gdb.current_objfile())
