# 
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
"""lsst.meas.algorithms
"""
from __future__ import absolute_import

import lsst.afw.image
import lsst.afw.math

from .binnedWcs import *
from .cr import *
from .coaddBoundedField import *
from .imagePsf import *
from .interp import *
from .kernelPsf import *
from .psf import *
from .pcaPsf import *
from .psfCandidate import * #python
from .singleGaussianPsf import *
from .spatialModelPsf import *
from .warpedPsf import *
from .coaddPsf import *
from .doubleGaussianPsf import *

from .defects import *
from .psfDeterminer import *
from .pcaPsfDeterminer import *
from .starSelector import *
from .findCosmicRaysConfig import *
from .detection import *
from .gaussianPsfFactory import *
from .loadReferenceObjects import *
from .secondMomentStarSelector import *
from .objectSizeStarSelector import *
from .makeCoaddApCorrMap import *
from .subtractBackground import *
from .measureApCorr import *
from .flaggedSourceSelector import *
from .sourceSelector import *
from .astrometrySourceSelector import *
from .matcherSourceSelector import *
from .ingestIndexReferenceTask import *
from .loadIndexedReferenceObjects import *
from .indexerRegistry import *

from .version import *

import lsst.utils

for name in dict(globals()):
    if name.endswith("_swigregister"):
        del globals()[name]
