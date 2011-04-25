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
import math
from lsst.pex.logging import Log, LogRec, endr


class SourceClassifier(object):
    """
    Base class for source classifiers. A SourceClassifier is initialized with a policy, and the positions of
    the flag bits it is allowed to set/clear. Once created, classify() is called on the SourceClassifier instance
    some number of times - each time one or more sources (typically a pair of measurements on the two exposures
    in an LSST visit) are supplied as arguments.

    Finally, the finish() method is called on the classifier so that it can log summary statistics and/or
    place SDQA ratings onto a stage clipboard.
    """
    def __init__(self, bits, policy):
        self._bits = bits
        self._policy = policy

    def classify(self, *args):
        """
        Classify one or more sources. To be overriden by subclasses.
        """
        pass

    def finish(self, log, clipboard):
        """
        Lifecycle method to be overriden by subclasses. Called once by SourceClassificationStage
        after the SourceClassifier instance has been used to classify one visits worth of data. Can
        be used to log visit summary statistics or place SDQA ratings onto the stage clipboard.
        """
        pass

    def getPolicy(self):
        """
        Return the policy containing classifier configuration parameters
        """
        return self._policy

    def getMask(self, n=0):
        """
        Return a mask for the n-th bit this classifier sets/clears
        """
        return 1 << self._bits[n]

    def getBit(self, n=0):
        """
        Return the index of the n-th bit this classifier sets/clears
        """
        return self._bits[n]

    def setBit(self, flag, n=0):
        """
        Set the n-th classifier bit in the given integer
        """
        return flag | self.getMask(n)

    def clearBit(self, flag, n=0):
        """
        Clear the n-th classifier bit in the given integer
        """
        return flag ^ (flag & self.getMask(n))


class PresentInBothExposuresClassifier(SourceClassifier):
    """
    Check whether the absolute values of the PSF fluxes of two sources are both above a threshold.

    Sets two bits in the classification flags for each source:
    - bit 0 is set when a source exceeds the flux threshold in the
      exposure it was measured on.
    - bit 1 is set when the sibling source (measured on the other
      exposure in the visit) exceeds the flux threshold.

    Configuration parameters:
        "psfFluxThreshold": if the absolute value of the PSF flux of a source is below this
                            threshold, the source is considered to be missing from an exposure
                            (this happens when source detection and measurement happen on
                            different exposures).
    """
    def __init__(self, bits, policy):
        SourceClassifier.__init__(self, bits, policy)
        self._numPresentInBoth = 0
        self._numPresentOnlyInFirst = 0
        self._numPresentOnlyInSecond = 0
        self._numMissingInBoth = 0
        self._psfFluxThreshold = policy.getDouble("psfFluxThreshold")
        assert self._psfFluxThreshold > 0

    def classify(self, *args):
        flag0 = args[0].getFlagClassification()
        flag1 = args[1].getFlagClassification()
        inFirst = abs(args[0].getPsfFlux()) > self._psfFluxThreshold
        inSecond = abs(args[1].getPsfFlux()) > self._psfFluxThreshold
        if inFirst:
            flag0 = self.setBit(flag0, 0)
            flag1 = self.setBit(flag1, 1)
            if inSecond: self._numPresentInBoth += 1
            else:        self._numPresentOnlyInFirst += 1
        else:
            flag0 = self.clearBit(flag0, 0)
            flag1 = self.clearBit(flag1, 1)
            if inSecond: self._numPresentOnlyInSecond += 1
            else:        self._numMissingInBoth += 1
        if inSecond:
            flag0 = self.setBit(flag0, 1)
            flag1 = self.setBit(flag1, 0)
        else:
            flag0 = self.clearBit(flag0, 1)
            flag1 = self.clearBit(flag1, 0)
        args[0].setFlagClassification(flag0)
        args[1].setFlagClassification(flag1)

    def finish(self, log=None, clipboard=None):
        if log:
            rec = LogRec(log, Log.INFO)
            rec << "PresentInBothExposuresClassifier visit statistics"
            rec << { "numPresentInBoth":       self._numPresentInBoth,
                     "numPresentOnlyInFirst":  self._numPresentOnlyInFirst,
                     "numPresentOnlyInSecond": self._numPresentOnlyInSecond,
                     "numMissingInBoth":       self._numMissingInBoth }
            rec << endr


class ShapeDiffersInExposuresClassifier(SourceClassifier):
    """
    Check whether the shapes of two difference sources differ significantly.
    Probably bogus and needs proof reading.

    Sets a single classification flag bit for both sources:
    - bit 0 is set when the two sources differ in shape

    Configuration parameters:
        "shapeNormDiffThreshold": if the absolute value of the difference of the shape
                                  parameter norms for two source is greater than this
                                  threshold, then the sources are considered to have
                                  different shape.
    """
    def __init__(self, bits, policy):
        SourceClassifier.__init__(self, bits, policy)
        self._numDifferentShape = 0
        self._numSimilarShape = 0
        self._shapeNormDiffThreshold = policy.getDouble("shapeNormDiffThreshold")
        assert self._shapeNormDiffThreshold > 0

    def _shapeNorm(self, ixx, iyy, ixy):
        """
        Computes norm of shape parameters
        """
        if ixx + iyy != 0.0:
            e1 = (ixx - iyy)/(ixx + iyy)
            e2 = 2*ixy/(ixx + iyy)
            return math.sqrt(e1*e1 + e2*e2)
        else: return 0.0

    def classify(self, *args):
        flag0 = args[0].getFlagClassification()
        flag1 = args[1].getFlagClassification()
        sn0 = self._shapeNorm(args[0].getIxx(), args[0].getIyy(), args[0].getIxy())
        sn1 = self._shapeNorm(args[1].getIxx(), args[1].getIyy(), args[1].getIxy())
        if abs(sn0 - sn1) > self._shapeNormDiffThreshold:
            args[0].setFlagClassification(self.setBit(flag0))
            args[1].setFlagClassification(self.setBit(flag1))
            self._numDifferentShape += 1
        else:
            args[0].setFlagClassification(self.clearBit(flag0))
            args[1].setFlagClassification(self.clearBit(flag1))
            self._numSimilarShape += 1

    def finish(self, log=None, clipboard=None):
        if log:
            LogRec(log, Log.INFO) << "ShapeDiffersInExposuresClassifier visit statistics" << \
                { "numDifferentShape": self._numDifferentShape, "numSimilarShape": self._numSimilarShape } << endr


class PositiveFluxExcursionClassifier(SourceClassifier):
    """
    Checks whether the flux excursion of a source is positive (for difference
    sources it can be negative).

    Sets two bits in the classification flags for each source:
    - bit 0 is set when a source has a positive flux excursion on the
      exposure it was measured on.
    - bit 1 is set when the sibling source (measured on the other
      exposure in the visit) has a positive flux excursion.

    Configuration parameters:
        "psfFluxThreshold": if the PSF flux of a source is above this threshold,
                            the source is considered to be a positive excursion
    """
    def __init__(self, bits, policy):
        SourceClassifier.__init__(self, bits, policy)
        self._numPositiveInBoth = 0
        self._numPositiveOnlyInFirst = 0
        self._numPositiveOnlyInSecond = 0
        self._numNegativeOrMissingInBoth = 0

        self._numMissingOrNegativeInFirst = 0
        self._psfFluxThreshold = policy.getDouble("psfFluxThreshold")
        assert self._psfFluxThreshold > 0

    def classify(self, *args):
        flag0 = args[0].getFlagClassification()
        flag1 = args[1].getFlagClassification()
        firstPos = args[0].getPsfFlux() > self._psfFluxThreshold
        secondPos = args[1].getPsfFlux() > self._psfFluxThreshold
        if firstPos:
            flag0 = self.setBit(flag0, 0)
            flag1 = self.setBit(flag1, 1)
            if secondPos: self._numPositiveInBoth += 1
            else:         self._numPositiveOnlyInFirst += 1
        else:
            flag0 = self.clearBit(flag0, 0)
            flag1 = self.clearBit(flag1, 1)
            if secondPos: self._numPositiveOnlyInSecond += 1
            else:         self._numNegativeOrMissingInBoth += 1
        if secondPos:
            flag0 = self.setBit(flag0, 1)
            flag1 = self.setBit(flag1, 0)
        else:
            flag0 = self.clearBit(flag0, 1)
            flag1 = self.clearBit(flag1, 0)
        args[0].setFlagClassification(flag0)
        args[1].setFlagClassification(flag1)

    def finish(self, log=None, clipboard=None):
        if log:
            rec = LogRec(log, Log.INFO)
            rec << "PositiveFluxExcursion visit statistics"
            rec << { "numPositiveInBoth":          self._numPositiveInBoth,
                     "numPositiveOnlyInFirst":     self._numPositiveOnlyInFirst,
                     "numPositiveOnlyInSecond":    self._numPositiveOnlyInSecond,
                     "numNegativeOrMissingInBoth": self._numNegativeOrMissingInBoth }
            rec << endr


class EllipticalAfterPSFDeconvolveClassifier(SourceClassifier):
    """
    Not available for DC3a. RHL says it is possible to estimate this from global
    PSF info and second moments, but it's not totally trivial, he's quite busy,
    and we really want to use local PSF information in the long run.
    """
    pass

