#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

import numpy

import lsst.pex.config
import lsst.pex.exceptions
import lsst.afw.image
import lsst.pipe.base

import math

from . import algorithmsLib

__all__ = ("MeasureApCorrConfig", "MeasureApCorrTask")

class KeyTuple(object):

    __slots__ = ("flux", "err", "flag", "used",)

    def __init__(self, name, schema):
        self.flux = schema.find(name).key
        self.err = schema.find(name + ".err").key
        self.flag = schema.find(name + ".flags").key
        self.used = schema.addField("apcorr." + name + ".used", type="Flag",
                                    doc="set if source was used in measuring aperture correction")

class MeasureApCorrConfig(lsst.pex.config.Config):
    reference = lsst.pex.config.Field(
        dtype=str, default="flux.naive",
        doc="Name of the flux field other measurements should be corrected to match"
    )
    inputFilterFlag = lsst.pex.config.Field(
        dtype=str, default="calib.psf.used",
        doc=("Name of a flag field that indicates that a source should be used to constrain the"
             " aperture corrections")
    )
    minDegreesOfFreedom = lsst.pex.config.Field(
        dtype=int, default=1, check=lambda x: x > 0,
        doc=("Minimum number of degrees of freedom (# of valid data points - # of parameters);"
             " if this is exceeded, the order of the fit is decreased (in both dimensions), and"
             " if we can't decrease it enough, we'll raise ValueError.")
        )
    fit = lsst.pex.config.ConfigField(
        dtype=lsst.afw.math.ChebyshevBoundedFieldConfig,
        doc="Configuration used in fitting the aperture correction fields"
    )
    numIter = lsst.pex.config.Field(
        dtype=int, default=4,
        doc="Number of iterations for sigma clipping"
    )
    numSigmaClip = lsst.pex.config.Field(
        dtype=float, default=3.0,
        doc="Number of standard devisations to clip at"
    )

class MeasureApCorrTask(lsst.pipe.base.Task):

    ConfigClass = MeasureApCorrConfig
    _DefaultName = "measureApCorr"

    def __init__(self, schema, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.reference = KeyTuple(self.config.reference, schema)
        self.toCorrect = {}
        for name in algorithmsLib.getApCorrRegistry():
            try:
                self.toCorrect[name] = KeyTuple(name, schema)
            except KeyError:
                # if a field in the registry is missing, just ignore it.
                pass
        self.inputFilterFlag = schema.find(self.config.inputFilterFlag).key

    def run(self, bbox, catalog):

        # First, create a subset of the catalog that contains only objects with inputFilterFlag set
        # and non-flagged reference fluxes.
        subset1 = [record for record in catalog
                   if record.get(self.inputFilterFlag) and not record.get(self.reference.flag)]

        apCorrMap = lsst.afw.image.ApCorrMap()

        # Outer loop over the fields we want to correct
        for name, keys in self.toCorrect.iteritems():

            # Create a more restricted subset with only the objects where the to-be-correct flux
            # is not flagged.
            subset2 = [record for record in subset1 if not record.get(keys.flag)]

            # Check that we have enough data points that we have at least the minimum of degrees of
            # freedom specified in the config.
            if len(subset2) - 1 < self.config.minDegreesOfFreedom:
                self.log.warn("Only %d sources for calculation of aperture correction for '%s'; "
                              "setting to 1.0" % (len(subset2), name,))
                apCorrMap[name] = lsst.afw.math.ChebyshevBoundedField(bbox, numpy.ones((1,1), dtype=float))
                apCorrMap[name + ".err"] = \
                    lsst.afw.math.ChebyshevBoundedField(bbox, numpy.zeros((1,1), dtype=float))
                continue

            # If we don't have enough data points to constrain the fit, reduce the order until we do
            ctrl = self.config.fit.makeControl()
            while len(subset2) - ctrl.computeSize() < self.config.minDegreesOfFreedom:
                if ctrl.orderX > 0:
                    ctrl.orderX -= 1
                if ctrl.orderY > 0:
                    ctrl.orderY -= 1

            # Fill numpy arrays with positions and the ratio of the reference flux to the to-correct flux
            x = numpy.zeros(len(subset2), dtype=float)
            y = numpy.zeros(len(subset2), dtype=float)
            apCorrData = numpy.zeros(len(subset2), dtype=float)
            indices = numpy.arange(len(subset2), dtype=int)
            for n, record in enumerate(subset2):
                x[n] = record.getX()
                y[n] = record.getY()
                apCorrData[n] = record.get(self.reference.flux)/record.get(keys.flux)

            for _i in range(self.config.numIter):

                # Do the fit, save it in the output map
                apCorrField = lsst.afw.math.ChebyshevBoundedField.fit(bbox, x, y, apCorrData, ctrl)

                # Compute errors empirically, using the RMS difference between the true reference flux and the
                # corrected to-be-corrected flux.
                apCorrDiffs = apCorrField.evaluate(x, y)
                apCorrDiffs -= apCorrData
                apCorrErr = numpy.mean(apCorrDiffs**2)**0.5

                # Clip bad data points
                apCorrDiffLim = self.config.numSigmaClip * apCorrErr
                keep = numpy.fabs(apCorrDiffs) < apCorrDiffLim
                x = x[keep]
                y = y[keep]
                apCorrData = apCorrData[keep]
                indices = indices[keep]

            # Final fit after clipping
            apCorrField = lsst.afw.math.ChebyshevBoundedField.fit(bbox, x, y, apCorrData, ctrl)

            self.log.info("Aperture correction for %s: RMS %f from %d" %
                          (name, numpy.mean((apCorrField.evaluate(x, y) - apCorrData)**2)**0.5, len(indices)))

            # Save the result in the output map
            # The error is constant spatially (we could imagine being
            # more clever, but we're not yet sure if it's worth the effort).
            # We save the errors as a 0th-order ChebyshevBoundedField
            apCorrMap[name] = apCorrField
            apCorrErrCoefficients = numpy.array([[apCorrErr]], dtype=float)
            apCorrMap[name + ".err"] = lsst.afw.math.ChebyshevBoundedField(bbox, apCorrErrCoefficients)

            # Record which sources were used
            for i in indices:
                subset2[i].set(keys.used, True)

        return apCorrMap
