# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

import lsst.pex.config as pexConfig
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase

from . import algorithmsLib
from .algorithmRegistry import *

__all__ = "SourceSlotConfig", "SourceMeasurementConfig", "SourceMeasurementTask"

class SourceSlotConfig(pexConf.Config):

    centroid = pexConf.Field(dtype=str, default="centroid.sdss", optional=True,
                             doc="the name of the centroiding algorithm used to set source x,y")
    shape = pexConf.Field(dtype=str, default="shape.sdss", optional=True,
                          doc="the name of the algorithm used to set source moments parameters")
    apFlux = pexConf.Field(dtype=str, default="flux.sinc", optional=True,
                           doc="the name of the algorithm used to set the source aperture flux slot")
    modelFlux = pexConf.Field(dtype=str, default="flux.gaussian", optional=True,
                           doc="the name of the algorithm used to set the source model flux slot")
    psfFlux = pexConf.Field(dtype=str, default="flux.psf", optional=True,
                            doc="the name of the algorithm used to set the source psf flux slot")
    instFlux = pexConf.Field(dtype=str, default="flux.gaussian", optional=True,
                             doc="the name of the algorithm used to set the source inst flux slot")

    def setupTable(self, table):
        """Convenience method to setup a table's slots according to the config definition.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        if self.centroid is not None: table.defineCentroid(self.centroid)
        if self.shape is not None: table.defineShape(self.shape)
        if self.apFlux is not None: table.defineApFlux(self.apFlux)
        if self.modelFlux is not None: table.defineModelFlux(self.modelFlux)
        if self.psfFlux is not None: table.definePsfFlux(self.psfFlux)
        if self.instFlux is not None: table.defineInstFlux(self.instFlux)

class SourceMeasurementConfig(pexConf.Config):
    """
    Configuration for SourceMeasurementTask.
    A configured instance of MeasureSources can be created using the
    makeMeasureSources method.
    """

    slots = pexConf.ConfigField(
        dtype = SourceSlotConfig,
        doc="Mapping from algorithms to special aliases in Source.\n"
        )

    algorithms = AlgorithmRegistry.all.makeField(
        multi=True,
        default=["flags.pixel",
                 "centroid.gaussian", "centroid.naive",
                 "shape.sdss",
                 "flux.gaussian", "flux.naive", "flux.psf", "flux.sinc",
                 "classification.extendedness",
                 ],
        doc="Configuration and selection of measurement algorithms."
        )
    
    centroider = AlgorithmRegistry.filter(CentroidConfig).makeField(
        multi=False, default="centroid.sdss", optional=True,
        doc="Configuration for the initial centroid algorithm used to\n"\
            "feed center points to other algorithms.\n\n"\
            "Note that this is in addition to the centroider listed in\n"\
            "the 'algorithms' field; the same name should not appear in\n"\
            "both.\n\n"\
            "This field DOES NOT set which field name will be used to define\n"\
            "the alias for source.getX(), source.getY(), etc.\n"
        )

    def __init__(self):
        pexConf.Config.__init__(self)
        self.slots.centroid = self.centroider.name
        self.slots.shape = "shape.sdss"
        self.slots.psfFlux = "flux.psf"
        self.slots.apFlux = "flux.naive"
        self.slots.modelFlux = "flux.gaussian"
        self.slots.instFlux = "flux.gaussian"

    def validate(self):
        pexConf.Config.validate(self)
        if self.centroider.name in self.algorithms.names:
            raise ValueError("The algorithm in the 'centroider' field must not also appear in the "\
                                 "'algorithms' field.")
        if self.slots.centroid is not None and (self.slots.centroid not in self.algorithms.names
                                                and self.slots.centroid != self.centroider.name):
            raise ValueError("source centroid slot algorithm '%s' is not being run." % self.slots.astrom)
        if self.slots.shape is not None and self.slots.shape not in self.algorithms.names:
            raise ValueError("source shape slot algorithm '%s' is not being run." % self.slots.shape)
        for slot in (self.slots.psfFlux, self.slots.apFlux, self.slots.modelFlux, self.slots.instFlux):
            if slot is not None and slot not in self.algorithms.names:
                raise ValueError("source flux slot algorithm '%s' is not being run." % slot)

    def makeMeasureSources(self, schema):
        """ Convenience method to make a MeasureSources instance and
        fill it with the configured algorithms.

        This is defined in the Config class to support use in unit tests without needing
        to construct a Task object.
        """
        builder = algorithmsLib.MeasureSourcesBuilder()
        if self.centroider is not None:
            builder.setCentroider(self.centroider.apply())
        builder.addAlgorithms(self.algorithms.apply())
        return builder.build(schema)

class SourceMeasurementTask(pipeBase.Task):
    """Measure the properties of sources on a single exposure.

    This task has no return value; it only modifies the SourceVector in-place.
    """
    ConfigClass = SourceMeasurementConfig

    def __init__(self, config, schema, **kwds):
        """Create the task, adding necessary fields to the given schema.
        """
        pipeBase.Task.__init__(self, config=config, **kwds)
        self.measurer = config.makeMeasureSources(schema)

    @pipeBase.timeMethod
    def run(self, exposure, sources):
        """Measure sources on an exposure.

        @param exposure Exposure to process
        @param sources  SourceVector containing sources detected on this exposure.
        @return None
        """
        assert exposure, "No exposure provided"
        self.config.slots.setupTable(sources.table)
        for record in sources:
            self.measurer.apply(record, exposure)
