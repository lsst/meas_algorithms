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
"""This is a useful module for debugging measurements.

It relies on the presence of a catalog already produced by running
measurements in the regular context.  You provide the image and
catalog on the command-line, a list of source identifiers and the
measurement configuration in the config; the module reads the inputs,
subsets the catalog to contain only the sources of interest, and measures
those.  This reduces the frustration of waiting for image processing
and the measurements to run on many other sources, greatly increasing
debugging efficiency.
"""

import sys
from argparse import ArgumentParser, Namespace
from lsst.log import Log
from lsst.pex.config import Config, ConfigurableField, Field, ListField
from lsst.pipe.base import CmdLineTask, ConfigValueAction, ConfigFileAction, TaskRunner, Struct
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable

from lsst.meas.algorithms.measurement import SourceMeasurementTask


class MeasurementDebuggerConfig(Config):
    sourceId = ListField(dtype=int, default=[], doc="List of source identifiers to measure")
    outputName = Field(dtype=str, default="catalog.fits", doc="Output name for source catalog")
    measurement = ConfigurableField(target=SourceMeasurementTask, doc="Measurements")


class MeasurementDebuggerRunner(TaskRunner):
    """Provide the image and catalog names to the Task

    We provide a dummy dataRef only to avoid further overrides
    of this class.
    """
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        kwargs["image"] = parsedCmd.image
        kwargs["catalog"] = parsedCmd.catalog
        return [(Struct(dataId="<none>"), kwargs)]


class MeasurementDebuggerArgumentParser(ArgumentParser):
    """A stripped down version of the pipe_base ArgumentParser

    We don't care about the butler, just the config, and log.
    """

    def __init__(self, *args, **kwargs):
        super(MeasurementDebuggerArgumentParser, self).__init__(*args, **kwargs)
        self.add_argument("image", help="Name of image to measure")
        self.add_argument("catalog", help="Name of catalog to measure")
        self.add_argument("-c", "--config", nargs="*", action=ConfigValueAction,
                          help="config override(s), e.g. -c foo=newfoo bar.baz=3", metavar="NAME=VALUE")
        self.add_argument("-C", "--configfile", dest="configfile", nargs="*", action=ConfigFileAction,
                          help="config override file(s)")
        self.add_argument("--doraise", action="store_true",
                          help="raise an exception on error (else log a message and continue)?")
        self.add_argument("--logdest", help="logging destination")

    def parse_args(self, config, args=None, log=None, override=None):
        if args is None:
            args = sys.argv[1:]
        namespace = Namespace()
        namespace.config = config
        namespace.clobberConfig = False
        namespace.butler = None
        namespace.log = log if log is not None else Log.getDefaultLogger()
        namespace = super(MeasurementDebuggerArgumentParser, self).parse_args(args=args, namespace=namespace)
        del namespace.configfile
        return namespace


class MeasurementDebuggerTask(CmdLineTask):
    _DefaultName = "debugger"
    ConfigClass = MeasurementDebuggerConfig
    RunnerClass = MeasurementDebuggerRunner

    def __init__(self, schema=None, **kwargs):
        super(MeasurementDebuggerTask, self).__init__(**kwargs)
        if schema is None:
            schema = afwTable.SourceTable.makeMinimalSchema()
        self.schema = schema
        self.makeSubtask("measurement", schema=schema)

    @classmethod
    def _makeArgumentParser(cls):
        return MeasurementDebuggerArgumentParser()

    def run(self, dataRef, image, catalog):
        exp = self.readImage(image)
        sources = self.readSources(catalog)
        sources = self.subsetSources(sources)
        sources = self.mapSchemas(sources)
        self.measurement.measure(exp, sources)
        self.writeSources(sources)
        return Struct(exp=exp, sources=sources)

    def readImage(self, image):
        exp = afwImage.ExposureF(image)
        self.log.info("Read %dx%d image", exp.getWidth(), exp.getHeight())
        return exp

    def readSources(self, catalog):
        sources = afwTable.SourceCatalog.readFits(catalog)
        self.log.info("Read %d sources", len(sources))
        return sources

    def mapSchemas(self, sources):
        catalog = afwTable.SourceCatalog(self.schema)
        for ss in sources:
            new = catalog.addNew()
            new.setFootprint(ss.getFootprint())
            for name in self.schema.getNames():
                if name in ss.schema:
                    new.set(name, ss.get(name))
        return catalog

    def subsetSources(self, sources):
        """Return a subset of the input catalog

        The full catalog is used if the 'sourceId' list is empty.

        Parent sources (in the deblending sense) are also added to the
        subset so that they can be removed (via replaceWithNoise).
        """
        if not self.config.sourceId:
            return sources

        identifiers = set(self.config.sourceId)
        subset = afwTable.SourceCatalog(sources.table)
        while len(identifiers) > 0:
            ident = identifiers.pop()
            ss = sources.find(ident)
            if ss is None:
                raise RuntimeError("Unable to find id=%d in catalog" % ident)
            subset.append(ss)
            parent = ss.getParent()
            if parent:
                identifiers.add(parent)
        self.log.info("Subset to %d sources", len(subset))
        return subset

    def writeSources(self, sources):
        sources.writeFits(self.config.outputName)
        self.log.info("Wrote %s", self.config.outputName)

    def writeConfig(self, *args, **kwargs):
        pass

    def writeMetadata(self, *args, **kwargs):
        pass

    def writeSchemas(self, *args, **kwargs):
        pass
