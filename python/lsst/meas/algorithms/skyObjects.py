# This file is part of meas_algorithms.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["SkyObjectsConfig", "SkyObjectsTask", "generateSkyObjects"]

from scipy.stats import qmc

from lsst.pex.config import Config, Field, ListField
from lsst.pipe.base import Task

import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.math


class SkyObjectsConfig(Config):
    """Configuration for generating sky objects"""
    avoidMask = ListField(
        dtype=str,
        default=["DETECTED", "DETECTED_NEGATIVE", "BAD", "NO_DATA"],
        doc="Avoid pixels masked with these mask planes."
    )
    growMask = Field(
        dtype=int,
        default=0,
        doc="Number of pixels to grow the masked pixels when adding sky sources."
    )
    sourceRadius = Field(
        dtype=float,
        default=8,
        doc="Radius, in pixels, of sky sources."
    )
    nSources = Field(
        dtype=int,
        default=100,
        doc="Try to add this many sky sources."
    )
    nTrialSources = Field(
        dtype=int,
        default=None,
        optional=True,
        doc="Maximum number of trial sky object positions "
            "(default: nSkySources*nTrialSkySourcesMultiplier)."
    )
    nTrialSourcesMultiplier = Field(
        dtype=int,
        default=5,
        doc="Set nTrialSkySources to nSkySources*nTrialSkySourcesMultiplier "
            "if nTrialSkySources is None."
    )


def generateSkyObjects(mask, seed, config):
    """Generate a list of Footprints of sky objects

    Sky objects don't overlap with other objects. This is determined
    through the provided `mask` (in which objects are typically flagged
    as `DETECTED`).

    Sky objects are positioned using a quasi-random Halton sequence number
    generator. This is a deterministic sequence that mimics a random trial and
    error approach whilst acting to minimize clustering of points for a given
    field of view. Up to `nTrialSources` points are generated, returning the
    first `nSources` that do not overlap with the mask.

    Parameters
    ----------
    mask : `lsst.afw.image.Mask`
        Input mask plane, which identifies pixels to avoid for the sky
        objects.
    seed : `int`
        Random number generator seed.
    config : `SkyObjectsConfig`
        Configuration for finding sky objects.

    Returns
    -------
    skyFootprints : `list` of `lsst.afw.detection.Footprint`
        Footprints of sky objects. Each will have a peak at the center
        of the sky object.
    """
    if config.nSources <= 0:
        return []

    skySourceRadius = config.sourceRadius
    nSkySources = config.nSources
    nTrialSkySources = config.nTrialSources
    if nTrialSkySources is None:
        nTrialSkySources = config.nTrialSourcesMultiplier*nSkySources

    box = mask.getBBox()
    box.grow(-(int(skySourceRadius) + 1))  # Avoid objects partially off the image
    xMin, yMin = box.getMin()
    xMax, yMax = box.getMax()

    avoid = lsst.afw.geom.SpanSet.fromMask(mask, mask.getPlaneBitMask(config.avoidMask))
    if config.growMask > 0:
        avoid = avoid.dilated(config.growMask)

    sampler = qmc.Halton(d=2, seed=seed).random(nTrialSkySources)
    sample = qmc.scale(sampler, [xMin, yMin], [xMax, yMax])

    skyFootprints = []
    for x, y in zip(sample[:, 0].astype(int), sample[:, 1].astype(int)):
        if len(skyFootprints) == nSkySources:
            break

        spans = lsst.afw.geom.SpanSet.fromShape(int(skySourceRadius), offset=(x, y))
        if spans.overlaps(avoid):
            continue

        fp = lsst.afw.detection.Footprint(spans, mask.getBBox())
        fp.addPeak(x, y, 0)
        skyFootprints.append(fp)

        # Add doubled-in-size sky object spanSet to the avoid mask.
        avoid = avoid.union(spans.dilated(int(skySourceRadius)))

    return skyFootprints


class SkyObjectsTask(Task):
    """Generate a list of Footprints of sky sources/objects (regions on the
    sky that do not otherwise have detections).

    Parameters
    ----------
    schema : `lsst.afw.table.Schema`
        Schema used to create the output `~lsst.afw.table.SourceCatalog`,
        updated with fields that will be written by this task.

    """
    ConfigClass = SkyObjectsConfig

    def __init__(self, schema=None, **kwargs):
        super().__init__(**kwargs)
        if schema is not None:
            self.skySourceKey = schema.addField("sky_source", type="Flag",
                                                doc="Region on image with no detections.")
        else:
            self.skySourceKey = None

    def run(self, mask, seed, catalog=None):
        """Generate a list of Footprints of sky sources/objects.

        Sky objects don't overlap with other objects. This is determined
        through the provided `mask` (in which objects are typically flagged
        as `DETECTED`).

        Sky objects are positioned using a quasi-random Halton sequence
        number generator. This is a deterministic sequence that mimics a random
        trial and error approach whilst acting to minimize clustering of points
        for a given field of view. Up to `nTrialSources` points are generated,
        returning the first `nSources` that do not overlap with the mask.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Input mask plane, which identifies pixels to avoid for the sky
            objects.
        seed : `int`
            Random number generator seed.
        catalog : `lsst.afw.table.SourceCatalog`, optional
            Catalog to add detected footprints to; modified in-place if any
            sky source/object footprints are created.

        Returns
        -------
        skyFootprints : `list` of `lsst.afw.detection.Footprint`
            Footprints of sky objects. Each will have a peak at the center
            of the sky object.
        """
        skyFootprints = generateSkyObjects(mask, seed, self.config)
        self.log.info("Added %d of %d requested sky sources (%.0f%%)", len(skyFootprints),
                      self.config.nSources, 100*len(skyFootprints)/self.config.nSources)
        self.metadata["sky_footprint_count"] = len(skyFootprints)

        if skyFootprints and self.skySourceKey is not None and catalog is not None:
            for footprint in skyFootprints:
                record = catalog.addNew()
                record.setFootprint(footprint)
                record.set(self.skySourceKey, True)

        return skyFootprints
