
__all__ = ["SkyObjectsConfig", "SkyObjectsTask", "generateSkyObjects"]

from lsst.pex.config import Config, Field, ListField
import lsst.pipe.base as pipeBase

import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.math


class SkyObjectsConfig(Config):
    """Configuration for generating sky objects"""
    avoidMask = ListField(dtype=str, default=["DETECTED", "DETECTED_NEGATIVE", "BAD", "NO_DATA"],
                          doc="Avoid pixels masked with these mask planes")
    growMask = Field(dtype=int, default=0,
                     doc="Number of pixels to grow the masked pixels when adding sky objects")
    sourceRadius = Field(dtype=float, default=8, doc="Radius, in pixels, of sky objects")
    nSources = Field(dtype=int, default=100, doc="Try to add this many sky objects")
    nTrialSources = Field(dtype=int, default=None, optional=True,
                          doc="Maximum number of trial sky object positions\n"
                              "(default: nSkySources*nTrialSkySourcesMultiplier)")
    nTrialSourcesMultiplier = Field(dtype=int, default=5,
                                    doc="Set nTrialSkySources to\n"
                                        "    nSkySources*nTrialSkySourcesMultiplier\n"
                                        "if nTrialSkySources is None")


def generateSkyObjects(mask, seed, config):
    """Generate a list of Footprints of sky objects

    Sky objects don't overlap with other objects. This is determined
    through the provided `mask` (in which objects are typically flagged
    as `DETECTED`).

    The algorithm for determining sky objects is random trial and error:
    we try up to `nTrialSkySources` random positions to find `nSources`
    sky objects.

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

    rng = lsst.afw.math.Random(seed=seed)

    skyFootprints = []
    for _ in range(nTrialSkySources):
        if len(skyFootprints) == nSkySources:
            break

        x = int(rng.flat(xMin, xMax))
        y = int(rng.flat(yMin, yMax))
        spans = lsst.afw.geom.SpanSet.fromShape(int(skySourceRadius), offset=(x, y))
        if spans.overlaps(avoid):
            continue

        fp = lsst.afw.detection.Footprint(spans, mask.getBBox())
        fp.addPeak(x, y, 0)
        skyFootprints.append(fp)

    return skyFootprints


class SkyObjectsTask(pipeBase.Task):
    ConfigClass = SkyObjectsConfig

    def __init__(self, schema=None, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.schema = schema

        if schema is not None:
            self.skySourceKey = self.schema.addField(
                "sky_source",
                type="Flag",
                doc="Sky objects.")

    def run(self, mask, seed, schema=None):
        """Generate a list of Footprints of sky objects

        Sky objects don't overlap with other objects. This is determined
        through the provided `mask` (in which objects are typically flagged
        as `DETECTED`).

        The algorithm for determining sky objects is random trial and error:
        we try up to `nTrialSkySources` random positions to find `nSources`
        sky objects.

        Parameters
        ----------
        mask : `lsst.afw.image.Mask`
            Input mask plane, which identifies pixels to avoid for the sky
            objects.
        seed : `int`
            Random number generator seed.

        Returns
        -------
        skyFootprints : `list` of `lsst.afw.detection.Footprint`
            Footprints of sky objects. Each will have a peak at the center
            of the sky object.
        """
        skyFootprints = generateSkyObjects(mask, seed, self.config)
        self.log.info("Added %d of %d requested sky sources (%.0f%%)", len(skyFootprints),
                      self.config.nSources, 100*len(skyFootprints)/self.config.nSources)
        return skyFootprints
