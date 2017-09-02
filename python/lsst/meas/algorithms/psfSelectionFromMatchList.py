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
import numpy
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.afw.display.ds9 as ds9

args = [None, "MatchList", None]        # allow the user to probe for this signature


def selectPsfSources(exposure, matches, psfPolicy):
    """Get a list of suitable stars to construct a PSF."""

    import lsstDebug
    display = lsstDebug.Info(__name__).display
    displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells
    #
    # Unpack policy
    #
    kernelSize = psfPolicy.get("kernelSize")
    borderWidth = psfPolicy.get("borderWidth")
    sizePsfCellX = psfPolicy.get("sizeCellX")
    sizePsfCellY = psfPolicy.get("sizeCellY")
    #
    mi = exposure.getMaskedImage()

    if display and displayExposure:
        frame = 0
        ds9.mtv(mi, frame=frame, title="PSF candidates")

    psfCellSet = afwMath.SpatialCellSet(mi.getBBox(), sizePsfCellX, sizePsfCellY)
    psfStars = []

    for val in matches:
        ref, source = val[0:2]
        if not (ref.getFlagForDetection() & measAlg.Flags.STAR) or \
               (source.getFlagForDetection() & measAlg.Flags.BAD):
            continue

        try:
            cand = measAlg.makePsfCandidate(source, mi)
            #
            # The setXXX methods are class static, but it's convenient to call them on
            # an instance as we don't know Exposure's pixel type (and hence cand's exact type)
            if cand.getWidth() == 0:
                cand.setBorderWidth(borderWidth)
                cand.setWidth(kernelSize + 2*borderWidth)
                cand.setHeight(kernelSize + 2*borderWidth)

            im = cand.getMaskedImage().getImage()
            max = afwMath.makeStatistics(im, afwMath.MAX).getValue()
            if not numpy.isfinite(max):
                continue

            psfCellSet.insertCandidate(cand)

            if display and displayExposure:
                ds9.dot("+", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(),
                        size=4, frame=frame, ctype=ds9.CYAN)
                ds9.dot("o", source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0(),
                        size=4, frame=frame, ctype=ds9.CYAN)
        except Exception:
            continue

        source.setFlagForDetection(source.getFlagForDetection() | measAlg.Flags.STAR)
        psfStars += [source]

    return psfStars, psfCellSet
