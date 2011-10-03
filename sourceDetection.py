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

import lsstDebug
from lsst.pex.logging import Log

import lsst.daf.persistence as dafPersist
import lsst.pex.policy as policy
import lsst.afw.detection as afwDet
import lsst.afw.display.ds9 as ds9
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg

def makePsf(psfPolicy):
    params = []        
    params.append(psfPolicy.getString("algorithm"))
    params.append(psfPolicy.getInt("width"))
    params.append(psfPolicy.getInt("height"))
    if psfPolicy.exists("parameter"):
        params += psfPolicy.getDoubleArray("parameter")
        
    return afwDet.createPsf(*params)

def addExposures(exposureList):
    """
    Add a set of exposures together. 
    Assumes that all exposures in set have the same dimensions
    """
    exposure0 = exposureList[0]
    image0 = exposure0.getMaskedImage()

    addedImage = image0.Factory(image0, True)
    addedImage.setXY0(image0.getXY0())

    for exposure in exposureList[1:]:
        image = exposure.getMaskedImage()
        addedImage += image

    addedExposure = exposure0.Factory(addedImage, exposure0.getWcs())
    return addedExposure

def getBackground(image, backgroundPolicy):
    """
    Make a new Exposure which is exposure - background
    """
    binsize = backgroundPolicy.get("binsize")
    statProp = backgroundPolicy.get("statisticsproperty")
    undersamplestyle = backgroundPolicy.get("undersamplestyle")

    nx = image.getWidth()/binsize + 1
    ny = image.getHeight()/binsize + 1

    sctrl = afwMath.StatisticsControl()
    try:
        sctrl.setAndMask(image.getMask().getPlaneBitMask("DETECTED"))
    except AttributeError:
        pass

    if False:                           # doesn't work as there's no bctrl.setStatisticsControl
        bctrl = afwMath.BackgroundControl(backgroundPolicy.get("algorithm"))

        # Set background control parameters
        bctrl.setNxSample(nx)
        bctrl.setNySample(ny)
        bctrl.setUndersampleStyle(undersamplestyle)
        bctrl.setStatisticsProperty(statProp)
        bctrl.setStatisticsControl(sctrl)
    else:
        bctrl = afwMath.BackgroundControl(backgroundPolicy.get("algorithm"), nx, ny, undersamplestyle,
                                          sctrl, statProp)

    #return a background object
    return afwMath.makeBackground(image, bctrl)
    
def estimateBackground(exposure, backgroundPolicy, subtract=True):
    """
    Estimate exposure's background using parameters in backgroundPolicy.  
    If subtract is true, make a copy of the exposure and subtract the background.  
    Return background, backgroundSubtractedExposure
    """
    displayEstimateBackground = lsstDebug.Info(__name__).displayEstimateBackground

    maskedImage = exposure.getMaskedImage()

    sctrl = afwMath.StatisticsControl()
    sctrl.setAndMask(maskedImage.getMask().getPlaneBitMask("DETECTED")) # ignore detected pixels

    background = getBackground(maskedImage, backgroundPolicy)

    if not background:
        raise RuntimeError, "Unable to estimate background for exposure"
    
    if displayEstimateBackground > 1:
        ds9.mtv(background.getImageF(), title="background", frame=1)

    if not subtract:
        return background, None

    bbox = maskedImage.getBBox(afwImage.PARENT)
    backgroundSubtractedExposure = exposure.Factory(exposure, bbox, afwImage.PARENT, True)
    copyImage = backgroundSubtractedExposure.getMaskedImage().getImage()
    copyImage -= background.getImageF()

    if displayEstimateBackground:
        ds9.mtv(backgroundSubtractedExposure, title="subtracted")

    return background, backgroundSubtractedExposure

def setEdgeBits(maskedImage, goodBBox, edgeBitmask):
    """Set the edgeBitmask bits for all of maskedImage outside goodBBox"""
    msk = maskedImage.getMask()

    mx0, my0 = maskedImage.getXY0()
    for x0, y0, w, h in ([0, 0,
                          msk.getWidth(), goodBBox.getBeginY() - my0],
                         [0, goodBBox.getEndY() - my0, msk.getWidth(),
                          maskedImage.getHeight() - (goodBBox.getEndY() - my0)],
                         [0, 0,
                          goodBBox.getBeginX() - mx0, msk.getHeight()],
                         [goodBBox.getEndX() - mx0, 0,
                          maskedImage.getWidth() - (goodBBox.getEndX() - mx0), msk.getHeight()],
                         ):
        edgeMask = msk.Factory(msk, afwGeom.BoxI(afwGeom.PointI(x0, y0),
                                                 afwGeom.ExtentI(w, h)), afwImage.LOCAL)
        edgeMask |= edgeBitmask

def thresholdImage(image, thresholdValue, thresholdType, thresholdParity, extraThreshold, minPixels):
    """Threshold the convolved image, returning a FootprintSet.
    Helper function for detectSources().

    @param image The (optionally convolved) MaskedImage to threshold
    @param thresholdValue Value for the threshold
    @param thresholdType Type of threshold
    @param thresholdParity Parity of threshold
    @param extraThreshold Threshold multiplier to apply (faint sources discarded, footprints unaffected)
    @param minPixels Minimum number of pixels in footprint
    @return FootprintSet
    """
    parity = False if thresholdParity == "negative" else True
    threshold = afwDet.createThreshold(thresholdValue, thresholdType, parity)
    footprints = afwDet.makeFootprintSet(image, threshold, "DETECTED", minPixels)

    if extraThreshold > 1.0:
        thresh = threshold.getValue(image) * extraThreshold
        sifted = afwDet.FootprintContainerT()
        for f in footprints.getFootprints():
            boxes = afwDet.footprintToBBoxList(f)
            peak = 0
            for b in boxes:
                subImage = image.Factory(image, b, afwImage.PARENT, False)
                value = afwMath.makeStatistics(subImage, afwMath.MAX if parity else afwMath.MIN).getValue()
                if (parity and value > peak) or (not parity and value < peak):
                    peak = value
            if (parity and peak > thresh) or (not parity and peak < thresh):
                sifted.push_back(f)
        footprints.setFootprints(sifted)
    
    return footprints



def detectSources(exposure, psf, detectionPolicy, extraThreshold=1.0):
    try:
        import lsstDebug

        display = lsstDebug.Info(__name__).display
    except ImportError, e:
        try:
            display
        except NameError:
            display = False

    minPixels = detectionPolicy.get("minPixels")
    
    thresholdValue = detectionPolicy.get("thresholdValue")
    thresholdType = detectionPolicy.get("thresholdType")
    thresholdPolarity = detectionPolicy.get("thresholdPolarity")
    nGrow = detectionPolicy.get("nGrow")

    if exposure is None:
        raise RuntimeException("No exposure for detection")
       
    #
    # Unpack variables
    #
    maskedImage = exposure.getMaskedImage()
    region = maskedImage.getBBox(afwImage.PARENT)

    mask = maskedImage.getMask()
    mask &= ~(mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE"))
    del mask

    if psf is None:
        convolvedImage = maskedImage.Factory(maskedImage)
        middle = convolvedImage
    else:
        # We may have a proxy;  if so instantiate it
        if isinstance(psf, dafPersist.readProxy.ReadProxy):
            psf = psf.__subject__

        ##########
        # use a separable psf for convolution ... the psf width for the center of the image will do
        
        xCen = maskedImage.getX0() + maskedImage.getWidth()/2
        yCen = maskedImage.getY0() + maskedImage.getHeight()/2

        # measure the 'sigma' of the psf we were given
        psfAttrib = measAlg.PsfAttributes(psf, xCen, yCen)
        sigma = psfAttrib.computeGaussianWidth()

        # make a SingleGaussian (separable) kernel with the 'sigma'
        gaussFunc = afwMath.GaussianFunction1D(sigma)
        gaussKernel = afwMath.SeparableKernel(psf.getKernel().getWidth(), psf.getKernel().getHeight(),
                                              gaussFunc, gaussFunc)
        
        convolvedImage = maskedImage.Factory(maskedImage.getDimensions())
        convolvedImage.setXY0(maskedImage.getXY0())

        afwMath.convolve(convolvedImage, maskedImage, gaussKernel, afwMath.ConvolutionControl())
        #
        # Only search psf-smooth part of frame
        #
        goodBBox = gaussKernel.shrinkBBox(convolvedImage.getBBox(afwImage.PARENT))
        middle = convolvedImage.Factory(convolvedImage, goodBBox, afwImage.PARENT, False)
        #
        # Mark the parts of the image outside goodBBox as EDGE
        #
        setEdgeBits(maskedImage, goodBBox, maskedImage.getMask().getPlaneBitMask("EDGE"))

    dsPositive, dsNegative = None, None
    if thresholdPolarity != "negative":
        dsPositive = thresholdImage(middle, thresholdValue, thresholdType,
                                    "positive", extraThreshold, minPixels)
    if thresholdPolarity != "positive":
        dsNegative = thresholdImage(middle, thresholdValue, thresholdType,
                                    "negative", extraThreshold, minPixels)

    for footprints in (dsPositive, dsNegative):
        if footprints is None:
            continue
        footprints.setRegion(region)
        if nGrow > 0:
            footprints = afwDet.FootprintSetF(footprints, nGrow, False)
        footprints.setMask(maskedImage.getMask(), "DETECTED")

    if display:
        ds9.mtv(exposure, frame=0, title="detection")

        if convolvedImage and display and display > 1:
            ds9.mtv(convolvedImage, frame=1, title="PSF smoothed")

        if middle and display and display > 1:
            ds9.mtv(middle, frame=2, title="middle")

    return dsPositive, dsNegative


