from lsst.pex.logging import Log

import lsst.pex.policy as policy
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImg
import lsst.afw.math as afwMath
import lsst.pex.exceptions as pexExcept
import algorithmsLib as measAlg

def makeBackgroundSubtractedExposure(exposure, backgroundPolicy):
    algorithm = backgroundPolicy.get("algorithm");
    binsize = backgroundPolicy.get("binsize")
    
    maskedImage = exposure.getMaskedImage()
    #make a deep copy
    deepCopy = maskedImage.Factory(maskedImage, True)
    deepCopy.setXY0(maskedImage.getXY0())
      
    #
    # Subtract background
    #
    if algorithm == "NATURAL_SPLINE":
        bctrl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
    else:
        raise RuntimeError, "Unknown background algorithm: %s" % (algorithm)

    bctrl.setNxSample(int(deepCopy.getWidth()/binsize) + 1)
    bctrl.setNySample(int(deepCopy.getHeight()/binsize) + 1)
    background = afwMath.makeBackground(deepCopy.getImage(), bctrl)

    image = deepCopy.getImage() 
    image -= background.getImageF()
    del image

    return exposure.Factory(deepCopy, afwImg.Wcs(exposure.getWcs()))
    
def detectSources(exposure, psf, detectionPolicy): 
    """
    Description:
        Blur the background subtracted exposure by a psf, and then detect 
        sources. Detection parameters must be specified in the policy.
    """
    log = Log(Log.getDefaultLog(), "lsst.meas.algorithms.detectSources")

    if exposure == None:
        log.log(Log.FATAL, "Cannot perform detection - missing an exposure")
        raise RuntimeException("No exposure for detection")
    if psf == None:
        log.log(Log.FATAL, "Cannot perform detection - missing a psf")
        raise RuntimeException("No psf for detection")
    if detectionPolicy == None:
        log.log(Log.FATAL, "Cannot perform detection - missing a policy")
        raise RuntimeException("No policy for detection")
    
    #
    # Unpack variables
    #
    minPix = detectionPolicy.get("minPixels")
    thresholdValue = detectionPolicy.get("thresholdValue")
    thresholdType = detectionPolicy.getString("thresholdType")
    thresholdPolarity = detectionPolicy.getString("thresholdPolarity")
    maskedImage = exposure.getMaskedImage()
    convolvedImage = maskedImage.Factory(maskedImage.getDimensions())
    convolvedImage.setXY0(maskedImage.getXY0())

    # 
    # Smooth the Image
    #
    psf.convolve(convolvedImage, maskedImage, 
            convolvedImage.getMask().getMaskPlane("EDGE"))
   




    #
    # Only search psf-smooth part of frame
    #
    llc = afwImg.PointI(psf.getKernel().getWidth()/2, psf.getKernel().getHeight()/2)
    urc = afwImg.PointI(convolvedImage.getWidth() - 1, convolvedImage.getHeight() - 1)
    urc -= llc
    detectionBbox = afwImg.BBox(llc, urc)
    imageRegion = afwImg.BBox(
            afwImg.PointI(maskedImage.getX0(), maskedImage.getY0()),
            maskedImage.getWidth(), maskedImage.getHeight())
    middle = convolvedImage.Factory(convolvedImage, detectionBbox)
    
    grow, isotropic = 1, False
    dsPositive = None
    dsNegative = None 
    if thresholdPolarity == "negative" or thresholdPolarity == "both":            
        #detect negative sources
        log.log(Log.DEBUG, "Do Negative Detection")
        #create a Threshold for negative detections
        threshold = afwDet.createThreshold(thresholdValue, thresholdType, False)
        #perfom detection
        dsNegative = afwDet.makeDetectionSet(middle, threshold, 
                "DETECTED_NEGATIVE", minPix)
        # ds only searched the middle but it belongs to the entire MaskedImage
        dsNegative.setRegion(imageRegion)
        dsNegative = afwDet.DetectionSetF(dsNegative, grow, isotropic)
        dsNegative.setMask(maskedImage.getMask(), "DETECTED_NEGATIVE")
    elif thresholdPolarity != "negative":
        log.log(Log.DEBUG, "Do Positive Detection")
        #create a Threshold for negative detections
        threshold = afwDet.createThreshold(thresholdValue, thresholdType, True)
        dsPositive = afwDet.makeDetectionSet(middle, threshold, 
                "DETECTED", minPix)
        # ds only searched the middle but it belongs to the entire MaskedImage
        dsPositive.setRegion(imageRegion)
        # grow the detections into the edge by at least one pixel
        #
        dsPositive = afwDet.DetectionSetF(dsPositive, grow, isotropic)
        dsPositive.setMask(maskedImage.getMask(), "DETECTED")

    #
    # clean up
    #
    del middle

    return dsPositive, dsNegative
