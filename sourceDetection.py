from lsst.pex.logging import Log

import lsst.pex.policy as policy
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImg
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg


def makePsf(psfPolicy):
    params = []        
    params.append(psfPolicy.getString("algorithm"))
    params.append(psfPolicy.getInt("width"))
    params.append(psfPolicy.getInt("height"))
    if psfPolicy.exists("parameter"):
        params += psfPolicy.getDoubleArray("parameter")
        
    return measAlg.createPSF(*params)

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
    bctrl = afwMath.BackgroundControl(backgroundPolicy.get("algorithm"))
    binsize = backgroundPolicy.get("binsize")

    # Set background control parameters
    bctrl.setNxSample(image.getWidth()/binsize + 1)
    bctrl.setNySample(image.getHeight()/binsize + 1)

    #return a background object
    return afwMath.makeBackground(image, bctrl)
    
def estimateBackground(exposure, backgroundPolicy, subtract=True):
    """
    Estimate exposure's background using parameters in backgroundPolicy.  
    If subtract is true, make a copy of the exposure and subtract the background.  
    Return background, backgroundSubtractedExposure
    """
    maskedImage = exposure.getMaskedImage()
    image = maskedImage.getImage()    
    background = getBackground(image, backgroundPolicy)


    if not background:
        raise RuntimeError, "Unable to estimate background for exposure"
    
    if not subtract:
        return background, None

    bbox = afwImg.BBox(
        afwImg.PointI(0,0), 
        maskedImage.getWidth(), 
        maskedImage.getHeight()
    )   
    backgroundSubtractedExposure = exposure.Factory(exposure, bbox, True)
    backgroundSubtractedExposure.getMaskedImage().setXY0(maskedImage.getXY0())
    copyImage = backgroundSubtractedExposure.getMaskedImage().getImage()
    copyImage -= background.getImageF()
    del copyImage
    del maskedImage
    del image

    return background, backgroundSubtractedExposure

def detectSources(exposure, psf, detectionPolicy):
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
    region = afwImg.BBox(
        afwImg.PointI(maskedImage.getX0(), maskedImage.getY0()),
        maskedImage.getWidth(), 
        maskedImage.getHeight()
    )
    
    if psf is None:
        middle = maskedImage.Factory(maskedImage)
    else:
        convolvedImage = maskedImage.Factory(maskedImage.getDimensions())
        convolvedImage.setXY0(maskedImage.getXY0())

        # 
        # Smooth the Image
        #
        psf.convolve(convolvedImage, 
            maskedImage, 
            convolvedImage.getMask().getMaskPlane("EDGE")
        )
        #
        # Only search psf-smooth part of frame
        #
        llc = afwImg.PointI(
            psf.getKernel().getWidth()/2, 
            psf.getKernel().getHeight()/2
        )
        urc = afwImg.PointI(
            convolvedImage.getWidth() - 1,
            convolvedImage.getHeight() - 1
        )
        urc -= llc
        bbox = afwImg.BBox(llc, urc)    
        middle = convolvedImage.Factory(convolvedImage, bbox)
	del convolvedImage

    dsPositive = None
    if thresholdPolarity != "negative":
        #do positive detections        
        threshold = afwDet.createThreshold(
            thresholdValue,
            thresholdType,
            True
        )
        dsPositive = afwDet.makeFootprintSet(
            middle,
            threshold,
            "DETECTED",
            minPixels
        )
        #set detection region to be the entire region
        dsPositive.setRegion(region)
        
        # We want to grow the detections into the edge by at least one pixel so 
        # that it sees the EDGE bit
        if nGrow > 0:
            dsPositive = afwDet.FootprintSetF(dsPositive, nGrow, False)
            dsPositive.setMask(maskedImage.getMask(), "DETECTED")

    dsNegative = None
    if thresholdPolarity != "positive":
        #do negative detections
        threshold = afwDet.createThreshold(
            thresholdValue,
            thresholdType,
            False
        )
        #detect negative sources
        dsNegative = afwDet.makeFootprintSet(
            middle,
            threshold,
            "DETECTED_NEGATIVE",
            minPixels
        )
        #set detection region to be the entire region
        dsNegative.setRegion(region)
        
        # We want to grow the detections into the edge by at least one pixel so 
        # that it sees the EDGE bit        
        if nGrow > 0:
            dsNegative = afwDet.FootprintSetF(dsNegative, nGrow, False)
            dsNegative.setMask(maskedImage.getMask(), "DETECTED_NEGATIVE")

    #
    # clean up
    #
    del middle    

    return dsPositive, dsNegative


