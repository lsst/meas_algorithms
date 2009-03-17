import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.policy as policy
import lsst.afw.detection as afwDetection
import algorithmsLib as algorithms

class PsfShapeHistogram(object):
    """A class to represent a histogram of (Ixx, Iyy)"""

    def __init__(self, xSize=20, ySize=20, xMax=15, yMax=15):
        self._xSize, self._ySize = xSize, ySize 
        self._xMax, self._yMax = xMax, yMax
        self._psfImage = afwImage.ImageF(self._xSize, self._ySize)
        self._psfImage.set(0)

    def getImage(self):
        return self._psfImage

    def insert(self, source):
        """Insert source into the histogram."""
        i = int(source.getIxx()*self._xSize/self._xMax + 0.5)
        j = int(source.getIyy()*self._ySize/self._yMax + 0.5)
        if i in range(0, self._xSize) and j in range(0, self._ySize):
            if i != 0 or j != 0:
                self._psfImage.set(i, j, self._psfImage.get(i, j) + 1)
    
    def peakToIxx(self, peakX, peakY):
        """
        Given a peak position in self._psfImage, 
        return the corresponding (Ixx, Iyy)
        """
        xx = peakX*self._xMax/self._xSize
        yy = peakY*self._yMax/self._ySize
        return xx, yy

    def getClump(self):
        psfImage = self.getImage()
        #
        # Embed psfImage into a larger image so we can smooth when measuring it
        #
        width, height = psfImage.getWidth(), psfImage.getHeight()
        largeImg = psfImage.Factory(2*width, 2*height)
        largeImg.set(0)

        bbox = afwImage.BBox(afwImage.PointI(width, height), width, height)
        subLargeImg = psfImage.Factory(largeImg, bbox)
        subLargeImg <<= psfImage
        del subLargeImg
        # 
        # Now measure that image, looking for the highest peak.  
        # Start by building an Exposure 
        #
        msk = afwImage.MaskU(largeImg.getDimensions()) 
        msk.set(0)
        var = afwImage.ImageF(largeImg.getDimensions()) 
        var.set(1)
        mpsfImage = afwImage.MaskedImageF(largeImg, msk, var)
        mpsfImage.setXY0(afwImage.PointI(-width, -height))
        del msk 
        del var                                       
        exposure = afwImage.makeExposure(mpsfImage)
        
        #
        # Next run an object detector
        #
        stats = afwMath.makeStatistics(psfImage, afwMath.STDEV)
        
        threshold = afwDetection.Threshold(2*stats.getValue(afwMath.STDEV))
        ds = afwDetection.DetectionSetF(mpsfImage, threshold, "DETECTED")
        objects = ds.getFootprints()
        
        #
        # And measure it.  This policy isn't the one we use to measure
        # Sources, it's only used to characterize this PSF histogram
        #
        psfImagePolicy = policy.Policy()
        psfImagePolicy.add("centroidAlgorithm", "NAIVE")
        psfImagePolicy.add("shapeAlgorithm", "SDSS")
        psfImagePolicy.add("photometryAlgorithm", "NAIVE")
        psfImagePolicy.add("apRadius", 3.0)
        
        sigma = 1; 
        psf = algorithms.createPSF("DoubleGaussian", 1, 1, sigma)
        measureSources = algorithms.makeMeasureSources(exposure, 
                                                       psfImagePolicy, 
                                                       psf)
                                                       
        sourceList = afwDetection.SourceSet()
        
        Imax = None                     # highest peak
        for i in range(len(objects)):
            source = afwDetection.Source()
            sourceList.append(source)
            source.setId(i)
            
            try:
                measureSources.apply(source, objects[i])
            except Exception, e:
                continue
                
            x, y = source.getXAstrom(), source.getYAstrom()
            val = mpsfImage.getImage().get(int(x), int(y))
            
            if Imax is None or val > Imax:
                Imax = val
                psfClumpX, psfClumpY = x, y
                psfClumpIxx = source.getIxx()
                psfClumpIxy = source.getIxy()
                psfClumpIyy = source.getIyy()


        IzzMin = 0.5
        if psfClumpIxx < IzzMin or psfClumpIyy < IzzMin:
                psfClumpIxx, psfClumpIxy, psfClumpIyy = IzzMin, 0, IzzMin
        #
        # Convert psfClump[XY] (in psfImage's coordinates) back to Ixx/Iyy
        #
        psfClumpX, psfClumpY = self.peakToIxx(psfClumpX, psfClumpY)
        return psfClumpX, psfClumpY, psfClumpIxx, psfClumIxy, psfClumpIyy 
