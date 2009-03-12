import lsst.afw.Image as afwImage

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
            if i == 0 and j == 0:
                return

            self._psfImage.set(i, j, self._psfImage.get(i, j) + 1)

            if False:
                print "Inserting %d at (%d, %d)" % (source.getId(), i, j),
                print "(%d, %d) (flux = %.0f), (%.1f %.1f)" % (source.getXAstrom(), source.getYAstrom(),
                                                                                    source.getPsfFlux(),
                                                                                    source.getIxx(), source.getIyy())

    def peakToIxx(self, peakX, peakY):
        """Given a peak position in self._psfImage, return the corresponding (Ixx, Iyy)"""
        xx = peakX*self._xMax/self._xSize
        yy = peakY*self._yMax/self._ySize
        return xx, yy


    def getCentroid(self):
        pass

    def getShape(self):
        pass
