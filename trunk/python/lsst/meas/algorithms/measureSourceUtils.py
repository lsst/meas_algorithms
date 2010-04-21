"""Support utilities for Measuring sources"""

import re
import lsst.meas.algorithms as measAlg
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

def explainDetectionFlags(flags):
    """Return a string explaining Source's detectionFlags"""

    result = []
    for k in measAlg.Flags.__dict__.keys():
        if not re.search(r"^[_A-Z0-9]+$", k): # flag names match this re
            continue

        if (flags & measAlg.Flags.__dict__[k]):
            result += [k]

    result.sort()
    return " ".join(result)
    
def showSourceSet(sSet, xy0=(0, 0), frame=0, ctype=ds9.GREEN, symb="+", size=2):
    """Draw the (XAstrom, YAstrom) positions of a set of Sources.  Image has the given XY0"""
    for s in sSet:
        ds9.dot(symb, s.getXAstrom() - xy0[0], s.getYAstrom() - xy0[1], frame=frame, ctype=ctype, size=size)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#
# PSF display utilities
#
def showPsfCandidates(exposure, psfCellSet, frame=None):
    """Display the PSF candidates"""
    #
    # Show us the ccandidates
    #
    mos = displayUtils.Mosaic()
    #
    # Instantiate a psfCandidate so we can use makePsfCandidate to determine the correct type
    #
    psfCandidate = measAlg.makePsfCandidate(afwDet.Source(), exposure.getMaskedImage())
    nu = psfCandidate.getWidth()*psfCandidate.getHeight() - 1 # number of dof/star for chi^2
    del psfCandidate

    stamps = []; stampInfo = []
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False): # include bad candidates
            cand = measAlg.cast_PsfCandidateF(cand)

            rchi2 = cand.getChi2()/nu

            if not cand.isBad():
                im = cand.getImage()
                stamps.append(im)
                stampInfo.append("%d %.1f" % (cand.getSource().getId(), rchi2))

    mos.makeMosaic(stamps, frame=frame, title="Psf Candidates")
    mos.drawLabels(stampInfo, frame=frame)

    return mos

def showPsf(psf, frame=None):
    """Display a PSF"""

    mos = displayUtils.Mosaic()
    eigenImages = []
    for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
        im = afwImage.ImageD(k.getDimensions())
        k.computeImage(im, False)
        eigenImages.append(im)

    mos.makeMosaic(eigenImages, frame=frame, title="Psf")

    return mos

def showPsfMosaic(exposure, psf, nx=7, ny=7, frame=None):
    mos = displayUtils.Mosaic()

    try:                                # maybe it's a real Exposure
        width, height = exposure.getWidth(), exposure.getHeight()
    except AttributeError:
        try:                            # OK, maybe a list [width, height]
            width, height = exposure[0], exposure[1]
        except TypeError:               # I guess not
            raise RuntimeError, ("Unable to extract width/height from object of type %s" % type(exposure))
        
    psfImages = []
    labels = []
    for ix in range(nx):
        for iy in range(ny):
            x = (ix + 0.5)*width/nx
            y = (iy + 0.5)*height/ny

            psfImages.append(psf.getImage(x, y))
            labels.append("PSF(%d,%d)" % (int(x), int(y)))

    mos.makeMosaic(psfImages, frame=frame, title="Psf", mode=nx)
    mos.drawLabels(labels, frame=frame)

    return mos
