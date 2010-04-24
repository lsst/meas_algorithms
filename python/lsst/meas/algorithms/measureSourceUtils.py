"""Support utilities for Measuring sources"""

import re
import lsst.pex.exceptions as pexExcept
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
def showPsfSpatialCells(exposure, psfCellSet, nMaxPerCell=-1, showChi2=False,
                        symb=None, ctype=None, size=2, frame=None):
    """Show the SpatialCells.  If symb is something that ds9.dot understands (e.g. "o"), the top nMaxPerCell candidates will be indicated with that symbol, using ctype and size"""

    origin = [-exposure.getMaskedImage().getX0(), -exposure.getMaskedImage().getY0()]
    for cell in psfCellSet.getCellList():
        displayUtils.drawBBox(cell.getBBox(), origin=origin, frame=frame)

        if nMaxPerCell < 0:
            nMaxPerCell = 0

        i = 0
        for cand in cell.begin(True):
            if nMaxPerCell > 0:
                i += 1

            cand = measAlg.cast_PsfCandidateF(cand)

            xc, yc = cand.getXCenter() + origin[0], cand.getYCenter() + origin[1]

            if i <= nMaxPerCell and symb:
                ds9.dot(symb, xc, yc, frame=frame, ctype=ctype, size=size)

            if showChi2:
                nu = cand.getWidth()*cand.getHeight() - 1 # number of dof/star for chi^2
                ds9.dot("%.1f" % (cand.getChi2()/nu), xc-size, yc-size, frame=frame, ctype=ctype, size=size)

def showPsfCandidates(exposure, psfCellSet, psf=None, frame=None):
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

    candidateCenters = []
    for cell in psfCellSet.getCellList():
        for cand in cell.begin(False): # include bad candidates
            cand = measAlg.cast_PsfCandidateF(cand)

            rchi2 = cand.getChi2()/nu

            if False and cand.isBad():
                continue

            if psf:
                im_resid = displayUtils.Mosaic(gutter=0, background=-5, mode="x")

                try:
                    im = cand.getImage()
                except:
                    continue
                im_resid.append(im.getImage())

                model = psf.getImage(cand.getXCenter(), cand.getYCenter()).convertF()
                model *= cand.getAmplitude()
                im_resid.append(model)

                resid = type(model)(model, True)
                if False:
                    try:
                        centroider = measAlg.createMeasureCentroid("GAUSSIAN")
                        c = centroider.apply(model,
                                             model.getX0() + model.getWidth()//2,
                                             model.getY0() + model.getHeight()//2, None, 0.0)
                        xmc, ymc = c.getX(), c.getY()
                        c = centroider.apply(im.getImage(),
                                             model.getX0() + model.getWidth()//2,
                                             model.getY0() + model.getHeight()//2, None, 0.0)
                        xc, yc = c.getX(), c.getY()
                        print "RHL %d %.2f %.2f  %.3f %.3f  %.3f %.3f" % (cand.getSource().getId(),
                                                               cand.getXCenter(), cand.getYCenter(),
                                                               cand.getXCenter() - xc, cand.getYCenter() - yc,
                                                               cand.getXCenter() -xmc, cand.getYCenter() - ymc)
                                                               
                        resid = afwMath.offsetImage(resid, xc - cand.getXCenter(), yc - cand.getYCenter())
                        #resid = afwMath.offsetImage(resid, xc - xmc, yc - ymc)
                    except pexExcept.LsstCppException, e:
                        print "RHL", e
                        pass

                resid *= -1
                resid += im.getImage()
                im_resid.append(resid)

                if False:
                    im = type(im)(im, True); im.setXY0(cand.getImage().getXY0())
                    chi2 = measAlg.subtractPsf(psf, im, cand.getXCenter(), cand.getYCenter())
                    im_resid.append(im.getImage())

                im = im_resid.makeMosaic()
            else:
                im = cand.getImage()

            mos.append(im, "%d %.1f" % (cand.getSource().getId(), rchi2),
                       ctype=ds9.RED if cand.isBad() else ds9.GREEN)

            im = cand.getImage()
            candidateCenters.append((cand.getXCenter() - im.getX0(), cand.getYCenter() - im.getY0()))

    mosaicImage = mos.makeMosaic(frame=frame, title="Psf Candidates")

    i = 0
    for cen in candidateCenters:
        bbox = mos.getBBox(i); i += 1
        ds9.dot("+", cen[0] + bbox.getX0(), cen[1] + bbox.getY0(), frame=frame)

    return mosaicImage

def showPsf(psf, frame=None):
    """Display a PSF"""

    mos = displayUtils.Mosaic()
    for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
        im = afwImage.ImageD(k.getDimensions())
        k.computeImage(im, False)
        mos.append(im)

    mos.makeMosaic(frame=frame, title="Eigen Images")

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

    centroider = measAlg.createMeasureCentroid("GAUSSIAN")

    centers = []
    for ix in range(nx):
        for iy in range(ny):
            x = int((ix + 0.5)*width/nx)
            y = int((iy + 0.5)*height/ny)

            im = psf.getImage(x, y).convertF()
            mos.append(im, "PSF(%d,%d)" % (x, y))

            w, h = im.getDimensions()
            c = centroider.apply(im, im.getX0() + w//2, im.getY0() + h//2)
            centers.append((c.getX() - im.getX0(), c.getY() - im.getY0()))

    mos.makeMosaic(frame=frame, title="Model Psf", mode=nx)

    if centers:
        i = 0
        for cen in centers:
            bbox = mos.getBBox(i); i += 1
            ds9.dot("+", cen[0] + bbox.getX0(), cen[1] + bbox.getY0(), frame=frame)

    return mos
