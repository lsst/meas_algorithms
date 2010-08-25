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

"""Support utilities for Measuring sources"""

import re, sys
import lsst.pex.exceptions as pexExcept
import lsst.meas.algorithms as measAlg
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

def explainDetectionFlags(flags):
    """Return a string explaining Source's detectionFlags"""

    result = []
    for k, v in getDetectionFlags().items():
        if (flags & v):
            result += [k]

    result.sort()
    return " ".join(result)
    
def getDetectionFlags():
    """Return a dictionary of Source's detectionFlags"""

    flags = {}
    for k in measAlg.Flags.__dict__.keys():
        if not re.search(r"^[_A-Z0-9]+$", k): # flag names match this re
            continue

        flags[k] = measAlg.Flags.__dict__[k]

    return flags
    
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

def showPsfCandidates(exposure, psfCellSet, psf=None, frame=None, normalize=True):
    """Display the PSF candidates.  If psf is provided include PSF model and residuals;  if normalize is true normalize the PSFs (and residuals)"""
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
                    im = type(im)(im, True)
                    im.setXY0(cand.getImage().getXY0())
                    if normalize:
                        im /= cand.getAmplitude()
                except:
                    continue

                im_resid.append(im.getImage())

                model = psf.computeImage(afwGeom.makePointD(cand.getXCenter(), cand.getYCenter())).convertF()
                if not normalize:
                    model *= cand.getAmplitude()
                im_resid.append(model)

                resid = type(model)(model, True)
                if not False:
                    try:
                        xc = model.getX0() + model.getWidth()//2
                        yc = model.getY0() + model.getHeight()//2

                        centroider = measAlg.makeMeasureAstrometry(afwImage.makeExposure(afwImage.makeMaskedImage(model)))
                        centroider.addAlgorithm("GAUSSIAN")
                        
                        c = centroider.measure(afwDet.Peak(xc, yc)).find()
                        xmc, ymc = c.getX(), c.getY()

                        centroider.setImage(afwImage.makeExposure(im))
                        c = centroider.measure(afwDet.Peak(xc, yc)).find()
                        xc, yc = c.getX(), c.getY()

                        if False:
                            print "RHL %d %.2f %.2f  %.3f %.3f  %.3f %.3f" % \
                                  (cand.getSource().getId(),
                                   cand.getXCenter(), cand.getYCenter(),
                                   cand.getXCenter() - xc, cand.getYCenter() - yc,
                                   cand.getXCenter() -xmc, cand.getYCenter() -ymc)

                        if False:
                            print "RHL %d (%.1f, %.1f)   %g %g %g  %s" % (
                                cand.getSource().getId(),
                                cand.getXCenter(), cand.getYCenter(),
                                afwMath.makeStatistics(im.getImage(), afwMath.SUM).getValue(),
                                cand.getSource().getPsfFlux(),
                                cand.getAmplitude(),
                                explainDetectionFlags(cand.getSource().getFlagForDetection())), \
                                "Im", cand.getImage().get(0,0)

                                                               
                        resid = afwMath.offsetImage(resid, xc - cand.getXCenter(), yc - cand.getYCenter())
                    except pexExcept.LsstCppException, e:
                        print "RHL", e
                        #import pdb; pdb.set_trace() 
                        pass

                resid *= -1
                resid += im.getImage()
                im_resid.append(resid)

                if False:
                    im = type(im)(im, True); im.setXY0(cand.getImage().getXY0())
                    chi2 = measAlg.subtractPsf(psf, im, cand.getXCenter(), cand.getYCenter())
                    im_resid.append(im.getImage())

                # Fit the PSF components directly to the data (i.e. ignoring the spatial model)
                im = cand.getImage()

                im = type(im)(im, True)
                im.setXY0(cand.getImage().getXY0())
                if normalize:
                    im /= cand.getAmplitude()

                pair = measAlg.fitKernelToImage(afwMath.cast_LinearCombinationKernel(psf.getKernel()), im,
                                                afwGeom.makePointD(cand.getXCenter(), cand.getYCenter()))
                outputKernel, chisq = pair

                outImage = afwImage.ImageD(outputKernel.getDimensions())
                outputKernel.computeImage(outImage, False)
                if not False:
                    im -= outImage.convertF()
                    
                    im_resid.append(im.getImage())
                else:
                    im_resid.append(outImage.convertF())                    

                im = im_resid.makeMosaic()
            else:
                im = cand.getImage()

            im /= afwMath.makeStatistics(im, afwMath.MAX).getValue()
            mos.append(im, "%d %.1f" % (cand.getSource().getId(), rchi2),
                       ctype=ds9.RED if cand.isBad() else ds9.GREEN)

            import math                 # XXX
            if False and math.isnan(rchi2):
                ds9.mtv(cand.getImage().getImage(), title="candidate", frame=1)
                print "amp",  cand.getAmplitude()
                #import pdb; pdb.set_trace() 

            im = cand.getImage()
            candidateCenters.append((cand.getXCenter() - im.getX0(), cand.getYCenter() - im.getY0()))

    mosaicImage = mos.makeMosaic(frame=frame, title="Psf Candidates")

    i = 0
    for cen in candidateCenters:
        bbox = mos.getBBox(i); i += 1
        ds9.dot("+", cen[0] + bbox.getX0(), cen[1] + bbox.getY0(), frame=frame)

    return mosaicImage

def showPsf(psf, eigenValues=None, frame=None):
    """Display a PSF"""

    mos = displayUtils.Mosaic()
    i = 0
    for k in afwMath.cast_LinearCombinationKernel(psf.getKernel()).getKernelList():
        im = afwImage.ImageD(k.getDimensions())
        k.computeImage(im, False)
        if eigenValues:
            mos.append(im, "%g" % eigenValues[i])
            i += 1
        else:
            mos.append(im)

    mos.makeMosaic(frame=frame, title="Eigen Images")

    return mos

def showPsfMosaic(exposure, psf, nx=7, ny=None, frame=None):
    mos = displayUtils.Mosaic()

    try:                                # maybe it's a real Exposure
        width, height = exposure.getWidth(), exposure.getHeight()
    except AttributeError:
        try:                            # OK, maybe a list [width, height]
            width, height = exposure[0], exposure[1]
        except TypeError:               # I guess not
            raise RuntimeError, ("Unable to extract width/height from object of type %s" % type(exposure))

    if not ny:
        ny = int(nx*float(height)/width + 0.5)
        if not ny:
            ny = 1

    centroider = measAlg.makeMeasureAstrometry(None)
    centroider.addAlgorithm("GAUSSIAN")

    centers = []
    for ix in range(nx):
        for iy in range(ny):
            x = int((ix + 0.5)*width/nx)
            y = int((iy + 0.5)*height/ny)

            im = psf.computeImage(afwGeom.makePointD(x, y)).convertF()
            mos.append(im, "PSF(%d,%d)" % (x, y))
    
            centroider.setImage(afwImage.makeExposure(afwImage.makeMaskedImage(im)))
            w, h = im.getDimensions()
            c = centroider.measure(afwDet.Peak(im.getX0() + w//2, im.getY0() + h//2)).find()

            centers.append((c.getX() - im.getX0(), c.getY() - im.getY0()))

    mos.makeMosaic(frame=frame, title="Model Psf", mode=nx)

    if centers:
        i = 0
        for cen in centers:
            bbox = mos.getBBox(i); i += 1
            ds9.dot("+", cen[0] + bbox.getX0(), cen[1] + bbox.getY0(), frame=frame)

    return mos

def writeSourceSetAsCsv(sourceSet, fd=sys.stdout):
    """Write a SourceSet as a CSV file"""

    if not sourceSet:
        raise RuntimeError, "Please provide at least one Source"

    source = sourceSet[0]

    measurementTypes = (("astrometry", source.getAstrometry),
                        ("photometry", source.getPhotometry),
                        ("shape", source.getShape),
                        )

    print >> fd, "#misc::id:int:1"

    for measureType, getWhat in measurementTypes:
        a = getWhat()

        for value in a.getValues():
            for s in value.getSchema():
                if s.getType() == s.LONG:
                    typeName = "long"
                else:
                    typeName = "double"

                if s.isArray():
                    n = s.getDimen()
                else:
                    n = 1
                print >> fd, "#%s:%s:%s:%s:%d" % (measureType, value.getAlgorithm(), s.getName(),
                                                  typeName, n)


    for source in sourceSet:
        out = "%d" % (source.getId())

        for a in (source.getAstrometry(),
                  source.getPhotometry(),
                  source.getShape(),
                  ):
            for value in a.getValues():
                for sch in value.getSchema():
                    if out:
                        out += ", "
                    out += str(value.get(sch.getName()))

        print >> fd, out

