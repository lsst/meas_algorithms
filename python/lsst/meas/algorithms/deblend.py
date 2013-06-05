#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
import math
import numpy

import lsst.pex.config as pexConf
import lsst.afw.table as afwTable
import lsst.pipe.base as pipeBase
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet

__all__ = 'SourceDeblendConfig', 'SourceDeblendTask'

class SourceDeblendConfig(pexConf.Config):
    psf_chisq_1 = pexConf.Field(dtype=float, default=1.5, optional=False,
                                doc=('Chi-squared per DOF cut for deciding a source is '+
                                     'a PSF during deblending (un-shifted PSF model)'))
    psf_chisq_2 = pexConf.Field(dtype=float, default=1.5, optional=False,
                                doc=('Chi-squared per DOF cut for deciding a source is '+
                                     'PSF during deblending (shifted PSF model)'))
    psf_chisq_2b = pexConf.Field(dtype=float, default=1.5, optional=False,
                                doc=('Chi-squared per DOF cut for deciding a source is '+
                                     'a PSF during deblending (shifted PSF model #2)'))
    maxNumberOfPeaks = pexConf.Field(dtype=int, default=0,
                                     doc=("Only deblend the brightest maxNumberOfPeaks peaks in the parent" +
                                          " (<= 0: unlimited)"))

class SourceDeblendTask(pipeBase.Task):
    """Split blended sources into individual sources.

    This task has no return value; it only modifies the SourceCatalog in-place.
    """
    ConfigClass = SourceDeblendConfig
    _DefaultName = "sourceDeblend"

    def __init__(self, schema, **kwargs):
        """Create the task, adding necessary fields to the given schema.

        @param[in,out] schema        Schema object for measurement fields; will be modified in-place.
        @param         **kwds        Passed to Task.__init__.
        """
        pipeBase.Task.__init__(self, **kwargs)

        self.nChildKey = schema.addField('deblend.nchild', type=int,
                                         doc='Number of children this object has (defaults to 0)')
        self.psfKey = schema.addField('deblend.deblended-as-psf', type='Flag',
                                      doc='Deblender thought this source looked like a PSF')
        self.psfCenterKey = schema.addField('deblend.psf-center', type='PointD',
                                         doc='If deblended-as-psf, the PSF centroid')
        self.psfFluxKey = schema.addField('deblend.psf-flux', type='D',
                                           doc='If deblended-as-psf, the PSF flux')
        self.tooManyPeaksKey = schema.addField('deblend.too-many-peaks', type='Flag',
                                               doc='Source had too many peaks; ' +
                                               'only the brightest were included')
        self.deblendFailedKey = schema.addField('deblend.failed', type='Flag',
                                                doc="Deblending failed on source")

        self.log.logdebug('Added keys to schema: %s' % ", ".join(str(x) for x in (
                    self.nChildKey, self.psfKey, self.psfCenterKey, self.psfFluxKey, self.tooManyPeaksKey)))

    @pipeBase.timeMethod
    def run(self, exposure, sources, psf):
        """Run deblend().

        @param[in]     exposure Exposure to process
        @param[in,out] sources  SourceCatalog containing sources detected on this exposure.
        @param[in]     psf      PSF

        @return None
        """
        self.deblend(exposure, sources, psf)

    @pipeBase.timeMethod
    def deblend(self, exposure, srcs, psf):
        """Deblend.
        
        @param[in]     exposure Exposure to process
        @param[in,out] srcs     SourceCatalog containing sources detected on this exposure.
        @param[in]     psf      PSF
                       
        @return None
        """
        self.log.info("Deblending %d sources" % len(srcs))

        from lsst.meas.deblender.baseline import deblend
        import lsst.meas.algorithms as measAlg

        # find the median stdev in the image...
        mi = exposure.getMaskedImage()
        stats = afwMath.makeStatistics(mi.getVariance(), mi.getMask(), afwMath.MEDIAN)
        sigma1 = math.sqrt(stats.getValue(afwMath.MEDIAN))

        schema = srcs.getSchema()

        n0 = len(srcs)
        nparents = 0
        for i,src in enumerate(srcs):
            fp = src.getFootprint()
            pks = fp.getPeaks()
            if len(pks) < 2:
                continue
            nparents += 1
            bb = fp.getBBox()
            xc = int((bb.getMinX() + bb.getMaxX()) / 2.)
            yc = int((bb.getMinY() + bb.getMaxY()) / 2.)
            if hasattr(psf, 'getFwhm'):
                psf_fwhm = psf.getFwhm(xc, yc)
            else:
                pa = measAlg.PsfAttributes(psf, xc, yc)
                psfw = pa.computeGaussianWidth(measAlg.PsfAttributes.ADAPTIVE_MOMENT)
                psf_fwhm = 2.35 * psfw

            self.log.logdebug('Parent %i: deblending %i peaks' % (int(src.getId()), len(pks)))

            self.preSingleDeblendHook(exposure, srcs, i, fp, psf, psf_fwhm, sigma1)
            npre = len(srcs)

            # This should really be set in deblend, but deblend doesn't have access to the src
            src.set(self.tooManyPeaksKey, len(fp.getPeaks()) > self.config.maxNumberOfPeaks)

            try:
                res = deblend(fp, mi, psf, psf_fwhm, sigma1=sigma1,
                              psf_chisq_cut1 = self.config.psf_chisq_1,
                              psf_chisq_cut2 = self.config.psf_chisq_2,
                              psf_chisq_cut2b= self.config.psf_chisq_2b,
                              maxNumberOfPeaks=self.config.maxNumberOfPeaks)
                src.set(self.deblendFailedKey, False)
            except Exception as e:
                self.log.warn("Error deblending source %d: %s" % (src.getId(), e))
                src.set(self.deblendFailedKey, True)
                continue

            kids = []
            nchild = 0
            for j,pkres in enumerate(res.peaks):
                if pkres.out_of_bounds:
                    # skip this source?
                    self.log.logdebug('Skipping out-of-bounds peak at (%i,%i)' %
                                      (pks[j].getIx(), pks[j].getIy()))
                    continue
                child = srcs.addNew(); nchild += 1
                child.setParent(src.getId())
                if hasattr(pkres, 'heavy'):
                    child.setFootprint(pkres.heavy)
                    #maskbits = pkres.heavy.getMaskBitsSet()
                    #print 'Mask bits set: 0x%x' % maskbits

                child.set(self.psfKey, pkres.deblend_as_psf)
                (cx,cy) = pkres.center
                child.set(self.psfCenterKey, afwGeom.Point2D(cx, cy))
                child.set(self.psfFluxKey, pkres.psfflux)
                kids.append(child)

            src.set(self.nChildKey, nchild)
            
            self.postSingleDeblendHook(exposure, srcs, i, npre, kids, fp, psf, psf_fwhm, sigma1, res)

        n1 = len(srcs)
        self.log.info('Deblended: of %i sources, %i were deblended, creating %i children, total %i sources' %
                      (n0, nparents, n1-n0, n1))

    def preSingleDeblendHook(self, exposure, srcs, i, fp, psf, psf_fwhm, sigma1):
        pass
    
    def postSingleDeblendHook(self, exposure, srcs, i, npre, kids, fp, psf, psf_fwhm, sigma1, res):
        pass

