# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

        self.psfkey = schema.addField('deblend.deblended-as-psf', type='Flag',
                                      doc='Deblender thought this source looked like a PSF')
        self.psf_xykey = schema.addField('deblend.psf-center', type='PointD',
                                         doc='If deblended-as-psf, the PSF centroid')
        self.psf_fluxkey = schema.addField('deblend.psf-flux', type='D',
                                           doc='If deblended-as-psf, the PSF flux')
        #self.deblended_at_edge = schema.addField('deblend.deblended-at-edge', type='Flag',
        #                                         doc='This source is near an edge so the deblender had to guess about the profiles.')

        # self.deblend_failed = ...

        self.log.logdebug('Added keys to schema: ' + str(self.psfkey) + ', ' + str(self.psf_xykey)
                          + ', ' + str(self.psf_fluxkey))

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
        # masksToFlag = [(maskname, flagname,
        #                 schema.find(flagname).key, mi.getMask().getPlaneBitMask(maskname))
        #                for maskname,flagname in [('EDGE', 'flags.pixel.edge'),
        #                                          ('INTRP', 'flags.pixel.interpolated.any'),
        #                                          ('SAT', 'flags.pixel.saturated.any'),
        #                                          ]]
        # maskBitsToFlag = sum([val for mn,kn,key,val in masksToFlag])

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

            # ekey = schema.find('flags.pixel.edge').key
            # if src.get(ekey):
            #     print 'Parent has EDGE flag set'
            # ph = afwDet.makeHeavyFootprint(fp, mi)
            # maskbits = ph.getMaskBitsSet()
            # #print 'fp', fp
            # #print 'ph', ph
            # #print 'mi', mi
            # m = ph.getMaskArray()
            # #print 'MaskArray m:', m
            # mb2 = 0
            # for mm in m:
            #     mb2 |= mm
            # print 'Maskbits', maskbits
            # #print 'Maskbits2', mb2
            # assert(maskbits == mb2)
            # if maskbits & mi.getMask().getPlaneBitMask('EDGE'):
            #     print 'Parent has EDGE mask pixels set'
            # #print 'Parent id %i: Mask bits set: 0x%x' % (src.getId() & 0xffff, maskbits)

            res = deblend(fp, mi, psf, psf_fwhm, sigma1=sigma1,
                          psf_chisq_cut1 = self.config.psf_chisq_1,
                          psf_chisq_cut2 = self.config.psf_chisq_2,
                          psf_chisq_cut2b= self.config.psf_chisq_2b)
            kids = []
            for j,pkres in enumerate(res.peaks):
                if pkres.out_of_bounds:
                    # skip this source?
                    self.log.logdebug('Skipping out-of-bounds peak at (%i,%i)' %
                                      (pks[j].getIx(), pks[j].getIy()))
                    continue
                child = srcs.addNew()
                child.setParent(src.getId())
                if hasattr(pkres, 'heavy'):
                    child.setFootprint(pkres.heavy)
                    maskbits = pkres.heavy.getMaskBitsSet()
                    #print 'Mask bits set: 0x%x' % maskbits

                    # This should get done by the usual method (meas_alg's PixelFlags.cc)
                    # if maskbits & maskBitsToFlag:
                    #     for maskname,keynm,key,bitval in masksToFlag:
                    #         if bitval & maskbits == 0:
                    #             continue
                    #         #print 'mask bit', maskname, 'is set; setting key', keynm
                    #         child.set(key, True)

                child.set(self.psfkey, pkres.deblend_as_psf)
                (cx,cy) = pkres.center
                child.set(self.psf_xykey, afwGeom.Point2D(cx, cy))
                child.set(self.psf_fluxkey, pkres.psfflux)
                kids.append(child)
                
            self.postSingleDeblendHook(exposure, srcs, i, npre, kids, fp, psf, psf_fwhm, sigma1, res)

        n1 = len(srcs)
        self.log.info('Deblended: of %i sources, %i were deblended, creating %i children, total %i sources' % (n0, nparents, n1-n0, n1))

    def preSingleDeblendHook(self, exposure, srcs, i, fp, psf, psf_fwhm, sigma1):
        pass
    
    def postSingleDeblendHook(self, exposure, srcs, i, npre, kids, fp, psf, psf_fwhm, sigma1, res):
        pass

