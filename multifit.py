import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult
import lsst.afw.image as afwImg
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPol
import lsst.afw.geom.ellipses as geomEllipses
import math
import numpy
import eups
import os.path

def makePsModel(source):
    if(numpy.isnan(source.getApFlux()) or \
            numpy.isnan(source.getXAstrom()) or numpy.isnan(source.getYAstrom())):
        flags = source.getFlagForDetection()
        source.setFlagForDetection(flags | measMult.Flags.FAIL_INITIALIZE_PS_NAN)
        raise Excetion("Initial PS parameters are NaN")

    pixel = afwGeom.makePointD(source.getXAstrom(), source.getYAstrom())
    astrometry = measMult.FixedAstrometry(pixel)
    morphology = measMult.createPointSourceMorphology(source.getApFlux())
    model = measMult.ComponentModel.create(astrometry, morphology)
    errs = [source.getApFluxErr()]
    if errs[0] == 0.0 or numpy.isnan(errs[0]):
        errs[0] = 1e-1

    return (model, errs)

def makeSgModel(source):
    if(numpy.isnan(source.getXAstrom()) or numpy.isnan(source.getYAstrom()) or \
            numpy.isnan(source.getYAstrom()) or numpy.isnan(source.getIxx()) or \
            numpy.isnan(source.getIxy()) or numpy.isnan(source.getIyy()) or \
            numpy.isnan(source.getApFlux())):
        flags= source.getFlagForDetection()
        source.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_SG_NAN)
        raise Exception("Initial SG parameters are NaN")

    pixel = afwGeom.makePointD(source.getXAstrom(), source.getYAstrom())
    quad = afwGeom.ellipses.Quadrupole(
            source.getIxx(), 
            source.getIyy(), 
            source.getIxy())
    if quad.getDeterminant() < 0:
        #moments are bad
        flags = source.getFlagForDetection()
        source.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_SG_MOMENTS)
        raise Exception("Initial SG moments are bad")

    #use fixed sersic Index
    #TODO compute this instead from source/fp instead
    n = 1.5

    astrometry = measMult.FixedAstrometry(pixel)
    morphology = measMult.createSersicMorphology(source.getApFlux(), quad, n)
    model = measMult.ComponentModel.create(astrometry, morphology)
    #use fixed ellipse and sersic errors
    #avoids issues of missing or invalid measurement errors
    errs = [source.getApFluxErr(), 1e-4, 1e-4, 1e-2, 1e-6]
    if errs[0] == 0.0 or numpy.isnan(errs[0]):
        errs[0] = 1e-1

    return (model, errs)

def fitSourceModels(mi, psf, sources, footprints, fitterPolicy, cacheName="sersicCache.dat"):    
    #some setup first. 
    #Make sure we have a cache.
    cacheDir = os.path.join(eups.productDir("meas_multifitData"), "cache")
    cachePath = os.path.join(cacheDir, cacheName)
    cache = measMult.Cache.load(cachePath, "SersicCache", False)
    measMult.SersicMorphology.setSersicCache(cache)

    fitter= measMult.MinuitNumericFitter(fitterPolicy)
    identityTransform = afwGeom.AffineTransform();
    for s, fp in zip(sources, footprints):
        flags = s.getFlagForDetection()
        try:
            psModel, psErrs = makePsModel(s)
            psEval = measMult.ModelEvaluator(psModel)            
            nPix = psEval.setData(mi, psf, identityTransform, fp)

            if nPix > 0:
                try:
                    psResult = fitter.apply(psEval, psErrs)
                except:
                    s.setFlagForDetection(flags | measMult.Flags.FAIL_FIT_PS)
                    print "fit fail ps"
                else:
                    if (psResult.flags & measMult.MinuitFitterResult.MAX_ITERATION_REACHED) != 0:                    
                        s.setFlagForDetection(flags | measMult.Flags.PS_MAX_ITERATIONS)
                    if (psResult.flags & measMult.MinuitFitterResult.CONVERGED) == 0:
                        s.setFlagForDetection(flags | measMult.Flags.FAIL_FIT_PS)
                    flux = psResult.parameters[0]
                    fluxErr = math.sqrt(psResult.covariance[0,0])
                    psPhotom = measMult.PointSourceModelPhotometry(flux, fluxErr)
                    s.getPhotometry().add(psPhotom)
            else:
                s.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_PS_NO_PIX)
                print "no pix ps"
        except Exception, e:
            #TODO: log this
            print "unknown err ps", e

        try:
            sgModel, sgErrs = makeSgModel(s)
            sgEval = measMult.ModelEvaluator(sgModel)            
            nPix= sgEval.setData(mi, psf, identityTransform, fp)

            if nPix > 0:
                try:
                    sgResult = fitter.apply(sgEval, sgErrs)                    
                except measMult.ParameterRangeException, e:
                    if e.notInRange(3):
                        s.setflagsForDetection(flags | measMult.Flags.FAIL_FIT_SG_SERSIC)
                    if e.notInRange(0) or e.notInRange(1):
                        s.setFlagForDetection(flags | measMult.Flags.FAIL_FIT_SG_DISTORTION)
                    print "param range sg"
                except Exception,e :
                    s.setFlagForDetection(flags | measMult.Flags.FAIL_FIT_SG_OPTIMIZER)
                    print "minuit sg", e
                else:
                    if (sgResult.flags & measMult.MinuitFitterResult.MAX_ITERATION_REACHED) != 0:                    
                        s.setFlagForDetection(flags | measMult.Flags.SG_MAX_ITERATIONS)
                    if (sgResult.flags & measMult.MinuitFitterResult.CONVERGED) == 0:
                        s.setFlagForDetection(flags | measMult.Flags.FAIL_FIT_SG_OPTIMIZER)
                    sgPhotom = measMult.SmallGalaxyModelPhotometry(sgResult.parameters, sgResult.covariance)
                    s.getPhotometry().add(sgPhotom)
            else:           
                s.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_SG_NO_PIX)
                print "no pix sg"
        except Exception, e:
            print "unknown err sg", e
            #TODO log this
