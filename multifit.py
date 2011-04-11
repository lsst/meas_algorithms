import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult
import lsst.afw.image as afwImg
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPol
import lsst.afw.geom.ellipses as geomEllipses
import lsst.pex.exceptions
import math
import numpy
import eups
import os.path

def makePsModel(source):
    if(numpy.isnan(source.getApFlux()) or \
            numpy.isnan(source.getXAstrom()) or numpy.isnan(source.getYAstrom())):
        flags = source.getFlagForDetection()
        source.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_PS_NAN)
        raise Excetion("Initial PS parameters are NaN")

    pixel = afwGeom.Point2D(source.getXAstrom(), source.getYAstrom())
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

    pixel = afwGeom.Point2D(source.getXAstrom(), source.getYAstrom())
    quad = afwGeom.ellipses.Quadrupole(
            source.getIxx(), 
            source.getIyy(), 
            source.getIxy())
    if quad.getDeterminant() < 0:
        #moments are bad
        flags= source.getFlagForDetection()
        source.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_SG_MOMENTS)
        raise Exception("Initial SG moments are bad")
    if afwGeom.ellipses.Distortion(quad).getE() > 0.95:
        flags= source.getFlagForDetection()
        source.setFlagForDetection(flags | measMult.Flags.FAIL_INIT_SG_DISTORTION)
        raise Exception("high ellipticity on construction")

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

def fitSourceModels(mi, psf, sources, footprints, policy, cacheName="sersicCache.boost"):    
    #some setup first. 
    #Make sure we have a cache.
    cacheDir = os.path.join(eups.productDir("meas_multifitData"), "cache")
    cachePath = os.path.join(cacheDir, cacheName)
    cache = measMult.SersicCache.load(cachePath)
    measMult.SersicMorphology.setSersicCache(cache)
   
    fitterPol = policy.getPolicy("fitterPolicy")
    minPixels = policy.getInt("minPixels")
    
    if policy.exists("fourierModelProjectionPolicy"):
        fourierPol = policy.getPolicy("fourierModelProjectionPolicy")
        measMult.FourierModelProjection.setPolicy(fourierPol)

    fitter= measMult.MinuitNumericFitter(fitterPol)
    identityTransform = afwGeom.AffineTransform();
    nSkipPs=0
    nSkipSg=0
    nFailPs=0
    nFailSg=0
    i = 0
    for s, fp in zip(sources, footprints):
        i +=1
        flags = s.getFlagForDetection()
        try:
            psModel, psErrs = makePsModel(s)
            psEval = measMult.ModelEvaluator(psModel, minPixels)            
            nPix = psEval.setData(mi, psf, identityTransform, fp)

            if nPix > 0:
                try:
                    psResult = fitter.apply(psEval, psErrs)
                except:
                    flags |= measMult.Flags.FAIL_FIT_PS
                    nFailPs += 1                    
                else:
                    if (psResult.flags & measMult.MinuitFitterResult.MAX_ITERATION_REACHED) != 0:                    
                        flags |= measMult.Flags.PS_MAX_ITERATIONS
                    if (psResult.flags & measMult.MinuitFitterResult.CONVERGED) == 0:
                        falgs |= measMult.Flags.FAIL_FIT_PS
                        nFailPs += 1
                    flux = psResult.parameters[0]
                    fluxErr = math.sqrt(psResult.covariance[0,0])
                    psPhotom = measMult.PointSourceModelPhotometry(flux, fluxErr)
                    s.getPhotometry().add(psPhotom)
            else:
                flags |= measMult.Flags.FAIL_INIT_PS_NO_PIX
                nSkipPs += 1
        except Exception, e:
            #errors in construction
            nSkipPs += 1          

        try:
            sgModel, sgErrs = makeSgModel(s)
            sgEval = measMult.ModelEvaluator(sgModel, minPixels)            
            nPix= sgEval.setData(mi, psf, identityTransform, fp)

            if nPix > 0:
                try:
                    sgResult = fitter.apply(sgEval, sgErrs)
                except lsst.pex.exceptions.LsstCppException, e:
                    if e.args[0].__class__ == measMult.ParameterRangeException:
                        if not e.args[0].isOutOfRange(0):
                            flags |= measMult.Flags.FAIL_FIT_SG_DISTORTION
                        if not e.args[0].isOutOfRange(3):
                            flags |= measMult.Flags.FAIL_FIT_SG_SERSIC
                    else:
                        flags |= measMult.Flags.FAIL_FIT_SG_FOURIER_AREA
                        nFailSg += 1
                except Exception,e :
                    falgs |= measMult.Flags.FAIL_FIT_SG_OPTIMIZER
                    nFailSg += 1
                else:
                    if (sgResult.flags & measMult.MinuitFitterResult.MAX_ITERATION_REACHED) != 0:                    
                        falgs |= measMult.Flags.SG_MAX_ITERATIONS
                    if (sgResult.flags & measMult.MinuitFitterResult.CONVERGED) == 0:
                        nFailSg+=1
                        flags |= measMult.Flags.FAIL_FIT_SG_OPTIMIZER
                    sgPhotom = measMult.SmallGalaxyModelPhotometry(sgResult.parameters, sgResult.covariance, \
                            cache.getInnerSersicRadius(), cache.getOuterSersicRadius())
                    s.getPhotometry().add(sgPhotom)
            else:           
                flags |= measMult.Flags.FAIL_INIT_SG_NO_PIX
                nSkipSg +=1
        except Exception, e:
            #errors in construction
            nSkipSg +=1
        
        s.setFlagForDetection(flags)
        
    #return some stats about fitting
    return (len(sources), nSkipPs, nFailPs, nSkipSg, nFailSg)
