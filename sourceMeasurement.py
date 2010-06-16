import lsst.pex.policy        as policy
import lsst.meas.algorithms   as measAlg
import lsst.afw.detection     as afwDetection
import lsst.afw.geom          as afwGeom
import lsst.afw.display.ds9  as ds9

def sourceMeasurement(
    exposure,                 # exposure to analyse
    psf,                      # psf
    footprintLists,           # footprints of the detected objects
    measObjPolicy,            # measureObjects policy
    ):
    """ Source Measurement """

    try:
        import lsstDebug

        display = lsstDebug.Info(__name__).display
    except ImportError, e:
        try:
            display
        except NameError:
            display = False

    if display:
        frame = 0
        ds9.mtv(exposure, title="input", frame=frame)

    # instantiate a measurement object for 
    # - instantiation only involves looking up the algorithms for centroid, shape, and photometry
    #   ... nothing actually gets measured yet.
    measureSources = measAlg.makeMeasureSources(exposure, measObjPolicy, psf)

    # create an empty list to contain the sources we found (as Source objects)
    sourceSet = afwDetection.SourceSet()
    
    for footprintList in footprintLists:
        footprints, isNegative = footprintList

        # loop over all the objects detected
        for i in range(len(footprints)):

            # create a new source, and add it to the list, initialize ...
            source = afwDetection.Source()
            sourceSet.append(source)
            source.setId(i)

            source.setFlagForDetection(source.getFlagForDetection() | measAlg.Flags.BINNED1);

            # actually try to "measure" this object
            # recall: footprints[i] is a footprint for an object, measured values will be recorded in 'source'
            try:
                measureSources.apply(source, footprints[i])
            except Exception, e:
                # logging might be nice here
                #self.log.log(Log.WARN, str(e))
                pass
            #
            # Set the time
            #
            if False and exposure.getDetector():
                pos = afwGeom.makePointI(int(source.getXAstrom()), int(source.getYAstrom()))
                midTime = exposure.getCalib().getMidTime(exposure.getDetector(), pos)
            else:
                midTime = exposure.getCalib().getMidTime()
                
            source.setTaiMidPoint(midTime.get())
            source.setTaiRange(exposure.getCalib().getExptime())

            if display:
                ds9.dot("+", source.getXAstrom(), source.getYAstrom(), size=3, ctype=ds9.RED)
                if display > 1:
                    ds9.dot(("@:%.1f,%.1f,%1f" % (source.getXAstromErr()**2, 0, source.getYAstromErr()**2)),
                            source.getXAstrom(), source.getYAstrom(), size=3, ctype=ds9.RED)

    return sourceSet
