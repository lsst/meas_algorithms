#
#
#
#

import lsst.pex.policy        as policy
import lsst.meas.algorithms   as measAlg
import lsst.afw.detection     as afwDetection

def sourceMeasurement(
    exposure,                 # exposure to analyse
    psf,                      # psf
    footprintLists,           # footprints of the detected objects
    measObjPolicy,            # measureObjects policy
    ):
    """ Source Measurement """

    # instantiate a measurement object for 
    # - instantiation only involves looking up the algorithms for centroid, shape, and photometry
    #   ... nothing actually gets measured yet.
    measureSources = measAlg.makeMeasureSources(exposure, measObjPolicy, psf)

    # create an empty list to contain the sources we found (as Source objects)
    sourceSet = afwDetection.SourceSet()
    
    for footprintList in footprintLists:
        #import pdb; pdb.set_trace() 
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

    return sourceSet
