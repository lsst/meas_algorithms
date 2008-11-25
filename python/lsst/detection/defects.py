"""Support for image defects"""

import lsst.afw.image.imageLib as afwImage
import lsst.pex.policy as policy
import lsst.detection.detectionLib as detection

def policyToBadRegionList(policyFile):
    """Given a Policy file describing a CCD's bad pixels, return a vector of BadRegion::Ptr""" 

    badPixelsPolicy = policy.Policy.createPolicy(policyFile)
    badPixels = detection.DefectListT()
    
    d = badPixelsPolicy.getArray("Defects")
    for reg in d:
        x0 = reg.get("x0")
        width = reg.get("width")
        if not width:
            x1 = reg.get("x1")
            width = x1 - x0 - 1

        y0 = reg.get("y0")
        height = reg.get("height")
        if not height:
            y1 = reg.get("y1")
            height = y1 - y0 - 1

        bbox = afwImage.BBox(afwImage.PointI(x0, y0), width, height)
        badPixels.push_back(detection.Defect(bbox))

    del badPixelsPolicy

    return badPixels
