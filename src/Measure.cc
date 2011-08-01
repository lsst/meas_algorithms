// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/// \file

#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/algorithms/Measure.h"

namespace lsst {
namespace meas {
namespace algorithms {
    
namespace pexLogging = lsst::pex::logging;
namespace afwImage = lsst::afw::image;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;

namespace {
    /*
     * Return the numeric value of name as double
     */
    double getNumeric(lsst::pex::policy::Policy const& policy, std::string const& name)
    {
        return policy.isDouble(name) ? policy.getDouble(name) : policy.getInt(name);
    }
}

/************************************************************************************************************/
/**
 * @brief Calculate a detected source's moments
 */
template <typename MaskedImageT>
class FootprintCentroid : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintCentroid(MaskedImageT const& mimage ///< The image the source lives in
                              ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                                  _n(0), _sum(0), _sumx(0), _sumy(0),
                                  _min( std::numeric_limits<double>::max()), _xmin(0), _ymin(0),
                                  _max(-std::numeric_limits<double>::max()), _xmax(0), _ymax(0),
                                  _bits(0) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _n = 0;
        _sum = _sumx = _sumy = 0.0;
        _min =  std::numeric_limits<double>::max();
        _xmin = _ymin = 0;
        _max = -std::numeric_limits<double>::max();
        _xmax = _ymax = 0;
        _bits = 0x0;
    }
    virtual void reset(afwDetection::Footprint const&) {}

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);

        _n++;
        _sum += val;
        _sumx += afwImage::indexToPosition(x)*val;
        _sumy += afwImage::indexToPosition(y)*val;
        _bits |= loc.mask(0, 0);

        if (val < _min) {
            _min = val;
            _xmin = x;
            _ymin = y;
        }
        if (val > _max) {
            _max = val;
            _xmax = x;
            _ymax = y;
        }
    }

    /// Return the number of pixels
    int getN() const { return _n; }
    /// Return the Footprint's flux
    double getSum() const { return _sum; }
    /// Return the Footprint's column centroid
    double getX() const { return _sumx/_sum; }
    /// Return the Footprint's row centroid
    double getY() const { return _sumy/_sum; }
    /// Return the Footprint's peak pixel
    PTR(afwDetection::Peak) makePeak(bool isNegative) const {
        return boost::make_shared<afwDetection::Peak>(isNegative ? afwDetection::Peak(_xmin, _ymin) :
                                                      afwDetection::Peak(_xmax, _ymax));
    }
    /// Return the union of the bits set anywhere in the Footprint
    typename MaskedImageT::Mask::Pixel getBits() const { return _bits; }
private:
    int _n;
    double _sum, _sumx, _sumy;
    double _min;
    int _xmin, _ymin;
    double _max;
    int _xmax, _ymax;
    typename MaskedImageT::Mask::Pixel _bits;
};


template <typename MaskedImageT>
class FootprintFlux : public afwDetection::FootprintFunctor<MaskedImageT> {
public:
    explicit FootprintFlux(MaskedImageT const& mimage ///< The image the source lives in
                 ) : afwDetection::FootprintFunctor<MaskedImageT>(mimage),
                     _sum(0) {}

    /// \brief Reset everything for a new Footprint
    void reset() {
        _sum = 0.0;
    }

    /// \brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator loc, ///< locator pointing at the pixel
                    int x,                                 ///< column-position of pixel
                    int y                                  ///< row-position of pixel
                   ) {
        typename MaskedImageT::Image::Pixel val = loc.image(0, 0);
        _sum += val;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

private:
    double _sum;
};

/************************************************************************************************************/
/**
 * Use *this to measure the Footprint foot, setting fields in src
 */
template<typename ExposureT>
void MeasureSources<ExposureT>::apply(
	PTR(lsst::afw::detection::Source) src,       ///< the Source to receive results
        CONST_PTR(lsst::afw::detection::Footprint) foot  ///< Footprint to measure
                                     ) {
    typedef typename ExposureT::MaskedImageT MaskedImageT;
    
    MaskedImageT const& mimage = getExposure()->getMaskedImage();
    if (foot) {
        src->setFootprint(foot);
    } else {
        foot = src->getFootprint();
    }

    bool const isNegative = (src->getFlagForDetection() & Flags::DETECT_NEGATIVE);
    //
    // Measure some properties of the Footprint
    //
    FootprintCentroid<MaskedImageT> centroidFunctor(mimage);
    centroidFunctor.apply(*foot);

    PTR(afwDetection::Peak) peak = centroidFunctor.makePeak(isNegative);
    //
    // Check for bits set in the Footprint
    //
    if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("EDGE")) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::EDGE);
    }
    if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::INTERP);
    }
    if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::SATUR);
    }
    //
    // Now run measure objects code (but not for edge objects)
    //
    typename MaskedImageT::Mask &mask = *mimage.getMask();
    if (mask(peak->getIx() - mask.getX0(), peak->getIy() - mask.getY0(),
             MaskedImageT::Mask::getMaskPlane("EDGE"))) {
        src->setFlagForDetection(src->getFlagForDetection() | Flags::EDGE);

        if (getMeasureAstrom()) {
            src->setAstrometry(getMeasureAstrom()->measure(getLog()));
        }
        if (getMeasureShape()) {
            src->setShape(getMeasureShape()->measure(getLog()));
        }
        if (getMeasurePhotom()) {
            src->setPhotometry(getMeasurePhotom()->measure(getLog()));
        }

        return;
    }
    //
    // Centroids
    //
    if (!getMeasureAstrom()) {
        src->setXAstrom(peak->getIx());
        src->setYAstrom(peak->getIy());
        src->setFlagForDetection(src->getFlagForDetection() | Flags::PEAKCENTER);
    } else {
        PTR(afwDetection::Measurement<afwDetection::Astrometry>) centroids =
            getMeasureAstrom()->measure(peak, src, getLog());
        src->setAstrometry(centroids);
        /*
         * Pack the answers into the Source
         */
        if (_policy.isString("source.astrom")) {
            std::string const& val = _policy.getString("source.astrom");
            if (val != "NONE") {
                afwDetection::Measurement<afwDetection::Astrometry>::TPtr astrom = centroids->find(val);
                
                double x = astrom->getX();
                double xErr = astrom->getXErr();
                double y = astrom->getY();
                double yErr = astrom->getYErr();

                if (lsst::utils::isnan(x) || lsst::utils::isnan(y)) {
                    // Everyone uses XAstrom and YAstrom, so we'll fix them up
                    src->setXAstrom(peak->getIx());
                    src->setYAstrom(peak->getIy());
                    src->setFlagForDetection(src->getFlagForDetection() | Flags::PEAKCENTER);
                } else {
                    peak->setFx(x);
                    peak->setFy(y);
                    
                    src->setXAstrom(x);
                    src->setYAstrom(y);
                    src->setXAstromErr(astrom->getXErr());
                    src->setYAstromErr(astrom->getYErr());
                }
                {                       // check if peak is off the image
                    int const ix = x;
                    int const iy = y;
                    if (ix < 0 || ix >= mimage.getWidth() || iy < 0 || iy >= mimage.getHeight()) {
                        x = peak->getIx();
                        y = peak->getIy();
                        // Fixup [xy]Err too?
                        src->setFlagForDetection(src->getFlagForDetection() | Flags::PEAKCENTER);
                    }
                }

                peak->setFx(x);
                peak->setFy(y);
                
                src->setXAstrom(x);
                src->setYAstrom(y);
                src->setXAstromErr(xErr);
                src->setYAstromErr(yErr);
            }
        }
    }
    //
    // Shapes
    //
    if (!getMeasureShape()) {
        ;
    } else {
        PTR(afwDetection::Measurement<afwDetection::Shape>) shapes =
            getMeasureShape()->measure(peak, src, getLog());
        src->setShape(shapes);
        /*
         * Pack the answers into the Source
         */
        if (_policy.isString("source.shape")) {
            std::string const& val = _policy.getString("source.shape");
            if (val != "NONE") {
                afwDetection::Measurement<afwDetection::Shape>::TPtr shape = shapes->find(val);

                src->setIxx(shape->getIxx());       // <xx>
                src->setIxxErr(shape->getIxxErr()); // sqrt(Var<xx>)
                src->setIxy(shape->getIxy());       // <xy>
                src->setIxyErr(shape->getIxyErr()); // sign(Covar(x, y))*sqrt(|Covar(x, y)|))        
                src->setIyy(shape->getIyy());       // <yy>
                src->setIyyErr(shape->getIyyErr()); // sqrt(Var<yy>)

                src->setPsfIxx(shape->getPsfIxx());       // <xx>
                src->setPsfIxxErr(shape->getPsfIxxErr()); // sqrt(Var<xx>)
                src->setPsfIxy(shape->getPsfIxy());       // <xy>
                src->setPsfIxyErr(shape->getPsfIxyErr()); // sign(Covar(x, y))*sqrt(|Covar(x, y)|))        
                src->setPsfIyy(shape->getPsfIyy());       // <yy>
                src->setPsfIyyErr(shape->getPsfIyyErr()); // sqrt(Var<yy>)

                src->setE1(shape->getE1());
                src->setE1Err(shape->getE1Err());
                src->setE2(shape->getE2());
                src->setE2Err(shape->getE2Err());
                src->setShear1(shape->getShear1());
                src->setShear1Err(shape->getShear1Err());
                src->setShear2(shape->getShear2());
                src->setShear2Err(shape->getShear2Err());

                src->setResolution(shape->getResolution());
                src->setShapeStatus(shape->getShapeStatus());
                src->setSigma(shape->getSigma());
                src->setSigmaErr(shape->getSigmaErr());
            }
        }
    }

    //
    // Photometry
    //
    if (!getMeasurePhotom()) {
        ;
    } else {
        PTR(afwDetection::Measurement<afwDetection::Photometry>) fluxes =
            getMeasurePhotom()->measure(peak, src, getLog());
        src->setPhotometry(fluxes);

        /*
         * Pack flags into the source
         */
        boost::int64_t flag = src->getFlagForDetection();
        for(afwDetection::Measurement<afwDetection::Photometry>::const_iterator i= fluxes->begin();
            i != fluxes->end(); ++i
        ) {
            flag |= (*i)->getFlag();
        }
        src->setFlagForDetection(flag);
        /*
         * Pack the answers into the Source
         */
        if (_policy.isString("source.apFlux")) {
            std::string const& val = _policy.getString("source.apFlux");
            if (val != "NONE") {
                afwDetection::Measurement<afwDetection::Photometry>::TPtr photom = fluxes->find(val);
                    
                src->setApFlux(photom->getFlux());
                src->setApFluxErr(photom->getFluxErr());
            }
        }
        
        if (_policy.isString("source.psfFlux")) {
            std::string const& val = _policy.getString("source.psfFlux");
            if (val != "NONE") {
                afwDetection::Measurement<afwDetection::Photometry>::TPtr photom = fluxes->find(val);
                    
                src->setPsfFlux(photom->getFlux());
                src->setPsfFluxErr(photom->getFluxErr());
            }
        }

        if (_policy.isString("source.modelFlux")) {
            std::string const& val = _policy.getString("source.modelFlux");
            if (val != "NONE") {
                afwDetection::Measurement<afwDetection::Photometry>::TPtr photom = fluxes->find(val);
                    
                src->setModelFlux(photom->getFlux());
                src->setModelFluxErr(photom->getFluxErr());
            }
        }

        if (_policy.isString("source.instFlux")) {
            std::string const& val = _policy.getString("source.instFlux");
            if (val != "NONE") {
                afwDetection::Measurement<afwDetection::Photometry>::TPtr photom = fluxes->find(val);
                    
                src->setInstFlux(photom->getFlux());
                src->setInstFluxErr(photom->getFluxErr());
            }
        }
    }
    
    //
    // Check for bits set near the centroid
    //
    {
        afwGeom::Point2I llc(afwImage::positionToIndex(src->getXAstrom()) - 1,
                             afwImage::positionToIndex(src->getYAstrom()) - 1);
        afwDetection::Footprint const middle(afwGeom::BoxI(llc, afwGeom::ExtentI(3))); // 3x3 centred at the the centroid
        centroidFunctor.apply(middle);
        if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("INTRP")) {
            src->setFlagForDetection(src->getFlagForDetection() | Flags::INTERP_CENTER);
        }
        if (centroidFunctor.getBits() & MaskedImageT::Mask::getPlaneBitMask("SAT")) {
            src->setFlagForDetection(src->getFlagForDetection() | Flags::SATUR_CENTER);
        }
    }

    //
    // Add some star/galaxy information.  The "extendedness" parameter is supposed to be the
    // probability of being extended
    //
    std::vector<float> fac(3);// Fiddle factors for star/galaxy separation
    fac[0] = getNumeric(_policy, "classification.sg_fac1");
    fac[1] = getNumeric(_policy, "classification.sg_fac2");
    fac[2] = getNumeric(_policy, "classification.sg_fac3");

    bool const isStar = ((fac[0]*src->getInstFlux() + fac[1]*src->getInstFluxErr()) <
                         (src->getPsfFlux() + fac[2]*src->getPsfFluxErr()) ? 0.0 : 1.0);

#if 0
    src->setExtendedness(isStar ? 0.0 : 1.0);
#else
    src->setApDia(isStar ? 0.0 : 1.0);
#endif
}

//
// Explicit instantiations
//
// \cond
template class MeasureSources<afwImage::Exposure<float> >;
template class MeasureSources<afwImage::Exposure<int> >;
// \endcond
}}}
