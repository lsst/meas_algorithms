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
 
#if !defined(LSST_MEAS_ALGORITHMS_PHOTOMETRY_H)
#define LSST_MEAS_ALGORITHMS_PHOTOMETRY_H 1
/**
 * @file
 */
#include <cmath>
#include <utility>
#include <string>
#include <typeinfo>
#include "boost/shared_ptr.hpp"
#include "boost/noncopyable.hpp"
#include "lsst/afw/image.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/meas/algorithms/detail/MeasureFactory.h"

namespace lsst {
namespace meas {
namespace algorithms {
/**
 * @brief Represent a position and its error
 */
class Photometry {
public:
    typedef boost::shared_ptr<Photometry> Ptr;
    typedef boost::shared_ptr<Photometry const> ConstPtr;
    typedef std::pair<double, double> xyAndError;

    Photometry(double apFlux = NAN, double psfFlux = NAN) :
        _apFlux(apFlux),
        _psfFlux(psfFlux) {}
    
    virtual ~Photometry() {}

    void setApFlux(double apFlux) { _apFlux = apFlux; }
    void setApFluxErr(double apFluxErr) { _apFluxErr = apFluxErr; }
    double getApFlux() const { return _apFlux; }
    double getApFluxErr() const { return _apFluxErr; }

    void setPsfFlux(double psfFlux) { _psfFlux = psfFlux; }
    void setPsfFluxErr(double psfFluxErr) { _psfFluxErr = psfFluxErr; }
    double getPsfFlux() const { return _psfFlux; }
    double getPsfFluxErr() const { return _psfFluxErr; }

private:

    double _apFlux, _apFluxErr;
    double _psfFlux, _psfFluxErr;

};

/************************************************************************************************************/
/**
 * @brief A pure virtual base class to calculate a photometry
 *
 * Different implementations will use different algorithms
 */
template<typename T>
class MeasurePhotometry : public MeasureProperty<MeasurePhotometry<T>, T> {
public:
    typedef T ImageT;
    
    typedef boost::shared_ptr<MeasurePhotometry> Ptr;
    typedef boost::shared_ptr<MeasurePhotometry const> ConstPtr;

    MeasurePhotometry(typename ImageT::ConstPtr image=typename ImageT::ConstPtr())
        : MeasureProperty<MeasurePhotometry<T>, T>(image) {}
    virtual ~MeasurePhotometry() {}

    Photometry apply(ImageT const& image, double xcen, double ycen,
                     lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                     double background = 0.0) const;
    Photometry apply(int x, int y,
                     lsst::meas::algorithms::PSF const* psf = NULL, // fully qualified to make swig happy
                     double background = 0.0) const {
        if (!this->getImage()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "You must provide an image to measure");
        }
        return apply(*this->getImage(), x, y, psf, background);
    }
    
    inline float getRadius() const {
        return _radius;
    }
    inline void setRadius(float const radius) const {
        _radius = radius;
    }

private:
    float mutable _radius;

    virtual Photometry doApply(ImageT const& image,
                               double xcen, double ycen, PSF const* psf, double background) const = 0;
};

/************************************************************************************************************/
/**
 * Provide a convenient wrapper round createMeasureProperty
 */
template<typename ImageT>
MeasurePhotometry<ImageT> *
createMeasurePhotometry(
        std::string const& type,        ///< Algorithm type (e.g. "NAIVE")
        boost::shared_ptr<ImageT const> image=boost::shared_ptr<ImageT const>() ///< The image to process
                       )
{
    MeasurePhotometry<ImageT> const* ptr = NULL;
    return createMeasureProperty(type, image, ptr);
}
    
}}}
#endif
