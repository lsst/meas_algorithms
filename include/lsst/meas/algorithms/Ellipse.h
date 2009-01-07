#ifndef LSST_MEAS_ALGORTIHMS_ELLIPSE_H
#define LSST_MEAS_ALGORTIHMS_ELLIPSE_H

#include <cmath>

namespace lsst {
namespace meas {
namespace algorithms {

class Ellipse {
public:
    /**
     * \brief Construct from second momemnts
     */
    Ellipse(double x, double y, double Ixx, double Iyy, double Ixy);  
    
    
    inline double getTheta() const     {return _theta; }
    inline double getMajorAxis() const {return _majorAxis; }
    inline double getMinorAxis() const {return _minorAxis; }
    
    /// Get axis ratio.
    inline double getAxisRatio() const { return _minorAxis/_majorAxis; }

    /// Get geometric mean radius.
    inline double getMeanRadius() const { return std::sqrt(_majorAxis*_minorAxis); }
    
    /// The two-dimensional size of the ellipse.
    inline double getArea() const { return _majorAxis*_minorAxis*M_PI; }
    
    /// The eccentricity of the ellipse
    double getEccentricity() const { return std::sqrt(1 - b2()/a2());}
    
    /**
     * @brief Test the inclusion of point (x,y)
     *
     * @return true if point lies on or inside the ellipse, false otherwise
     */
    bool contains(double x, double y);
    
private:
	double _x;
	double _y;
    double _theta;
    double _majorAxis;
    double _minorAxis;    
    
    double a2() const {return _majorAxis*_majorAxis;}
    double b2() const {return _minorAxis*_minorAxis;}
};

}}} //namespace lsst::meas::algorithms

#endif //LSST_MEAS_ALGORTIHMS_FWHM_H
