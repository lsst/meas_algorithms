#include "lsst/meas/algorithms/Ellipse.h"

namespace algorithms = lsst::meas::algorithms;

algorithms::Ellipse::Ellipse(double x, double y, double Ixx, double Iyy, double Ixy)
	: _x(x), _y(y)
{
    double halfXxMinusYy = 0.5 * (Ixx - Iyy);
    double halfXxPlusYy  = 0.5 * (Ixx + Iyy);
            
	_theta = 0.5 * std::atan2(Ixy, halfXxMinusYy);       
            
    double tmp = (halfXxMinusYy*halfXxMinusYy) + Ixy*Ixy;
    tmp = std::sqrt(tmp);        
    
    _majorAxis = std::sqrt(halfXxPlusYy + tmp);
    _minorAxis = std::sqrt(halfXxPlusYy - tmp);
} 


bool algorithms::Ellipse::contains(double x, double y)
{
	//Find the two foci of the ellipse: f1 and f2
	double f1_x, f1_y;
	double f2_x, f2_y;	
	
	//c is the distance from the origin to the Foci along the major axis
	double c = std::sqrt(a2() - b2());
	
	double cos = std::cos(_theta)*c;
	double sin = std::sin(_theta)*c;
	
	//f1 is the positive foci
	f1_x = _x + cos;
	f1_y = _y + sin;
	
	//f2 is the negative foci
	f2_x = _x - cos;
	f2_y = _y - sin;
	
	//calculate the distance to each foci from point (x,y)
	double dx, dy;
	
	dx = x - f1_x;
	dy = y - f1_y;	
	double distanceToF1 = std::sqrt(dx*dx + dy*dy);
	
	dx = x - f2_x;
	dy = y - f2_y;
	double distanceToF2 = std::sqrt(dx*dx + dy*dy);
	
	//calculate the sum of distance from foci to the curve
	//sum == r1 + r2 == sqrt(b^2 + c^2) == 2*a
	//where a is the major axis
	double distanceToCurve = 2*_majorAxis;
	
	return (distanceToF2+distanceToF1) <= distanceToCurve;
}
