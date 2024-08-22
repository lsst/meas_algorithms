#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace std;

#define DSQUARE(x) (double(x)*double(x))
#define DEGPRAD (180.0L/M_PI) /*Degrees per radian*/
#define LARGERR 1e30L // Large number supposed to be a safe initialization
                      // for most minimum-finding problems.

class longpair{ // Pair of long integers
public:
  long i1;
  long i2;
  longpair(long i1, long i2) :i1(i1), i2(i2) { }
  longpair() = default;
};

class long_index{   // Pairs a long integer with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a long integer element,
                    // without actually re-sorting the vector.
public:
  long lelem;
  long index;
  long_index(long lelem, long index) :lelem(lelem), index(index) {}
};

class lower_long_index{
public:
  inline bool operator() (const long_index& i1, const long_index& i2) {
    return(i1.lelem < i2.lelem);
  }
};

class double_index{ // Pairs a double with an index, intended for use
                    // accessing elements in a vector of a more complex
                    // class in an order sorted by a double-precision component,
                    // without actually re-sorting the vector.
public:
  double delem;
  long index;
  double_index(double delem, long index) :delem(delem), index(index) {}
};

class lower_double_index{
public:
  inline bool operator() (const double_index& i1, const double_index& i2) {
    return(i1.delem < i2.delem);
  }
};

class xy_index{ // Double-precision x,y point plus long index
public:
  double x;
  double y;
  long index;
  xy_index(double x, double y, long index) :x(x), y(y), index(index) { }
  xy_index() = default;
};
  
class xyind_lower_x{ // Function for sorting vectors of type xy_index by x
public:
  inline bool operator() (const xy_index& p1, const xy_index& p2) {
    return(p1.x < p2.x);
  }
};

class xyind_lower_y{ // Function for sorting vectors of type xy_index by y
public:
  inline bool operator() (const xy_index& p1, const xy_index& p2) {
    return(p1.y < p2.y);
  }
};

class kdpoint { // Point in an xy_index k-d tree designed as a vector,
                // with left and right giving the index of the left and
                // right branches, and the dimension of splitting
                // specified at each node.
public:
  xy_index point;
  long left;
  long right;
  int dim;
  kdpoint(xy_index point, long left, long right, int dim) :point(point), left(left), right(right), dim(dim) {}
  kdpoint() = default;
};

class point3d_index{ // Double-precision 3-D point with long-integer idex
public:
  double x;
  double y;
  double z;
  long index;
  point3d_index(double x, double y, double z, long index) :x(x), y(y), z(z), index(index) { }
  point3d_index() = default;
};

class lower_point3d_index_x{ // Sort point3d_index by x
public:
  inline bool operator() (const point3d_index& p1, const point3d_index& p2) {
    return(p1.x < p2.x);
  }
};

class lower_point3d_index_y{ // Sort point3d_index by y
public:
  inline bool operator() (const point3d_index& p1, const point3d_index& p2) {
    return(p1.y < p2.y);
  }
};

class lower_point3d_index_z{ // Sort point3d_index by z
public:
  inline bool operator() (const point3d_index& p1, const point3d_index& p2) {
    return(p1.z < p2.z);
  }
};

struct FindGlintsConfig {       // Note that the pixel-xy and RA/Dec versions use the same config struct.
  long minpoints = 5;          // Minimum number of points that counts as a glint-trail
  double maxgcr = 1.0;         // Maximum RMS deviation from a straight line (pixels or arcsec)
  double maxrange = 500.0;     // Range of the k-d tree query (pixels or arcsec) used to identify groups of sources
                               // that could make up a glint-trail.
  // The next three parameters are relevant only for the RA/Dec version.
  int centerknown = 0;         // Were RA and Dec of the image center supplied by the calling function?
  double incenRA = 0.0;        // RA of image center, supplied by calling function.
  double incenDec = 0.0;       // Dec of image center, supplied by calling function.
  // The next three parameters apply to the check for equal spacing of glints along the trail.
  // The reference frequency is set by the number of detected sources (npt) divided by the
  // length of the trail (traillen).
  double freq_downscale = 0.9; // Lowest frequency probed, in units of the reference frequency.
  double freq_upscale = 3.0;   // Highest frequency probed, in units of the reference frequency.
                               // The default setting of 3.0 enables scanning of the correct
                               // frequency even if only 1/3 of the flashes were detected.
  double max_phase_err = 0.02; // Determines the sampling interval for the frequency search, which is
                               // set to max_phase_err/traillen
};

class glint_trail{ // Trail of glints within a single image
public:
  double x;        // RA in decimal degrees, or pixel x coodinate
  double y;        // Dec in decimal degrees, or pixel y coordinate
  double length;   // trail length in arcsec or pixels
  double PA;       // position angle, east from north or CCW from +y, degrees
  double linrms;   // RMS deviation from best-fit line, arcsec or pixels
  double eqrms;    // RMS deviation from best line with constant step size, arcsec or pixels
  double magmean;  // Mean magnitude of individual flashes
  double magrms;   // RMS scatter of magnitudes about mean
  double stepsize; // Distance between adjacent flashes, arcsec or pixels
  double qc1;      // Additional parameter reserved for flexibility
  long npt;        // Number of detections
  long flashnum;    // Inferred total number of flashes, allowing for missing detections
  glint_trail(double x, double y, double length, double PA, double linrms, double eqrms, double magmean, double magrms, double stepsize, double qc1, long npt, long flashnum) :x(x), y(y), length(length), PA(PA), linrms(linrms), eqrms(eqrms), magmean(magmean), magrms(magrms), stepsize(stepsize), qc1(qc1), npt(npt), flashnum(flashnum) { }
  glint_trail() = default;
};



void make_ivec(long nx, vector <int> &ivec);
int dmeanrms01(const vector <double> &invec, double *mean, double *rms);
int linfituw01(const vector <double> &x, const vector <double> &y, double &slope, double &intercept);
long medindex(const vector <xy_index> &xyvec, int dim);
int splitxy(const vector <xy_index> &xyvec, int dim, long unsigned int splitpoint, vector <xy_index> &left, vector <xy_index> &right);
int kdtree01(const vector <xy_index> &xyvec, int dim, long unsigned int rootptxy, long unsigned int rootptkd, vector <kdpoint> &kdvec);
int kdrange01(const vector <kdpoint> &kdvec,double x,double y,double range,vector <long> &indexvec);
int find_glints_xypix(const vector <point3d_index> &detvec, FindGlintsConfig config, vector <glint_trail> &trailvec, vector <longpair> &trail2det);
