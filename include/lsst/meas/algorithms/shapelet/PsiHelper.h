#ifndef MeasAlgoShapeletPsiHelper_H
#define MeasAlgoShapeletPsiHelper_H

#include "lsst/meas/algorithms/shapelet/MyMatrix.h"

// Here is the order of p,q along the indices of psi:
//
// k 0 1,2 3,4 5 6,7 8,9 10,11 12,13 14 15,16 17,18 19,20 21,22 23,24 25,26 27
// p 0  1   2  1  3   2    4     3    2   5     4     3     6     5     4    3
// q 0  0   0  1  0   1    0     1    2   0     1     2     0     1     2    3
// n 0  1   2  2  3   3    4     4    4   5     5     5     6     6     6    6
// m 0  1   2  0  3   1    4     2    0   5     3     1     6     4     2    0

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    // Make the psi matrix.
    // Each row is a pixel.
    // The columns correspond to each psi index (see above).
    // The input z vector has the position of each pixel.
    // The optional coeff parameter specifies a coefficient to multiply each row.
    void makePsi(
        DMatrix& psi, CDVectorView z, int order, const DVector* coeff=0);

    // Same thing, but for a single pixel.
    void makePsi(DVector& psi, std::complex<double> z, int order);

    // Add some more columns to a psi matrix for two more radial orders.
    // The order is the original order of the matrix.
    // So the result is a psi matrix for order+2.
    void augmentPsi(DMatrix& psi, CDVectorView z, int order);

    // Gx = d(psi)/dx.  Likewise for the other parameters.
    void setupGx(DMatrix& Gx, int order1, int order2);
    void setupGy(DMatrix& Gy, int order1, int order2);
    void setupGg1(DMatrix& Gg1, int order1, int order2);
    void setupGg2(DMatrix& Gg2, int order1, int order2);
    void setupGmu(DMatrix& Gmu, int order1, int order2);
    void setupGth(DMatrix& Gth, int order1, int order2);

}}}}

#endif
