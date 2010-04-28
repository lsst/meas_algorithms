#ifndef MeasAlgoShapeletMyMatrix_H
#define MeasAlgoShapeletMyMatrix_H

#ifdef USE_TMV

#include "TMV.h"
#include "TMV_Sym.h"
#include "TMV_Small.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {


    typedef tmv::Vector<double> DVector;
    typedef tmv::VectorView<double> DVectorView;
    typedef tmv::ConstVectorView<double> DConstVectorView;
    typedef tmv::Matrix<double> DMatrix;
    typedef tmv::MatrixView<double> DMatrixView;
    typedef tmv::ConstMatrixView<double> DConstMatrixView;
    typedef tmv::SmallMatrix<double,2,2> DSmallMatrix22;
    typedef tmv::Matrix<double> DRowVector;

    typedef tmv::Vector<std::complex<double> > CDVector;
    typedef tmv::VectorView<std::complex<double> > CDVectorView;
    typedef tmv::ConstVectorView<std::complex<double> > CDConstVectorView;
    typedef tmv::Matrix<std::complex<double> > CDMatrix;
    typedef tmv::MatrixView<std::complex<double> > CDMatrixView;
    typedef tmv::ConstMatrixView<std::complex<double> > CDConstMatrixView;

    typedef tmv::Vector<float> FVector;
    typedef tmv::VectorView<float> FVectorView;
    typedef tmv::ConstVectorView<float> FConstVectorView;
    typedef tmv::Matrix<float> FMatrix;
    typedef tmv::MatrixView<float> FMatrixView;
    typedef tmv::ConstMatrixView<float> FConstMatrixView;
    typedef tmv::SmallMatrix<float,2,2> FSmallMatrix22;

    typedef tmv::Vector<std::complex<float> > CFVector;
    typedef tmv::VectorView<std::complex<float> > CFVectorView;
    typedef tmv::ConstVectorView<std::complex<float> > CFConstVectorView;
    typedef tmv::Matrix<std::complex<float> > CFMatrix;
    typedef tmv::MatrixView<std::complex<float> > CFMatrixView;
    typedef tmv::ConstMatrixView<std::complex<float> > CFConstMatrixView;

    typedef tmv::DiagMatrix<double> DDiagMatrix;
    typedef tmv::SymMatrix<double> DSymMatrix;

}}}}

#define TVector(T) tmv::Vector<T>
#define TVectorView(T) tmv::VectorView<T>
#define TConstVectorView(T) tmv::ConstVectorView<T>
#define TMatrix(T) tmv::Matrix<T>
#define TMatrixView(T) tmv::MatrixView<T>
#define TConstMatrixView(T) tmv::ConstMatrixView<T>

// Macros that go after a . or ->
#define TMV_rowsize() rowsize()
#define TMV_colsize() colsize()
#define TMV_cref(i,j) cref(i,j)
#define TMV_ref(i,j) ref(i,j)
#define TMV_diag() diag()
#define EIGEN_diag() 
#define EIGEN_asDiag()
#define TMV_setToIdentity() setToIdentity()
#define TMV_subVector(i1,i2) subVector(i1,i2)
#define TMV_subMatrix(i1,i2,j1,j2) subMatrix(i1,i2,j1,j2)
#define TMV_realPart() realPart()
#define TMV_imagPart() imagPart()
#define TMV_normSq() normSq()
#define TMV_normInf() normInf()
#define TMV_maxAbsElement() maxAbsElement()
#define TMV_addToAll(x) addToAll(x)
#define TMV_setAllTo(x) setAllTo(x)
#define TMV_colpart(j,i1,i2) col(j,i1,i2)
#define TMV_rowpart(i,j1,j2) row(i,j1,j2)
#define TMV_det() det()
#define TMV_transposeSelf() transposeSelf()
#define TMV_sumElements() sumElements()

// Macros that require the matrix or vector as an argument:
#define TMV_view(m) (m).view()
#define TMV_vview(v) (v).view()
#define TMV_colRange(m,j1,j2) (m).colRange(j1,j2)
#define TMV_rowRange(m,i1,i2) (m).rowRange(i1,i2)
#define TMV_ptr(m) (m).ptr()
#define TMV_cptr(m) (m).cptr()
#define TMV_stepi(m) (m).stepi()
#define TMV_stepj(m) (m).stepj()
#define TMV_DiagMatrixViewOf(v) DiagMatrixViewOf(v)
#define TMV_conjugateSelf(m) (m).conjugateSelf();
#define EIGEN_Transpose(m) m
#define EIGEN_ToScalar(m) m

// Standalone macros:
#define TMV_const const
#define EIGEN_mutable 
#define EIGEN_twice(x) x

// QR solver:
#define TMV_QR(m) (m).divideUsing(tmv::QR); (m).saveDiv()
#define TMV_QR1(m) (m).divideUsing(tmv::QR); (m).saveDiv()
#define TMV_QR2(m) (m).resetDiv()
#define TMV_QRisSingular(m) (m).isSingular()
#define TMV_throwSingular throw tmv::Singular()
#define TMV_QR_Solve(m,x,b) x = b/m
#define TMV_QR_InverseATA(m,cov) (m).makeInverseATA(cov)

#else

#include "Eigen/Core"
#include "Eigen/QR"
#include "Eigen/SVD"
#include "Eigen/LU"
#include "Eigen/Cholesky"
#include "Eigen/Array"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {


    typedef Eigen::VectorXd DVector;
    typedef Eigen::Block<DVector,Eigen::Dynamic,1> DVectorView;
    typedef const DVectorView DConstVectorView;
    typedef Eigen::MatrixXd DMatrix;
    typedef Eigen::Block<DMatrix> DMatrixView;
    typedef const DMatrixView DConstMatrixView;
    typedef Eigen::Matrix<double,2,2> DSmallMatrix22;
    typedef Eigen::RowVectorXd DRowVector;

    typedef Eigen::VectorXcd CDVector;
    typedef Eigen::Block<CDVector,Eigen::Dynamic,1> CDVectorView;
    typedef const CDVectorView CDConstVectorView;
    typedef Eigen::MatrixXcd CDMatrix;
    typedef Eigen::Block<CDMatrix> CDMatrixView;
    typedef const CDMatrixView CDConstMatrixView;

    typedef Eigen::VectorXf FVector;
    typedef Eigen::Block<FVector,Eigen::Dynamic,1> FVectorView;
    typedef const FVectorView FConstVectorView;
    typedef Eigen::MatrixXf FMatrix;
    typedef Eigen::Block<FMatrix> FMatrixView;
    typedef const FMatrixView FConstMatrixView;
    typedef Eigen::Matrix<float,2,2> FSmallMatrix22;

    typedef Eigen::VectorXcf CFVector;
    typedef Eigen::Block<CFVector,Eigen::Dynamic,1> CFVectorView;
    typedef const CFVectorView CFConstVectorView;
    typedef Eigen::MatrixXcf CFMatrix;
    typedef Eigen::Block<CFMatrix> CFMatrixView;
    typedef const CFMatrixView CFConstMatrixView;

    typedef DVector DDiagMatrix;
    typedef DMatrix DSymMatrix;

}}}}

#define TVector(T) Eigen::Matrix<T,Eigen::Dynamic,1>
#define TVectorView(T) Eigen::Block<TVector(T) >
#define TConstVectorView(T) const Eigen::Block<TVector(T) >
#define TMatrix(T) Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>
#define TMatrixView(T) Eigen::Block<TMatrix(T) >
#define TConstMatrixView(T) const Eigen::Block<TMatrix(T) >

// Macros that go after a . or ->
#define TMV_rowsize() cols()
#define TMV_colsize() rows()
#define TMV_cref(i,j) coeff(i,j)
#define TMV_ref(i,j) coeffRef(i,j)
#define TMV_diag() diagonal()
#define EIGEN_diag() diagonal().cwise()
#define EIGEN_asDiag() .asDiagonal()
#define TMV_setToIdentity() setIdentity()
#define TMV_subVector(i1,i2) segment(i1,((i2)-(i1)))
#define TMV_subMatrix(i1,i2,j1,j2) block(i1,j1,((i2)-(i1)),((j2)-(j1)))
#define TMV_realPart() real()
#define TMV_imagPart() imag()
#define TMV_normSq() squaredNorm()
#define TMV_normInf() lpNorm<Eigen::Infinity>()
#define TMV_maxAbsElement() cwise().abs().maxCoeff()
#define TMV_addToAll(x) cwise() += (x)
#define TMV_setAllTo(x) setConstant(x)
#define TMV_colpart(j,i1,i2) col(j).segment(i1,i2)
#define TMV_rowpart(i,j1,j2) row(i).segment(j1,j2)
#define TMV_det() determinant()
#define TMV_transposeSelf() transposeInPlace()
#define TMV_sumElements() sum()

// Macros that require the matrix or vector as an argument:
#define TMV_view(m) (m).block(0,0,(m).rows(),(m).cols())
#define TMV_vview(v) (v).segment(0,(v).size())
#define TMV_colRange(m,j1,j2) (m).block(0,j1,(m).rows(),((j2)-(j1)))
#define TMV_rowRange(m,i1,i2) (m).block(i1,0,((i2)-(i1)),(m).cols())
#define TMV_ptr(m) &((m).coeffRef(0,0))
#define TMV_cptr(m) &((m).coeffRef(0,0))
#define TMV_stepi(m) ((m).Flags & Eigen::RowMajorBit ? (m).stride() : 1)
#define TMV_stepj(m) ((m).Flags & Eigen::RowMajorBit ? 1 :  (m).stride())
#define TMV_DiagMatrixViewOf(v) (v).asDiagonal()
#define TMV_conjugateSelf(m) (m) = (m).conjugate()
#define EIGEN_Transpose(m) (m).transpose()
#define EIGEN_ToScalar(m) (m)(0,0)

// Standalone macros:
#define TMV_const 
#define EIGEN_mutable mutable
#define EIGEN_twice(x) x,x

// QR solver:
#define TMV_QR(m) Eigen::QR<DMatrix> QR_Solver_ ## m  = (m).qr();
#define TMV_QR1(m)
#define TMV_QR2(m) Eigen::QR<DMatrix> QR_Solver_ ## m  = (m).qr();
// Eigen changed how the QR module works sometime between 2.0.0 and 2.0.10
#if (EIGEN_MINOR_VERSION >= 10) || (EIGEN_MAJOR_VERSION > 0)
#define TMV_QRisSingular(m) QR_Solver_ ## m .isInjective()
#define TMV_QR_Solve(m,x,b) QR_Solver_ ## m .solve(b,&(x));
#else
#define TMV_QRisSingular(m) QR_Solver_ ## m .isFullRank()
// 2.0.0 doesn't even have a QR.solve() function.
#define TMV_QR_Solve(m,x,b) \
    (x) = QR_Solver_ ## m .matrixQ().transpose() * (b); \
    QR_Solver_ ## m .matrixR().solveTriangularInPlace(x)
#endif
#define TMV_throwSingular throw std::runtime_error("Singular")
#define TMV_QR_InverseATA(m,cov) \
    do { \
        (cov).setIdentity(); \
        QR_Solver_ ## m .matrixR().transpose().solveTriangularInPlace(cov); \
        QR_Solver_ ## m .matrixR().solveTriangularInPlace(cov); \
    } while (false)

#endif

#endif
