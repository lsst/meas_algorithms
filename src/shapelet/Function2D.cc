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
#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>
#include <stdexcept>

#include "lsst/meas/algorithms/shapelet/Function2D.h"
#include "lsst/meas/algorithms/shapelet/dbg.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    Constant2D::Constant2D(std::istream& fin) : Function2D()
    {
        if(!(fin >> (*_coeffs)(0,0))) throw std::runtime_error("reading constant");
    }

    void Constant2D::write(std::ostream& fout) const
    {
        int oldPrec = fout.precision(6);
        std::ios_base::fmtflags oldf =
            fout.setf(std::ios_base::scientific,std::ios_base::floatfield);
        fout << "C " << (*_coeffs)(0,0) << std::endl;
        if (!fout) throw std::runtime_error("writing (constant) function");
        fout.precision(oldPrec);
        fout.flags(oldf);
    }

    void Constant2D::operator+=(const Function2D& rhs)
    {
        const Constant2D* crhs = dynamic_cast<const Constant2D*>(&rhs);
        Assert(crhs);
        (*_coeffs)(0,0) += (*crhs->_coeffs)(0,0);
    }


    void Polynomial2D::setFunction(int xo, int yo, const DVector& fVect)
    {
        if(_xOrder != xo || _yOrder != yo) {
            _xOrder = xo; _yOrder = yo;
            _coeffs.reset(new DMatrix(xo+1,yo+1));
            _coeffs->setZero();
        }
        int k=0;

        int maxOrder = std::max(xo,yo);
        for(int m=0;m<=maxOrder;++m) {
            int i0 = std::min(xo,m);
            int len = std::min(yo,i0)+1;
#ifdef USE_TMV
            DVectorView coeffDiag = _coeffs->subVector(i0,m-i0,-1,1,len);
            tmv::ConstVectorView<double> subf = fVect.subVector(k,k+len);
            coeffDiag = subf;
#else
            for(int i=0;i<len;++i)
                (*_coeffs)(i0-i,m-i0+i) = fVect(k+i);
#endif
            k += len;
        }

        Assert(k==(int)fVect.size());
    }

    Polynomial2D::Polynomial2D(std::istream& fin) : 
        Function2D()
    {
        // Order of parameters:  (example is for xorder = 2, yorder = 3
        // xorder(2) yorder(3) a00 a10 a01 a20 a11 a02 a21 a12 a03
        // where f = a00 + a10 x + a01 y + a20 x^2 + a11 xy + a02 y^2
        //           + a21 x^2 y + a12 xy^2 + a03 y^3
        // Note that total order doesn't go past the max of xorder and yorder.
        // Also note that a30 is not listed since xorder is only 2.
        // Note that aij are complex numbers so each is listed as 
        // real_part imag_part.
        int xo,yo;
        if (!(fin >> xo >> yo >> _scale)) 
            throw std::runtime_error("reading xorder,yorder,scale");
        _xOrder = xo;
        _yOrder = yo;
        _coeffs.reset(new DMatrix(xo+1,yo+1));
        _coeffs->setZero();
        int maxOrder = std::max(xo,yo);
        for(int m=0;m<=maxOrder;++m) {
            int i0 = std::min(xo,m);
            int len = std::min(yo,i0)+1;
#ifdef USE_TMV
            DVectorView coeffDiag = _coeffs->subVector(i0,m-i0,-1,1,len);
            for(int i=0;i<len;++i) fin >> coeffDiag(i);
#else
            for(int i=0;i<len;++i) fin >> (*_coeffs)(i0-i,m-i0+i);
#endif
        }
        if (!fin) throw std::runtime_error("reading (polynomial)");
    }

    void Polynomial2D::write(std::ostream& fout) const
    {
        int oldPrec = fout.precision(6);
        std::ios_base::fmtflags oldf = 
            fout.setf(std::ios_base::scientific,std::ios_base::floatfield);
        int maxOrder = std::max(_xOrder,_yOrder);
        if (maxOrder == 0) {
            fout << "C " << (*_coeffs)(0,0) << std::endl;
        } else {
            fout << "P " << _xOrder << ' ' << _yOrder << ' ' << _scale << ' ';
            for(int m=0;m<=maxOrder;++m) {
                int i0 = std::min(_xOrder,m);
                int len = std::min(_yOrder,i0)+1;
#ifdef USE_TMV
                DVectorView coeffDiag = _coeffs->subVector(i0,m-i0,-1,1,len);
                for(int i=0;i<len;++i) fout << coeffDiag(i) << ' ';
#else
                for(int i=0;i<len;++i) fout << (*_coeffs)(i0-i,m-i0+i);
#endif
            }
        }
        fout << std::endl;
        if (!fout) throw std::runtime_error("writing (polynomial) function");
        fout.flags(oldf);
        fout.precision(oldPrec);
    }

    void Polynomial2D::addLinear(double a, double b, double c)
    {
        (*_coeffs)(0,0) += a;
        (*_coeffs)(1,0) += b*_scale;
        (*_coeffs)(0,1) += c*_scale;
    }

    void Polynomial2D::linearPreTransform(
        double a, double b, double c, double d, double e, double f)
    {
        // F(x,y) = Sum_i,j a(i,j) x^i y^j
        // F'(x,y) = F(a+bx+cy,d+ex+fy)
        //         = Sum_i,j a(i,j) (a+bx+cy)^i (d+ex+fy)^j
        //         = Sum_i,j a(i,j) (Sum_kl iCk kCl a^i-k (bx)^k-l (cy)^l) *
        //                Sum_mn jCm mCn d^j-m (ex)^m-n (fy)^n
        int maxOrder = std::max(_xOrder,_yOrder);
        std::vector<double> scaleToThe(maxOrder+1);
        scaleToThe[0] = 1.0; scaleToThe[1] = _scale;
        for(int i=2;i<=maxOrder;++i) scaleToThe[i] = scaleToThe[i-1]*_scale;
        std::vector<double> aToThe(maxOrder+1);
        std::vector<double> bToThe(maxOrder+1);
        std::vector<double> cToThe(maxOrder+1);
        std::vector<double> dToThe(maxOrder+1);
        std::vector<double> eToThe(maxOrder+1);
        std::vector<double> fToThe(maxOrder+1);
        aToThe[0] = 1.; aToThe[1] = a;
        bToThe[0] = 1.; bToThe[1] = b;
        cToThe[0] = 1.; cToThe[1] = c;
        dToThe[0] = 1.; dToThe[1] = d;
        eToThe[0] = 1.; eToThe[1] = e;
        fToThe[0] = 1.; fToThe[1] = f;
        for(int i=2;i<=maxOrder;++i) aToThe[i] = aToThe[i-1]*a;
        for(int i=2;i<=maxOrder;++i) bToThe[i] = bToThe[i-1]*b;
        for(int i=2;i<=maxOrder;++i) cToThe[i] = cToThe[i-1]*c;
        for(int i=2;i<=maxOrder;++i) dToThe[i] = dToThe[i-1]*d;
        for(int i=2;i<=maxOrder;++i) eToThe[i] = eToThe[i-1]*e;
        for(int i=2;i<=maxOrder;++i) fToThe[i] = fToThe[i-1]*f;
        _xOrder = maxOrder; _yOrder = maxOrder;

        DMatrix binom(maxOrder+1,maxOrder+1);
        binom(0,0) = 1.0;
        for(int n=1;n<=maxOrder;++n) {
            binom(n,0) = 1.0;
            binom(n,n) = 1.0;
            for(int m=1;m<n;++m) {
                binom(n,m) = binom(n-1,m-1) + binom(n-1,m);
            }
        }

        std::auto_ptr<DMatrix> oldCoeffs = _coeffs;
        _coeffs.reset(new DMatrix(_xOrder+1,_yOrder+1));
        _coeffs->setZero();
        for(int i=0;i<=_xOrder;++i) for(int j=0;j<=_yOrder&&i+j<=maxOrder;++j) {
            for(int k=0;k<=i;++k) for(int l=0;l<=k;++l) {
                for(int m=0;m<=j;++m) for(int n=0;n<=m;++n) {
                    (*_coeffs)(k-l+m-n,l+n) += 
                        binom(i,k)*binom(k,l)*binom(j,m)*binom(m,n)*
                        aToThe[i-k]*bToThe[k-l]*cToThe[l]*
                        dToThe[j-m]*eToThe[m-n]*fToThe[n]*
                        (*oldCoeffs)(i,j)/scaleToThe[i+j-k-m];
                }
            }
        }
    }

    void Polynomial2D::operator+=(const Function2D& rhs)
    {
        const Polynomial2D* prhs = dynamic_cast<const Polynomial2D*>(&rhs);
        Assert(prhs);
        Assert(_scale == prhs->_scale);
        if (_xOrder == prhs->_xOrder && _yOrder == prhs->_yOrder) {
            *_coeffs += *prhs->_coeffs;
        } else {
            int newXOrder = std::max(_xOrder,prhs->_xOrder);
            int newYOrder = std::max(_yOrder,prhs->_yOrder);
            std::auto_ptr<DMatrix > newCoeffs(
                new DMatrix(newXOrder+1,newYOrder+1));
            newCoeffs->setZero();
            newCoeffs->TMV_subMatrix(0,_xOrder+1,0,_yOrder+1) = *_coeffs;
            newCoeffs->TMV_subMatrix(0,prhs->_xOrder+1,0,prhs->_yOrder+1) += 
                *prhs->_coeffs;
            _coeffs = newCoeffs;
            _xOrder = newXOrder;
            _yOrder = newYOrder;
        }
    }

    void Polynomial2D::makeProductOf(
        const Polynomial2D& f, const Polynomial2D& g)
    {
        // h(x,y) = f(x,y) * g(x,y)
        //        = (Sum_ij f_ij x^i y^j) (Sum_mn g_mn x^m y^n)
        //        = Sum_ijmj f_ij g_mn x^(i+m) y^(j+n)
        Assert(_scale == f._scale);
        Assert(_scale == g._scale);
        int newXOrder = f.getXOrder() + g.getXOrder();
        int newYOrder = f.getYOrder() + g.getYOrder();
        if (_xOrder != newXOrder || _yOrder != newYOrder) {
            _xOrder = newXOrder;
            _yOrder = newYOrder;
            _coeffs.reset(new DMatrix(_xOrder+1,_yOrder+1));
        }
        _coeffs->setZero();
        for(int i=0;i<=f.getXOrder();++i) for(int j=0;j<=f.getYOrder();++j) {
            for(int m=0;m<=g.getXOrder();++m) for(int n=0;n<=g.getYOrder();++n) {
                Assert (i+m <= _xOrder && j+n <= _yOrder);
                (*_coeffs)(i+m,j+n) += (*f._coeffs)(i,j) * (*g._coeffs)(m,n);
            }
        }
    }

    std::auto_ptr<Function2D> Polynomial2D::dFdX() const
    {
        if (_xOrder == 0) {
            return std::auto_ptr<Function2D>(new Constant2D());
        }
        if (_xOrder == 1 && _yOrder == 0) {
            return std::auto_ptr<Function2D>(
                new Constant2D((*_coeffs)(1,0)));
        }

        int newXOrder = _xOrder-1;
        int newYOrder = _xOrder > _yOrder ? _yOrder : _yOrder-1;

        std::auto_ptr<Polynomial2D> temp(
            new Polynomial2D(newXOrder,newYOrder));

        int maxOrder = std::max(newXOrder,newYOrder);
        for(int i=newXOrder;i>=0;--i) {
            for(int j=std::min(maxOrder-i,newYOrder);j>=0;--j) {
                Assert(i+1<=_xOrder);
                Assert(j<=_yOrder);
                Assert(i+1+j<=std::max(_xOrder,_yOrder));
                (*temp->_coeffs)(i,j) = (*_coeffs)(i+1,j)*(i+1.)/_scale;
            }
        }
        return std::auto_ptr<Function2D>(temp);
    }

    std::auto_ptr<Function2D> Polynomial2D::dFdY() const 
    {
        if (_yOrder == 0) {
            return std::auto_ptr<Function2D>(new Constant2D());
        }
        if (_yOrder == 1 && _xOrder == 0) {
            return std::auto_ptr<Function2D>(
                new Constant2D((*_coeffs)(0,1)));
        }

        int newXOrder = _yOrder > _xOrder ? _xOrder : _xOrder-1;
        int newYOrder = _yOrder-1;

        std::auto_ptr<Polynomial2D> temp(
            new Polynomial2D(newXOrder,newYOrder));

        int maxOrder = std::max(newXOrder,newYOrder);
        for(int i=newXOrder;i>=0;--i) 
            for(int j=std::min(maxOrder-i,newYOrder);j>=0;--j) {
                Assert(i<=_xOrder);
                Assert(j+1<=_yOrder);
                Assert(i+j+1<=std::max(_xOrder,_yOrder));
                (*temp->_coeffs)(i,j) = (*_coeffs)(i,j+1)*(j+1.)/_scale;
            }
        return std::auto_ptr<Function2D>(temp);
    }

    std::auto_ptr<Function2D> Function2D::conj() const
    {
        std::auto_ptr<Function2D> temp = copy();
        TMV_conjugateSelf(*(temp->_coeffs));
        return temp;
    }

    double Function2D::operator()(double x,double y) const
    {
        DVector px = definePX(_xOrder,x);
        DVector py = definePY(_yOrder,y);
        double result = EIGEN_ToScalar(EIGEN_Transpose(px) * (*_coeffs) * py);
        return result;
    }

    std::auto_ptr<Function2D> Function2D::read(std::istream& fin) 
    {
        char fc,tc;

        fin >> fc >> tc;
        if (tc != 'D' && tc != 'C') fin.putback(tc);
        std::auto_ptr<Function2D> ret;
        switch(fc) {
          case 'C' : ret.reset(new Constant2D(fin));
                     break;
          case 'P' : ret.reset(new Polynomial2D(fin));
                     break;
#ifdef LEGENDRE2D_H
          case 'L' : ret.reset(new Legendre2D(fin));
                     break;
#endif
#ifdef CHEBY2D_H
          case 'X' : ret.reset(new Cheby2D(fin));
                     break;
#endif
          default: throw std::runtime_error("invalid type");
        }
        return ret;
    }

    void Function2D::linearTransform(
        double a, double b, double c,
        const Function2D& f, const Function2D& g)
    {
        if(dynamic_cast<Constant2D*>(this)) {
            Assert(dynamic_cast<const Constant2D*>(&f));
            Assert(dynamic_cast<const Constant2D*>(&g));
        }
        if(dynamic_cast<Polynomial2D*>(this)) {
            Assert(dynamic_cast<const Polynomial2D*>(&f));
            Assert(dynamic_cast<const Polynomial2D*>(&g));
        }
#ifdef LEGENDRE2D_H
        if(dynamic_cast<Legendre2D*>(this)) {
            Assert(dynamic_cast<const Legendre2D*>(&f));
            Assert(dynamic_cast<const Legendre2D*>(&g));
            Assert(dynamic_cast<const Legendre2D*>(&f)->getBounds() ==
                   dynamic_cast<const Legendre2D*>(&g)->getBounds());
        }
#endif
#ifdef CHEBY2D_H
        if(dynamic_cast<Cheby2D*>(this)) {
            Assert(dynamic_cast<const Cheby2D*>(&f));
            Assert(dynamic_cast<const Cheby2D*>(&g));
            Assert(dynamic_cast<const Cheby2D*>(&f)->getBounds() ==
                   dynamic_cast<const Cheby2D*>(&g)->getBounds());
        }
#endif
        Assert(f.getXOrder() == g.getXOrder());
        Assert(f.getYOrder() == g.getYOrder());

        if (_xOrder != f.getXOrder() || _yOrder != f.getYOrder()) {
            _xOrder = f.getXOrder();
            _yOrder = f.getYOrder();
            _coeffs.reset(new DMatrix(_xOrder+1,_yOrder+1));
            _coeffs->setZero();
        } else _coeffs->setZero();
        for(int i=0;i<=_xOrder;++i) for(int j=0;j<=_yOrder;++j) {
            (*_coeffs)(i,j) = a + b*f.getCoeffs()(i,j) + c*g.getCoeffs()(i,j);
        }
    }

    inline int fitSize(const int xOrder, const int yOrder)
    {
        int lowOrder = std::min(xOrder,yOrder);
        int highOrder = std::max(xOrder,yOrder);
        return (lowOrder+1)*(lowOrder+2)/2 + (lowOrder+1)*(highOrder-lowOrder);
    }

    void Function2D::doSimpleFit(
        int xOrder, int yOrder, 
        const std::vector<Position>& pos, const std::vector<double>& vals, 
        const std::vector<bool>& shouldUse, DVector *f,
        const std::vector<double>* sigList,
        int *dof, DVector *diff, DMatrix* cov)
    {
        // f(x,y) = Sum_pq k_pq px_p py_q
        //        = P . K (where each is a vector in pq)
        // chisq = Sum_n ((v_n - f(x,y))/s_n)^2
        // minchisq => 0 = Sum_n (v_n - f(x,y)) P_pq/s_n^2
        // => [Sum_n P_pq P/s_n^2].K = [Sum_n v_n P_pq/s_n^2]
        // Using ^ to represent an outer product, this can be written:
        //
        //    [Sum_n (P/s_n)^(P/s_n)] . K = [Sum_n (v_n/s_n) (P/s_n)]
        //
        // Or if P' is a matrix with n index for rows, and pq index for columns,
        // with each element P'(n,pq) = P_n,pq/s_n
        // and v' is a vector with v'(n) = v_n/s_n
        //
        // Then, this can be written:
        //
        // P' K = v'
        //
        // The solution to this equation, then, gives our answer for K.
        //
        Assert(pos.size() == vals.size());
        Assert(shouldUse.size() == vals.size());
        Assert((sigList==0) || (sigList->size() == vals.size()));
        xdbg<<"Start SimpleFit: size = "<<pos.size()<<std::endl;
        xdbg<<"order = "<<xOrder<<','<<yOrder<<std::endl;

        int highOrder = std::max(xOrder,yOrder);
        int size = fitSize(xOrder,yOrder);
        xdbg<<"size = "<<size<<std::endl;

        const int nVals = vals.size();

        Assert(int(f->size()) == size);
        Assert(!diff || int(diff->size()) == nVals);
        Assert(!cov || 
               (int(cov->TMV_colsize()) == size && int(cov->TMV_rowsize()) == size));

        int nUse = std::count(shouldUse.begin(),shouldUse.end(),true);
        xdbg<<"nuse = "<<nUse<<std::endl;

        DMatrix P(nUse,size);
        P.setZero();
        DVector V(nUse);

        int ii=0;
        for(int i=0;i<nVals;++i) if (shouldUse[i]) {
            if (sigList) {
                Assert((*sigList)[i] > 0.);
                V(ii) = vals[i]/(*sigList)[i];
            } else {
                V(ii) = vals[i];
            }

            DVector px = definePX(xOrder,pos[i].getX());
            DVector py = definePY(yOrder,pos[i].getY());
            int pq=0;
            for(int pplusq=0;pplusq<=highOrder;++pplusq) { 
                for(int p=std::min(pplusq,xOrder),q=pplusq-p;
                    q<=std::min(pplusq,yOrder);--p,++q,++pq) {
                    Assert(p<int(px.size()));
                    Assert(q<int(py.size()));
                    Assert(pq<int(P.TMV_rowsize()));
                    P(ii,pq) = px[p]*py[q];
                }
            }
            Assert(pq == int(P.TMV_rowsize()));
            if (sigList) P.row(ii) /= (*sigList)[i];
            ++ii;
        }
        Assert(ii==nUse);

        xdbg<<"Done make V,P\n";
        xdbg<<"V = "<<EIGEN_Transpose(V)<<std::endl;
        //P.divideUsing(tmv::QR);
        //P.saveDiv();
        //*f = V/P;
        TMV_QR(P);
        TMV_QR_Solve(P,*f,V);
        xdbg<<"*f = "<<EIGEN_Transpose(*f)<<std::endl;

        if (diff) {
            int k=0;
            for(int i=0;i<nVals;++i) {
                if (shouldUse[i]) {
                    (*diff)(i) = 
                        V(k) - EIGEN_ToScalar(P.row(k) * (*f));
                    ++k;
                } else {
                    (*diff)(i) = 0.;
                }
            }
        }
        if (dof) {
            *dof = P.TMV_colsize() - P.TMV_rowsize();
            if (*dof < 0) *dof = 0;
        }
        if (cov) {
            //P.makeInverseATA(*cov);
            TMV_QR_InverseATA(P,*cov);
        }
        xdbg<<"Done simple fit\n";
    }

    void Function2D::simpleFit(
        int order, const std::vector<Position>& pos, 
        const std::vector<double>& vals, const std::vector<bool>& shouldUse, 
        const std::vector<double>* sig,
        double *chisqOut, int *dofOut, DMatrix* cov) 
    {
        DVector fVect(fitSize(order,order));
        if (chisqOut) {
            DVector diff(vals.size());
            doSimpleFit(order,order,pos,vals,shouldUse,&fVect,sig,dofOut,&diff,cov);
            *chisqOut = diff.TMV_normSq();
        } else {
            doSimpleFit(order,order,pos,vals,shouldUse,&fVect,sig);
        }
        setFunction(order,order,fVect);
    }

    inline double absSq(const double& x) { return x*x; }

    //inline double absSq(const std::complex<double>& x) { return std::norm(x); }

    void Function2D::outlierFit(
        int order,double nSig,
        const std::vector<Position>& pos, const std::vector<double>& vals,
        std::vector<bool> *shouldUse,
        const std::vector<double>* sig, double *chisqOut, int *dofOut, 
        DMatrix *cov)
    {
        xdbg<<"start outlier fit\n";
        const int nVals = vals.size();
        bool isDone=false;
        DVector fVect(fitSize(order,order));
        int dof;
        double chisq=0.;
        double nSigSq = nSig*nSig;
        while (!isDone) {
            DVector diff(nVals);
            xdbg<<"before dosimple\n";
            doSimpleFit(order,order,pos,vals,*shouldUse,&fVect,sig,&dof,&diff,cov);
            xdbg<<"after dosimple\n";
            // Caclulate chisq, keeping the vector diffsq for later when
            // looking for outliers
            chisq = diff.TMV_normSq();
            // If sigmas are given, then chisq should be 1, since diff is
            // already normalized by sig_i.  But if not, then this effectively
            // assumes that all the given errors are off by some uniform factor.

            // Look for outliers, setting isDone = false if any are found
            isDone = true;
            if (dof <= 0) break;
            double thresh=nSigSq*chisq/dof;
            xdbg<<"chisq = "<<chisq<<", thresh = "<<thresh<<std::endl;
            for(int i=0;i<nVals;++i) if( (*shouldUse)[i]) { 
                xdbg<<"pos ="<<pos[i]<<", v = "<<vals[i]<<
                    " diff = "<<diff[i]<<std::endl;
                if (absSq(diff[i]) > thresh) {
                    isDone = false; 
                    (*shouldUse)[i] = false;
                    xdbg<<i<<" ";
                }
            }
            if (!isDone) xdbg<<" are outliers\n";
        }
        setFunction(order,order,fVect);
        if (chisqOut) *chisqOut = chisq;
        if (dofOut) *dofOut = dof;
        xdbg<<"done outlier fit\n";
    }

    inline double betai(double a,double b,double x);

    inline bool Equivalent(double chisq1,double chisq2, int n1, int n2,
                           double equivProb)
    {
        if (chisq2 <= chisq1) return true;
        // should only happpen when essentially equal but have rounding errors
        if (chisq1 <= 0.) return (chisq2 <= 0.);

        Assert(n1 < n2);
        if (n1 <= 0) return (n2 <= 0);
        double f = (chisq2-chisq1)/(n2-n1) / (chisq1/n1);

        double prob = betai(0.5*n2,0.5*n1,n2/(n2+n1*f));
        // = probability that these chisq would happen by chance if equiv
        // (technically if underlying chisq2 really smaller or = to chisq1)
        // so equiv if prob is large, not equiv if prob is small.
        return (prob > 1.-equivProb);
    }

    void Function2D::orderFit(
        int maxOrder, double equivProb,
        const std::vector<Position>& pos,const std::vector<double>& vals,
        const std::vector<bool>& shouldUse, const std::vector<double>* sig,
        double *chisqOut, int *dofOut, DMatrix* cov) 
    {
        xdbg<<"Start OrderFit\n";
        DVector fVectMax(fitSize(maxOrder,maxOrder));
        DVector diff(vals.size());
        int dofmax;
        doSimpleFit(maxOrder,maxOrder,pos,vals,shouldUse,
                    &fVectMax,sig,&dofmax,&diff,cov);
        double chisqmax = diff.TMV_normSq();
        xdbg<<"chisq,dof(n="<<maxOrder<<") = "<<chisqmax<<','<<dofmax<<std::endl;
        int tryOrder;
        double chisq=-1.;
        int dof=-1;
        std::auto_ptr<DVector > fVect(0);
        for(tryOrder=0;tryOrder<maxOrder;++tryOrder) {
            fVect.reset(new DVector(fitSize(tryOrder,tryOrder)));
            doSimpleFit(tryOrder,tryOrder,pos,vals,shouldUse,
                        fVect.get(),sig,&dof,&diff,cov);
            chisq = diff.TMV_normSq();
            xdbg<<"chisq,dof(n="<<tryOrder<<") = "<<chisq<<','<<dof<<"....  ";
            if (Equivalent(chisqmax,chisq,dofmax,dof,equivProb)) {
                xdbg<<"equiv\n";
                break;
            }
            xdbg<<"not equiv\n";
        }
        if (tryOrder == maxOrder) {
            setFunction(tryOrder,tryOrder,fVectMax);
            if(chisqOut) *chisqOut = chisqmax;
            if(dofOut) *dofOut = dofmax;
        } else {
            Assert(fVect.get());
            setFunction(tryOrder,tryOrder,*fVect);
            if(chisqOut) *chisqOut = chisq;
            if(dofOut) *dofOut = dof;
        }
    }

    // These are numerical recipes routines:

#include <cmath>
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

    inline double betacf(double a, double b, double x)
    {
        double qab = a+b;
        double qap = a+1.0;
        double qam = a-1.0;
        double c = 1.0;
        double d = 1.0 - qab*x/qap;
        if (std::abs(d) < FPMIN) d = FPMIN;
        d = 1.0/d;
        double h = d;
        int m;
        for(m=1;m<=MAXIT;++m) {
            int m2 = 2*m;
            double aa = m*(b-m)*x/((qam+m2)*(a+m2));
            d = 1.0+aa*d;
            if (std::abs(d) < FPMIN) d = FPMIN;
            c = 1.0+aa/c;
            if (std::abs(c) < FPMIN) c = FPMIN;
            d = 1.0/d;
            h *= d*c;
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
            d = 1.0+aa*d;
            if (std::abs(d) < FPMIN) d = FPMIN;
            c = 1.0+aa/c;
            if (std::abs(c) < FPMIN) c = FPMIN;
            d = 1.0/d;
            double del = d*c;
            h *= del;
            if (std::abs(del-1.0) < EPS) break;
        }
        if (m>MAXIT) {
            throw std::runtime_error(
                "a or b too big in betacf, or MAXIT too small");
        }
        return h;
    }

#undef MAXIT
#undef EPS
#undef FPMIN

    inline double gammln(double x)
    {
        const double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
            -0.5395239384953e-5};
        double temp = x+5.5;
        temp -= (x+0.5)*log(temp);
        double ser = 1.000000000190015;
        double y=x;
        for(int j=0;j<6;++j) ser += cof[j]/(y+=1.0);
        return -temp+log(2.5066282746310005*ser/x);
    }

    inline double betai(double a,double b,double x)
    {
        if (x<0.0 || x>1.0) throw std::runtime_error("Bad x in betai");
        if (x==0.0) return 0.;
        if (x==1.0) return 1.;
        double bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
        if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
        else return 1.0-bt*betacf(b,a,1.0-x)/b;
    }

}}}}
