
#include "lsst/meas/algorithms/shapelet/EllipseSolver.h"
#include "lsst/meas/algorithms/shapelet/PsiHelper.h"

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    static int calculateMaxPsfOrder(const std::vector<BVec>& psf)
    {
        int maxpsforder = 0;
        for(size_t k=0;k<psf.size();++k) 
            if (psf[k].getOrder() > maxpsforder)
                maxpsforder = psf[k].getOrder();
        return maxpsforder;
    }

    static size_t calculateSumSize(const std::vector<PixelList>& v)
    {
        size_t sum = 0;
        for(size_t i=0;i<v.size();++i) sum += v[i].size();
        return sum;
    }

    //
    // EllipseSolver2
    //

    struct EllipseSolver2::ESImpl2 
    {

    public :

        ESImpl2(const std::vector<PixelList>& _pix,
                int _order, double _sigma, double _pixScale,
                bool _isFixedCen, bool _isFixedGamma, bool _isFixedMu,
                bool _shouldUseFlux);

        ESImpl2(const std::vector<PixelList>& _pix,
                const std::vector<BVec>& _psf, double _fp,
                int _order, double _sigma, double _pixScale,
                bool _isFixedCen, bool _isFixedGamma, bool _isFixedMu,
                bool _shouldUseFlux);

        void calculateF(const DVector& x, DVector& f) const;

        void calculateJ(
            const DVector& x, const DVector& f,
            DMatrix& J) const;

        int nsize, np2size;
        mutable BVec b;
        const std::vector<PixelList>& pix;
        const std::vector<BVec>* psf;
        double f_psf;
        std::vector<double> sigma_obs;
        DVector I;
        CDVector Z;
        mutable CDVector Z1;
        mutable std::vector<DMatrix > A_aux;
        mutable std::vector<DMatrixView > A;
        mutable std::vector<DMatrix > C;
        mutable std::vector<DMatrix > S;
        mutable std::vector<DMatrix > D;
        bool isFixedCen, isFixedGamma, isFixedMu, shouldUseFlux, numeric_j;
        DMatrix U;
        mutable DVector xinit;
        mutable DVector xx;
        mutable DVector x_short;
        mutable DVector f_short;
        mutable DVector b1;
        mutable DMatrix dbdE;
        mutable DMatrix dbdE1;
        mutable DVector AtI;
        mutable DMatrix dCdE;
        mutable double fixuc,fixvc,fixg1,fixg2,fixm,flux;
        double pixScale;
        DDiagMatrix normD;
        int maxpsforder;
        DMatrix _Gx;
        DMatrix _Gy;
        DMatrix _Gg1;
        DMatrix _Gg2;
        DMatrix _Gth;
        DMatrix _Gmu;
        mutable DMatrix _dTdE;
    };

    EllipseSolver2::EllipseSolver2(
        const std::vector<PixelList>& pix,
        int order, double sigma, double pixScale,
        bool isFixedCen, bool isFixedGamma, bool isFixedMu, bool shouldUseFlux) :
        _pimpl(
            new EllipseSolver2::ESImpl2(
                pix,order,sigma,pixScale,
                isFixedCen,isFixedGamma,isFixedMu,shouldUseFlux))
    {}

    EllipseSolver2::EllipseSolver2(
        const std::vector<PixelList>& pix,
        const std::vector<BVec>& psf, double fp,
        int order, double sigma, double pixScale,
        bool isFixedCen, bool isFixedGamma, bool isFixedMu, bool shouldUseFlux) :
        _pimpl(
            new EllipseSolver2::ESImpl2(
                pix,psf,fp,order,sigma,pixScale,
                isFixedCen,isFixedGamma,isFixedMu,shouldUseFlux))
    {}

    EllipseSolver2::~EllipseSolver2()
    { delete _pimpl; }

    void EllipseSolver2::calculateF(const DVector& x, DVector& f) const
    { _pimpl->calculateF(x,f); }

    void EllipseSolver2::calculateJ(
        const DVector& x, const DVector& f, DMatrix& J) const
    {
        if (_pimpl->numeric_j) NLSolver::calculateJ(x,f,J);
        else _pimpl->calculateJ(x,f,J); 
    }

    EllipseSolver2::ESImpl2::ESImpl2(
        const std::vector<PixelList>& _pix,
        int _order, double _sigma, double _pixScale,
        bool _isFixedCen, bool _isFixedGamma, bool _isFixedMu,
        bool _shouldUseFlux) :

        nsize((_order+1)*(_order+2)/2),
        np2size((_order+3)*(_order+4)/2), b(_order,_sigma),
        pix(_pix), psf(0), f_psf(1.), sigma_obs(pix.size(),_sigma), 
        I(calculateSumSize(pix)), Z(I.size()), Z1(I.size()),
        isFixedCen(_isFixedCen), isFixedGamma(_isFixedGamma),
        isFixedMu(_isFixedMu), shouldUseFlux(_shouldUseFlux),
        numeric_j(false),
        U((isFixedCen?0:2)+(isFixedGamma?0:2)+(isFixedMu?0:1),5),
        xinit(5), xx(5), x_short(U.TMV_colsize()),
        f_short(U.TMV_colsize()+(shouldUseFlux?1:0)),
        b1(nsize), dbdE(nsize,5), dbdE1(nsize,5), AtI(np2size),
        dCdE(1,1), pixScale(_pixScale), normD(b.size()), maxpsforder(0),
        _Gx(np2size,nsize), _Gy(np2size,nsize),
        _Gg1(np2size,nsize), _Gg2(np2size,nsize),
        _Gth(np2size,nsize), _Gmu(np2size,nsize),
        _dTdE(1,1)
        {
            //dbg<<"Start ESImpl2"<<std::endl;
            U.setZero();
            xinit.setZero();

            A_aux.reserve(pix.size());
            A.reserve(pix.size());
            for(size_t k=0,n=0;k<pix.size();++k) {
                for(size_t i=0;i<pix[k].size();++i,++n) {
                    I(n) = pix[k][i].getFlux();
                    Z(n) = pix[k][i].getPos();
                }
                A_aux.push_back(DMatrix(pix[k].size(),np2size));
                A.push_back(TMV_colRange(A_aux[k],0,nsize));
            }

            int j=0;
            if (!isFixedCen) { U(j++,0) = 1.; U(j++,1) = 1.; }
            if (!isFixedGamma) { U(j++,2) = 1.; U(j++,3) = 1.; }
            if (!isFixedMu) { U(j,4) = 1.; }

            setupGx(_Gx,_order+2,_order);
            setupGy(_Gy,_order+2,_order);
            setupGg1(_Gg1,_order+2,_order);
            setupGg2(_Gg2,_order+2,_order);
            setupGth(_Gth,_order+2,_order);
            setupGmu(_Gmu,_order+2,_order);

            double val = pixScale*pixScale/(_sigma*_sigma)/pix.size();
            for(int n=0,k=0;n<=b.getOrder();++n) {
                for(int p=n,q=0;p>=q;--p,++q,++k) {
                    if (p!=q) { normD(k) = val/2.; normD(++k) = val/2.; }
                    else normD(k) = val;
                }
            }
            //dbg<<"Done ESImpl2"<<std::endl;
        }

#define ORDER_1 (std::max(maxpsforder,_order))
#define ORDER_2 (std::max(maxpsforder,_order+2))
#define ORDER_3 (std::max(maxpsforder+2,_order))
#define SIZE_1 (ORDER_1+1)*(ORDER_1+2)/2
#define SIZE_1b (ORDER_1+3)*(ORDER_1+4)/2
#define SIZE_2 (ORDER_2+1)*(ORDER_2+2)/2
#define SIZE_3 (ORDER_3+1)*(ORDER_3+2)/2
#define SIZE_4 (maxpsforder+1)*(maxpsforder+2)/2

    EllipseSolver2::ESImpl2::ESImpl2(
        const std::vector<PixelList>& _pix,
        const std::vector<BVec>& _psf, double _fp, int _order, double _sigma, 
        double _pixScale, bool _isFixedCen, bool _isFixedGamma,
        bool _isFixedMu, bool _shouldUseFlux) :

        nsize((_order+1)*(_order+2)/2),
        np2size((_order+3)*(_order+4)/2), b(_order,_sigma),
        pix(_pix), psf(&_psf), f_psf(_fp), sigma_obs(pix.size()), 
        I(calculateSumSize(pix)), Z(I.size()), Z1(I.size()),
        isFixedCen(_isFixedCen), isFixedGamma(_isFixedGamma),
        isFixedMu(_isFixedMu), shouldUseFlux(_shouldUseFlux),
        numeric_j(false),
        U((isFixedCen?0:2)+(isFixedGamma?0:2)+(isFixedMu?0:1),5),
        xinit(5), xx(5), x_short(U.TMV_colsize()),
        f_short(U.TMV_colsize()+(shouldUseFlux?1:0)),
        b1(nsize), dbdE(nsize,5), dbdE1(nsize,5), AtI(np2size),
        dCdE(nsize,nsize), pixScale(_pixScale), normD(b.size()),
        maxpsforder(calculateMaxPsfOrder(_psf)),
        _Gx(np2size,nsize), _Gy(np2size,nsize),
        _Gg1(SIZE_1b,SIZE_1), _Gg2(SIZE_1b,SIZE_1),
        _Gth(SIZE_1b,SIZE_1), _Gmu(SIZE_2,SIZE_3),
        _dTdE(SIZE_4,SIZE_4)

        {
            U.setZero();
            xinit.setZero();

            int maxpsforder = 0;
            A_aux.reserve(pix.size());
            A.reserve(pix.size());
            C.reserve(pix.size());
            S.reserve(pix.size());
            D.reserve(pix.size());
            for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
                nx = n+pix[k].size();

                for(size_t i=0,n1=n;i<pix[k].size();++i,++n1) {
                    I(n1) = pix[k][i].getFlux();
                    Z(n1) = pix[k][i].getPos();
                }

                int psfsize = (*psf)[k].size();
                int psfsize2 = ((*psf)[k].getOrder()+3)*((*psf)[k].getOrder()+4)/2;
                A_aux.push_back(
                    DMatrix(pix[k].size(),np2size));
                A.push_back(TMV_colRange(A_aux[k],0,nsize));
                C.push_back(DMatrix(nsize,nsize));
                S.push_back(DMatrix(psfsize,psfsize2));
                D.push_back(DMatrix(psfsize2,psfsize));
                if ((*psf)[k].getOrder() > maxpsforder)
                    maxpsforder = (*psf)[k].getOrder();
                sigma_obs[k] = sqrt(pow(_sigma,2)+f_psf*pow((*psf)[k].getSigma(),2));
#ifdef USE_TMV
                C[k].saveDiv();
#endif

                I.TMV_subVector(n,nx) *= pow(_sigma,2) / pow(sigma_obs[k],2);
            }

            int j=0;
            if (!isFixedCen) { U(j++,0) = 1.; U(j++,1) = 1.; }
            if (!isFixedGamma) { U(j++,2) = 1.; U(j++,3) = 1.; }
            if (!isFixedMu) { U(j,4) = 1.; }

            setupGx(_Gx,_order+2,_order);
            setupGy(_Gy,_order+2,_order);
            setupGg1(_Gg1,ORDER_1+2,ORDER_1);
            setupGg2(_Gg2,ORDER_1+2,ORDER_1);
            setupGth(_Gth,ORDER_1+2,ORDER_1);
            setupGmu(_Gmu,ORDER_2,ORDER_3);

            double val = pixScale*pixScale/(_sigma*_sigma)/pix.size();
            for(int n=0,k=0;n<=b.getOrder();++n) {
                for(int p=n,q=0;p>=q;--p,++q,++k) {
                    if (p!=q) { normD(k) = val/2.; normD(++k) = val/2.; }
                    else normD(k) = val;
                }
            }
        }
#undef ORDER_1
#undef ORDER_2
#undef ORDER_3
#undef SIZE_1
#undef SIZE_1b
#undef SIZE_2
#undef SIZE_3

    void EllipseSolver2::ESImpl2::calculateF(const DVector& x, DVector& f) const
    {
        Assert(x.size() == U.TMV_colsize());
        if (shouldUseFlux) Assert(f.size() == U.TMV_colsize()+1);
        else Assert(f.size() == U.TMV_colsize());
        // x(0),x(1) = u_cen, v_cen
        // x(2),x(3) = gamma_1, gamma_2
        // x(4) = mu
        //
        // ( u )' = exp(-mu)/sqrt(1-gsq) ( 1-g1  -g2  ) ( u-uc )
        // ( v )                         (  -g2  1+g1 ) ( v-vc )
        // 

        EIGEN_Transpose(xx) = EIGEN_Transpose(x)*U;
        if (isFixedCen) { xx[0] = fixuc;  xx[1] = fixvc; }
        if (isFixedGamma) { xx[2] = fixg1;  xx[3] = fixg2; }
        if (isFixedMu) { xx[4] = fixm; }

        if ((xx-xinit).TMV_normInf() > 4.) {
            if (flux == 0.) flux = b(0);
            if (shouldUseFlux) {
                f(0) = 2.*(b(0)/flux-1.);
                f.TMV_subVector(1,f.size()) = 2.*U*b.vec().TMV_subVector(1,6)/flux; 
            }
            else f = 2.*U*b.vec().TMV_subVector(1,6)/flux;
            return; 
        }

        std::complex<double> zc(xx[0],xx[1]);
        std::complex<double> g(xx[2],xx[3]);
        double mu = xx[4];
        double gsq = std::norm(g);

        if (gsq > 0.99 || (mu < -2. && (xx-xinit).norm() > 0.3)) {
            if (flux == 0.) flux = b(0);
            if (flux == 0.) flux = 1.;
            if (shouldUseFlux) {
                f(0) = 2.*(b(0)/flux-1.);
                f.TMV_subVector(1,f.size()) = 2.*U*b.vec().TMV_subVector(1,6)/flux; 
            } else {
                f = 2.*U*b.vec().TMV_subVector(1,6)/flux; 
            }
            return; 
        }


        double m0 = exp(-mu)/sqrt(1.-std::norm(g));
        b.vec().setZero();

        // z' = m*(z-zc) - m*g*conj(z-zc)
        //    = m*z - m*g*conj(z) - m*zc + m*g*conj(zc);
        Z1 = m0*Z;
        Z1 -= m0*g*Z.conjugate();
        Z1.TMV_addToAll(m0*(g*std::conj(zc)-zc));
        for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
            nx = n+pix[k].size();
            //dbg<<"k,n,nx = "<<k<<','<<n<<','<<nx<<std::endl;

            Z1.TMV_subVector(n,nx) /= sigma_obs[k];

            makePsi(A_aux[k],Z1.TMV_subVector(n,nx),b.getOrder());
            // b = C^-1 normD AT I
            b1 = A[k].transpose() * I.TMV_subVector(n,nx);
            b1 = normD EIGEN_asDiag() * b1;

            if (psf) {
                int psfsize = (*psf)[k].size();

                BVec newpsf((*psf)[k].getOrder(),(*psf)[k].getSigma());
                DMatrixView S1 = TMV_colRange(S[k],0,psfsize);
                calculateGTransform(g,(*psf)[k].getOrder(),S[k]);
                newpsf.vec() = S1 * (*psf)[k].vec();
                DMatrixView D1 = TMV_rowRange(D[k],0,psfsize);
                calculateMuTransform(mu,(*psf)[k].getOrder(),D[k]);
                newpsf.vec() = D1 * newpsf.vec();
                calculatePsfConvolve(newpsf,b.getOrder(),b.getSigma(),C[k]);
#ifdef USE_TMV
                C[k].resetDiv();
                b1 /= C[k];
#else
                C[k].lu().solve(b1,&b1);
#endif
            }
            b.vec() += b1;
        }

        // First time through, set flux
        if (flux == 0.) {
            flux = b(0);
        }
        if (shouldUseFlux) {
            f(0) = b(0)/flux-1.;
            f.TMV_subVector(1,U.TMV_colsize()+1) = U*b.vec().TMV_subVector(1,6)/flux;
        } else {
            f = U*b.vec().TMV_subVector(1,6)/flux;
        }
        // Leave this out of f, but get it right in case getB is called.
        b.vec() *= exp(-2.*mu);
    }

    void EllipseSolver2::ESImpl2::calculateJ(
        const DVector& x, const DVector& f, DMatrix& J) const
    {
        // b = C^-1 normD AT I
        // db/dE = C^-1 normD (dA/dE)T I
        //         - C^-1 (dC/dE) C^-1 normD AT I
        Assert(x.size() == U.TMV_colsize());
        Assert(J.TMV_rowsize() == U.TMV_colsize());
        if (shouldUseFlux) {
            Assert(f.size() == U.TMV_colsize()+1);
            Assert(J.TMV_colsize() == U.TMV_colsize()+1);
        } else {
            Assert(f.size() == U.TMV_colsize());
            Assert(J.TMV_colsize() == U.TMV_colsize());
        }

        EIGEN_Transpose(xx) = EIGEN_Transpose(x)*U;
        if (isFixedCen) { xx[0] = fixuc;  xx[1] = fixvc; }
        if (isFixedGamma) { xx[2] = fixg1;  xx[3] = fixg2; }
        if (isFixedMu) { xx[4] = fixm; }
        std::complex<double> zc(xx[0],xx[1]);
        std::complex<double> g(xx[2],xx[3]);
        double mu = xx[4];
        double fact = 1./(1.-std::norm(g));

        double m0 = exp(-mu)*sqrt(fact);
        double g1 = std::real(g);
        double g2 = std::imag(g);

        DConstMatrixView Gx = _Gx.TMV_subMatrix(0,np2size,0,nsize);
        DConstMatrixView Gy = _Gy.TMV_subMatrix(0,np2size,0,nsize);
        DConstMatrixView Gg1 = _Gg1.TMV_subMatrix(0,np2size,0,nsize);
        DConstMatrixView Gg2 = _Gg2.TMV_subMatrix(0,np2size,0,nsize);
        DConstMatrixView Gmu = _Gmu.TMV_subMatrix(0,np2size,0,nsize);
        DConstMatrixView Gth = _Gth.TMV_subMatrix(0,np2size,0,nsize);

        dbdE.setZero();
        for(size_t k=0,n=0,nx;k<pix.size();++k,n=nx) {
            nx = n+pix[k].size();
            double m = m0/sigma_obs[k];

            augmentPsi(A_aux[k],Z1.TMV_subVector(n,nx),b.getOrder());

            AtI = A_aux[k].transpose() * I.TMV_subVector(n,nx);

            dbdE1.col(0) = -m*(1.-g1) * Gx.transpose() * AtI;
            dbdE1.col(0) += m*g2 * Gy.transpose() * AtI;

            dbdE1.col(1) = -m*(1.+g1) * Gy.transpose() * AtI;
            dbdE1.col(1) += m*g2 * Gx.transpose() * AtI;

            dbdE1.col(2) = -fact * Gg1.transpose() * AtI;
            dbdE1.col(2) += fact*g2 * Gth.transpose() * AtI;

            dbdE1.col(3) = -fact * Gg2.transpose() * AtI;
            dbdE1.col(3) -= fact*g1 * Gth.transpose() * AtI;

            dbdE1.col(4) = -Gmu.transpose() * AtI;

            dbdE1 = normD EIGEN_asDiag() * dbdE1;

            // So far dbdE1 = normD * (dA/dE)T * I

            if (psf) {

                // dbdE -= dCdE C^-1 normD AT I
                int n2 = ((*psf)[k].getOrder()+3)*((*psf)[k].getOrder()+4)/2;
                int psfsize = (*psf)[k].size();

                b1 = normD EIGEN_asDiag() * AtI.TMV_subVector(0,nsize);
#ifdef USE_TMV
                Assert(C[k].divIsSet());
                b1 /= C[k];
#else
                Eigen::LU<DMatrix> LU_Solver_C = C[k].lu();
                LU_Solver_C.solve(b1,&b1);
#endif


                DMatrixView dTdE = _dTdE.TMV_subMatrix(0,psfsize,0,psfsize);
                DConstMatrixView D1 = TMV_rowRange(D[k],0,psfsize);
                DConstMatrixView S1 = TMV_colRange(S[k],0,psfsize);

                augmentGTransformCols(g,(*psf)[k].getOrder(),S[k]);
                augmentMuTransformRows(mu,(*psf)[k].getOrder(),D[k]);

                // E = g1
                dTdE = fact*D1*S[k]*_Gg1.TMV_subMatrix(0,n2,0,psfsize);
                dTdE += fact*g2*D1*S[k]*_Gth.TMV_subMatrix(0,n2,0,psfsize);
                BVec dbpsfdE((*psf)[k].getOrder(),(*psf)[k].getSigma());
                dbpsfdE.vec() = dTdE*(*psf)[k].vec();
                calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
                dbdE1.col(2) -= dCdE*b1;

                // E = g2
                dTdE = fact*D1*S[k]*_Gg2.TMV_subMatrix(0,n2,0,psfsize);
                dTdE -= fact*g1*D1*S[k]*_Gth.TMV_subMatrix(0,n2,0,psfsize);
                dbpsfdE.vec() = dTdE*(*psf)[k].vec();
                calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
                dbdE1.col(3) -= dCdE*b1;

                // E = mu
                dTdE = _Gmu.TMV_subMatrix(0,psfsize,0,n2)*D[k]*S1;
                dbpsfdE.vec() = dTdE*(*psf)[k].vec();
                calculatePsfConvolve(dbpsfdE,b.getOrder(),b.getSigma(),dCdE);
                dbdE1.col(4) -= dCdE*b1;

#ifdef USE_TMV
                Assert(C[k].divIsSet());
                dbdE1 /= C[k];
#else
                LU_Solver_C.solve(dbdE1,&dbdE1);
#endif
                // db/dE = C^-1 normD (dA/dE)T I
                //         - C^-1 (dC/dE) C^-1 normD AT I
            }

            dbdE += dbdE1;
        }

        if (shouldUseFlux) {
            J.row(0) = dbdE.row(0)*U.transpose();
            TMV_rowRange(J,1,U.TMV_colsize()+1) = U*TMV_rowRange(dbdE,1,6)*U.transpose();
        } else {
            J = U*TMV_rowRange(dbdE,1,6)*U.transpose();
        }
        J /= flux;
    }

    void EllipseSolver2::callF(
        const DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);
        _pimpl->xinit = x;
        if (_pimpl->isFixedCen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
        if (_pimpl->isFixedGamma) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
        if (_pimpl->isFixedMu) { _pimpl->fixm = x[4]; }

        _pimpl->x_short = _pimpl->U * x;

        _pimpl->flux = 0;
        calculateF(_pimpl->x_short,_pimpl->f_short);

        if (_pimpl->shouldUseFlux) 
            EIGEN_Transpose(f) = 
                EIGEN_Transpose(_pimpl->f_short.TMV_subVector(1,_pimpl->f_short.size())) * 
                _pimpl->U;
        else EIGEN_Transpose(f) = EIGEN_Transpose(_pimpl->f_short)*_pimpl->U;
    }

    bool EllipseSolver2::solve(DVector& x, DVector& f) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);
        _pimpl->xinit = x;
        if (_pimpl->isFixedCen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
        if (_pimpl->isFixedGamma) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
        if (_pimpl->isFixedMu) { _pimpl->fixm = x[4]; }

        _pimpl->x_short = _pimpl->U * x;

        _pimpl->flux = 0;
        bool ret = NLSolver::solve(_pimpl->x_short,_pimpl->f_short);

        EIGEN_Transpose(x) = EIGEN_Transpose(_pimpl->x_short)*_pimpl->U;
        if (_pimpl->isFixedCen) { x[0] = _pimpl->fixuc; x[1] = _pimpl->fixvc; }
        if (_pimpl->isFixedGamma) { x[2] = _pimpl->fixg1; x[3] = _pimpl->fixg2; }
        if (_pimpl->isFixedMu) { x[4] = _pimpl->fixm; }
        if (_pimpl->shouldUseFlux) 
            EIGEN_Transpose(f) = 
                EIGEN_Transpose(_pimpl->f_short.TMV_subVector(1,_pimpl->f_short.size())) * 
                _pimpl->U;
        else EIGEN_Transpose(f) = EIGEN_Transpose(_pimpl->f_short)*_pimpl->U;

        return ret;
    }

    bool EllipseSolver2::testJ(
        const DVector& x, DVector& f,
        std::ostream* os, double relErr) const
    {
        Assert(x.size() == 5);
        Assert(f.size() == 5);
        _pimpl->xinit = x;
        if (_pimpl->isFixedCen) { _pimpl->fixuc = x[0]; _pimpl->fixvc = x[1]; }
        if (_pimpl->isFixedGamma) { _pimpl->fixg1 = x[2]; _pimpl->fixg2 = x[3]; }
        if (_pimpl->isFixedMu) { _pimpl->fixm = x[4]; }

        _pimpl->x_short = _pimpl->U * x;

        _pimpl->flux = 0;
        return NLSolver::testJ(_pimpl->x_short,_pimpl->f_short,os,relErr);
    }

    void EllipseSolver2::useNumericJ() { _pimpl->numeric_j = true; }

    const BVec& EllipseSolver2::getB() const { return _pimpl->b; }

    void EllipseSolver2::getBCov(DMatrix&) const 
    { Assert(false); }

}}}}
