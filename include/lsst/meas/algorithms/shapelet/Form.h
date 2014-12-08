#ifndef MeasAlgoShapeletForm_H
#define MeasAlgoShapeletForm_H

// Formatting for simplified stream output
// see Stroustrup(2000), p. 635

#include <complex>
#include <sstream>
#include <iostream>
#include <string>

namespace lsst {
namespace meas {
namespace algorithms {
namespace shapelet {

    template <typename T> class BoundForm;

    class Form 
    {
    public:
        Form() : 
            _prc(6), _wdt(0), _fmt(), _base(std::ios_base::dec),
            _just(std::ios_base::left), _newFillCh(0),
            _doUpper(0), _doPlus(0), _doTrail(0), _doBoolAlpha(0),
            _nTrail(1), _trailCh(' ')  
        {}

        template <typename T> BoundForm<T> operator()(T val) const;

        Form& prec(int p) { _prc = p; return *this; }

        Form& sci() { _fmt = std::ios_base::scientific; return *this; }
        Form& fix() { _fmt = std::ios_base::fixed; return *this; }
        Form& gen() { _fmt = ~std::ios_base::floatfield; return *this; }

        Form& width(int w) { _wdt = w; return *this; }
        Form& fill(char c) { _newFillCh = c; return *this; }

        Form& dec() { _base = std::ios_base::dec; return *this; }
        Form& oct() { _base = std::ios_base::oct; return *this; }
        Form& hex() { _base = std::ios_base::hex; return *this; }

        Form& left() { _just = std::ios_base::left; return *this; }
        Form& right() { _just = std::ios_base::right; return *this; }
        Form& internal() { _just = std::ios_base::internal; return *this; }

        Form& uppercase(bool b=true) { _doUpper = b?1:-1; return *this; }
        Form& showpos(bool b=true) { _doPlus = b?1:-1; return *this; }
        Form& showpoint(bool b=true) { _doTrail = b?1:-1; return *this; }
        Form& boolalpha(bool b=true) { _doBoolAlpha = b?1:-1; return *this; }

        Form& trail(int n, char ch=' ') 
        { _nTrail = n; _trailCh = ch; return *this; }

    private:
        template <typename T> 
        friend std::ostream& operator<<(std::ostream&, const BoundForm<T>&);

        friend void setupFloat(std::ostream&, const Form&);
        friend void setupInt(std::ostream&, const Form&);

        int _prc; // precision
        int _wdt; // width, 0 means as wide as necessary
        std::ios_base::fmtflags _fmt; // general sci, or fixed
        std::ios_base::fmtflags _base; // dec, hex, oct
        std::ios_base::fmtflags _just; // left, right, internal fill
        char _newFillCh; // fill character
        int _doUpper; // +1 to have uppercase E,X, -1 turn off, (0 leave as is)
        int _doPlus; // +1 to have explicit plus for positive values, -1 off, 0 same
        int _doTrail; // +1 to write trailing zeros, -1 off, 0 same
        int _doBoolAlpha; // +1 to write "true","false", -1 off, 0 same
        int _nTrail; // number of spaces after output
        char _trailCh; // character of trailing "spaces"

    };

    template <typename T>
    class BoundForm 
    {
    public:
        const Form& f;
        T val;
        BoundForm(const Form& f_, T val_) : f(f_), val(val_) {}
    };

    template <typename T>
    inline BoundForm<T> Form::operator()(T val) const 
    { return BoundForm<T>(*this,val); }

    inline void setupFloat(std::ostream& s, const Form& f)
    {
        s.precision(f._prc);
        s.setf(f._fmt,std::ios_base::floatfield);
        s.setf(f._just,std::ios_base::adjustfield);
        if (f._wdt) s.width(f._wdt);
        if (f._newFillCh) s.fill(f._newFillCh);
        if (f._doUpper && f._fmt == std::ios_base::scientific) {
            if (f._doUpper>0) s.setf(std::ios_base::uppercase);
            else s.unsetf(std::ios_base::uppercase); 
        }
        if (f._doPlus) {
            if (f._doPlus>0) s.setf(std::ios_base::showpos); 
            else s.unsetf(std::ios_base::showpos); 
        }
        if (f._doTrail) {
            if (f._doTrail>0) s.setf(std::ios_base::showpoint); 
            else s.unsetf(std::ios_base::showpoint); 
        }
    }

    inline void setupInt(std::ostream& s, const Form& f)
    {
        s.setf(f._just,std::ios_base::adjustfield);
        s.setf(f._base,std::ios_base::basefield);
        if (f._wdt) s.width(f._wdt);
        if (f._newFillCh) s.fill(f._newFillCh);
        if (f._doUpper && f._base == std::ios_base::hex) {
            if (f._doUpper>0) s.setf(std::ios_base::uppercase); 
            else s.unsetf(std::ios_base::uppercase); 
        }
        if (f._doPlus) {
            if (f._doPlus>0) s.setf(std::ios_base::showpos); 
            else s.unsetf(std::ios_base::showpos); 
        }
        if (f._base != std::ios_base::dec) s.setf(std::ios_base::showbase);
    }

    inline void setup(std::ostream& os, const BoundForm<double>& bf)
    { setupFloat(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<long double>& bf)
    { setupFloat(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<float>& bf)
    { setupFloat(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<std::complex<double> >& bf)
    { setupFloat(os,bf.f); }

    inline void setup(
        std::ostream& os, const BoundForm<std::complex<long double> >& bf)
    { setupFloat(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<std::complex<float> >& bf)
    { setupFloat(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<int>& bf)
    { setupInt(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<short>& bf)
    { setupInt(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<long>& bf)
    { setupInt(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<unsigned int>& bf)
    { setupInt(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<unsigned short>& bf)
    { setupInt(os,bf.f); }

    inline void setup(std::ostream& os, const BoundForm<unsigned long>& bf)
    { setupInt(os,bf.f); }

    template <typename T>
    inline void setup(std::ostream& os, const BoundForm<T>& bf)
    { setupFloat(os,bf.f); }

    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, const BoundForm<T>& bf)
    {
        std::ostringstream s;
        setup(s,bf);
        s << bf.val;
        if (bf.f._nTrail>0) s << std::string(bf.f._nTrail,bf.f._trailCh);
        os << s.str();
        return os;
    }

}}}}

#endif
