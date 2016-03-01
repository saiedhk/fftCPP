/**********************************************************************
   Project: C++ Class for Complex Numbers

   Language: C++ 2007	   
   Author: Saied H. Khayat
   Date:   Nov 2011
   URL: https://github.com/saiedhk/ComplexCPP
   
   Copyright Notice: Free use of this library is permitted under
   the guidelines and in accordance with the MIT License (MIT).
   http://opensource.org/licenses/MIT

**********************************************************************/

#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
using std::ostream;

#ifndef SHK_COMPLEX
#define SHK_COMPLEX

namespace shk
{

// This class define complex numbers and basic complex arithmetic operations.
// All class functions are included in this header file as inline functions.
class Complex
{
    public:
        // constructor
        Complex(double=0.0, double=0.0);
        // operators
        Complex operator+(const Complex &) const;     // addition
        Complex operator-(const Complex &) const;     // subtraction
        Complex operator*(const Complex &) const;     // multiplication
        Complex operator/(const Complex &) const;     // division
        bool operator==(const Complex &) const;       // equality check
        bool operator!=(const Complex &) const;       // inequality check
        const Complex &operator=(const Complex &);    // assignment
        // friend functions
        friend Complex Conj(const Complex &);         // get complex conjugate
        friend double Re(const Complex &);            // get real part
        friend double Im(const Complex &);            // get imaginary part
        friend double Arg(const Complex &);           // get argument (phase)
        friend double Mod(const Complex &);           // get modulus (absolute value)
        friend Complex Exp(const Complex &);          // get exp(z)
        friend ostream &operator<<(ostream &, const Complex &);
    private:
        double real;  // real part
        double imag;  // imaginary part
};


//---------------------- inline functions ----------------------------
// Constructor
inline Complex::Complex(double r, double i) : real(r), imag(i) {};


// addition operator
inline Complex Complex::operator+(const Complex &z) const
{
    return Complex(real + z.real, imag + z.imag);
}


// subtraction operator
inline Complex Complex::operator -(const Complex &z) const
{
    return Complex(real - z.real, imag - z.imag);
}


// multiplication operator
inline Complex Complex::operator *(const Complex &z) const
{
    double r = real * z.real - imag * z.imag;
    double i = real * z.imag + imag * z.real;
    return Complex(r,i);
}


// division operator
inline Complex Complex::operator /(const Complex &z) const
{
    double denum = (z.real * z.real) + (z.imag * z.imag);
    assert(!(denum==0.0));
    double num_real = real * z.real + imag * z.imag;
    double num_imag = imag * z.real - real * z.imag;
    return Complex(num_real / denum , num_imag / denum);
}


// equality operator
inline bool Complex::operator==(const Complex &z) const
{
    if ( (real==z.real) && (imag==z.imag) )
        return true;
    else
        return false;
}


// inequality operator
inline bool Complex::operator!=(const Complex &z) const
{
    if ( (real != z.real) || (imag != z.imag) )
        return true;
    else
        return false;
}


// assignment operator
inline const Complex& Complex::operator=(const Complex & z)
{
    real = z.real;
    imag = z.imag;
    return *this;
}

// Exp(z)
inline Complex Exp(const Complex &z)
{
    double r = exp(z.real) * cos(z.imag);
    double i = exp(z.real) * sin(z.imag);
    return Complex(r,i);
}

// complex conjugate
inline Complex Conj(const Complex& z)
{
    return Complex(z.real,-z.imag);
}


inline ostream& operator<<(ostream &output, const Complex& z)
{
    output << "Z( " << z.real << " , " << z.imag << " )";
    //output << setiosflags(ios::showpos) << z.real << z.imag << "i";
    return output;
}


// useful friend functions
inline double Re(const Complex &z)  {return z.real;}
inline double Im(const Complex &z)  {return z.imag;}
inline double Mod(const Complex &z) {return sqrt(z.real * z.real + z.imag * z.imag);}
inline double Arg(const Complex &z) {return atan(z.imag/z.real);}

// useful constants
const Complex C_I(0.,1.);
const Complex C_ZERO(0.,0.);
const Complex C_ONE(1.,0.);

} // namespace shk
#endif
