/**********************************************************************
   Project: C++ Functions for Fast Fourier Transform

   Language: C++ 2007	   
   Author: Saied H. Khayat
   Date:   Nov 2011
   URL: https://github.com/saiedhk/fftCPP
   
   Copyright Notice: Free use of this library is permitted under
   the guidelines and in accordance with the MIT License (MIT).
   http://opensource.org/licenses/MIT

**********************************************************************/

#ifndef SHK_FFT
#define SHK_FFT

#include <cassert>
using namespace std;

// requires shk_complex.h header file
#include "shk_complex.h"

namespace shk
{

// function used by fft()
unsigned bit_reverse(unsigned index, int width);

const double PI = 4 * atan(1.0);

// the fft_recursive() function
void fft_recursive
(
    const Complex in[],  // input array of complex numbers
    Complex out[],       // output array of complex numbers 
    int length           // length of array (must be a power of 2)
);


// the fft_iterative() function
void fft_iterative
(
    const Complex in[],  // input array of complex numbers
    Complex out[],       // output array of complex numbers 
    int length           // length of array (must be a power of 2)
);

}
#endif // SHK_FFT
