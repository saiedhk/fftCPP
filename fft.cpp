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

#include "fft.h"

namespace shk
{


/**
  Computes FFT of an input array of complex numbers (recursive implementation)
  @param in: an input array of Complex type
  @param out: an array holding the FFT of the input array
  @param length: the length of input and output arrays; must be a power of two
*/
void fft_recursive
    (
    const Complex in[],  // input array
    Complex out[],       // output array
    int length           // array length
    )
{

    if (length==2) // end of recursion (compute FFT-2)
    {
        out[0] = in[0] + in[1];
        out[1] = in[0] - in[1];
    }
    else // if (length != 2)
    {
        // fatal exit if length is too short.
        assert(length > 2);

        int Q = 0; // There are Q stages for an N-point FFT.
        int N = 1; // N=2^Q  (N-point FFT is to be computed.)

        // compute actual Q and N based on length argument provided
        while (N<length) { N *= 2; Q++; }
        
        // fatal exit if length is not a power of 2.
        assert(N==length);

        // construct the array of twiddle factors (WN^k values)
        Complex  WN = Exp(Complex(0.0,-2*PI/N)); // WN = e^{-2*PI*j/N}
        Complex* WNk = new Complex[N/2];         // array to hold all WN^k

        // compute all required WN^k values and store in WNk array
        WNk[0] = Complex(1.,0.);
        for (int k=1; k<N/2; k++) { WNk[k] = WN * WNk[k-1]; }

        // allocate four arrays for scratch work
        Complex* XL = new Complex[N/2];
        Complex* XH = new Complex[N/2];
        Complex* G  = new Complex[N/2];
        Complex* H  = new Complex[N/2];

        // copy input array into array X (after bit reversing the indices)
        for (int i=0; i<N/2; i++)
        {
            XL[i] = in[i*2];
            XH[i] = in[i*2+1];
        }

        // recursive calls
        fft_recursive(XL,G,N/2);
        fft_recursive(XH,H,N/2);

        for (int i=0; i<N/2; i++)
        {
            Complex temp = WNk[i] * H[i];
            out[i] = G[i] + temp;
            out[(N/2)+i] = G[i] - temp;
        }

        // deallocate internal arrays
        delete [] WNk;
        delete [] XL;
        delete [] XH;
        delete [] G;
        delete [] H;
    }

    /* print for debugging */
    // cout << "----------\n";
    // for (int i=0; i<length; i++)
    // {
        // cout << in[i] << "\t" << out[i] << endl;
    // }
    // cout << "----------\n";

} // end fft_recursive()


//-----------------------------------------------------------------------------
/**
  Computes FFT of an input array of complex numbers (iterative implementation)
  @param in: an input array of Complex type
  @param out: an array holding the FFT of the input array
  @param length: the length of input and output arrays; must be a power of two
*/
void fft_iterative
    (
    const Complex in[],  // input array
    Complex out[],       // output array
    int length           // array length
    )
{
    // fatal exit if length is too short.
    assert(length >= 2);

    int Q = 0; // There are Q stages for an N-point FFT.
    int N = 1; // N=2^Q  (N-point FFT is to be computed.)

    // compute actual Q and N based on length argument provided
    while (N<length) { N *= 2; Q++; }
    
    // throw an exception if length is not a power of 2.
    assert(N==length); 
    
    // construct the array of twiddle factors (WN^k values)
    Complex  WN = Exp(Complex(0.0,-2*PI/N)); // WN = e^{-2*PI*j/N}
    Complex* WNk = new Complex[N/2];         // array to hold all WN^k

    // compute all required WN^k values and store in WNk array
    WNk[0] = Complex(1.,0.);
    for (int k=1; k<N/2; k++) { WNk[k] = WN * WNk[k-1]; }

    // allocate two arrays X and Y for scratch work
    Complex* X = new Complex[N];
    Complex* Y = new Complex[N];

    // copy input array into array X (after bit reversing the indices)
    for (int i=0; i<N; i++) { X[i] = in[bit_reverse(i,Q)]; }

    int p = 1;
    int k = 1;
    int m = N/2;
    for (int q=0; q<Q; q++) // for each stage q=0 to Q-1
    {
        p *= 2;
        for (int i=0; i<N; i+=p)
        {
            for (int j=0; j<k; j++)
            {
                int n = i+j;
                // butterfly operation (next 3 lines)
                Complex XWNk = X[n+k] * WNk[j*m];
                Y[n]   = X[n] + XWNk;
                Y[n+k] = X[n] - XWNk;
            }
        }
        k *= 2;
        m /= 2;

        if (q<Q-1) // if not last stage, copy Y to X
            { for (int i=0; i<N; i++) X[i] = Y[i]; }
        else // if last stage, copy Y to output array
            { for (int i=0; i<N; i++) out[i] = Y[i]; }

    }

    // deallocate internal arrays
    delete [] X;
    delete [] Y;
    delete [] WNk;

} // end fft_iterative()



/**
  Computes the bit-reversed version of an integer
  @param index: an input integer; must be unsigned
  @param width: the number of LSB bits in index to be included in reversion
  @return bit-reversed version of index
*/
unsigned bit_reverse(unsigned index, int width)
{
    unsigned result=0;

    for (int i=0; i<width; i++)
    {
        result <<= 1;
        result |= (index & 0x00000001);
        index >>= 1;
    }
    return result;
}


} // namespace shk
