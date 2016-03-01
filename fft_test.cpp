// Sample program to use fft functions

#include "fft.h"
#include <ctime>

#define max 100
#define min -100

using namespace shk;

int main(void)
{
    int N;

    cout << "N="; cin >> N;
    //srand((unsigned)time(NULL));
    srand((unsigned)82954);

    Complex* in  = new Complex[N];
    Complex* out = new Complex[N];

    for (int i=0; i<N; i++)
    {
        double r = (static_cast<double>(rand())/RAND_MAX)*(max-min)+min;
        in[i] = Complex(r,0.0);
    }

    fft_recursive(in,out,N);

    cout << "input\t\toutput\n";
    for (int i=0; i<N; i++)
    {
        cout << setprecision(14) << in[i] << "\t\t" << out[i] << endl;
    }

    fft_iterative(in,out,N);

    cout << "input\t\toutput\n";
    for (int i=0; i<N; i++)
    {
        cout << setprecision(14) << in[i] << "\t\t" << out[i] << endl;
    }

    return 0;
}
