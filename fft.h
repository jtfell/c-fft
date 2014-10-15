#include "complex.h"

complex* FFT_CooleyTukey(complex* x, int N, int N1, int N2);
complex* FFT_GoodThomas(complex* x, int N, int N1, int N2);
complex* DFT_naive(complex* x, int N);