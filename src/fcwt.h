//
//  fcwt.h
//  fCWT
//
//  Created by Lukas Arts on 21/12/2020.
//  Copyright © 2021 Lukas Arts.
/*Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#ifndef FCWT_H
#define FCWT_H

#ifdef _WIN32
  #ifdef FCWT_LIBRARY_DLL_BUILDING
    #define FCWT_LIBRARY_API __declspec(dllexport)
  #else
    #if FCWT_LIBRARY_DLL
      #define FCWT_LIBRARY_API __declspec(dllimport)
    #else /* static or header-only library on Windows */
      #define FCWT_LIBRARY_API
    #endif
  #endif
#else /* Unix */
  #define FCWT_LIBRARY_API
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <stdbool.h>
#include <vector>
#include <chrono>
#include <cassert>
#include <math.h>
#include <complex>

#include <iostream>
#include <sstream>
#include <omp.h>
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif
#include <fftw3.h>
#include <memory>
//check if avx is supported and include the header
#if defined(__AVX__)
    #include <immintrin.h>
    #define AVX
#endif

#define PI                    3.14159265358979323846264338f
#define sqrt2PI               2.50662827463100050241576528f
#define IPI4                  0.75112554446f

using namespace std;

enum SCALETYPE {FCWT_LINSCALES,FCWT_LOGSCALES,FCWT_LINFREQS};

union U256f {
	__m256 v;
	float a[8];
};

class Wavelet {
public:
    Wavelet() {};
    virtual void generate(float* real, float* imag, int size, float scale) { printf("ERROR [generate time complex]: Override this virtual class"); };
    virtual void generate(int size) { printf("ERROR [generate freq]: Override this virtual class"); };
    virtual int getSupport(float scale) { printf("ERROR [getsupport]: Override this virtual class"); return 0; };
    
    int width;
    float four_wavelen;
    bool imag_frequency, doublesided;
    float *mother;
};

class Morlet : public Wavelet {
public:
    FCWT_LIBRARY_API Morlet(float bandwidth); //frequency domain
    ~Morlet() { free(mother); };
    
    void generate(int size); //frequency domain
    void generate(float* real, float* imag, int size, float scale); //time domain
    int getSupport(float scale) { return (int)(fb*scale*3.0f); };
    float fb;
    
private:
    float ifb, fb2;
};

// class DOG : public Wavelet {
// public:
//     FCWT_LIBRARY_API DOG(int order);
//     ~DOG() { free(mother); };
    
//     void generate(int size); //frequency domain
//     void generate(float* real, float* imag, int size, float scale); //time domain
//     int getSupport(float scale) { return (int)(fb*scale*3.0f+1); };
//     int order;
    
// private:
//     float fb, ifb, fb2;
// };

// class Paul : public Wavelet {
// public:
//     FCWT_LIBRARY_API Paul(int porder);
//     ~Paul() { free(mother); };
    
//     void generate(int size); //frequency domain
//     void generate(float* real, float* imag, int size, float scale); //time domain
//     int getSupport(float scale) { return (int)(fb*scale*fmax(3.0,(6.0-order*1.0))); };
//     int order;
    
// private:
//     float fb;
// };

class Scales {
public:
    FCWT_LIBRARY_API Scales(Wavelet *pwav, SCALETYPE st, int fs, float f0, float f1, int fn);

    void FCWT_LIBRARY_API getScales(float *pfreqs, int pnf);
    void FCWT_LIBRARY_API getFrequencies(float *pfreqs, int pnf);

    float* scales;
    int fs;
    float fourwavl;
    int nscales;

private:
    void calculate_logscale_array(float base, float four_wavl, int fs, float f0, float f1, int fn);
    void calculate_linscale_array(float four_wavl, int fs, float f0, float f1, int fn);
    void calculate_linfreq_array(float four_wavl, int fs, float f0, float f1, int fn);
};

class FCWT {
public:
    FCWT_LIBRARY_API FCWT(Wavelet *pwav, int pthreads, bool puse_optimalization_schemes, bool puse_normalization):  
        wavelet(pwav), 
        threads(pthreads), 
        use_optimalization_schemes(puse_optimalization_schemes),
        use_normalization(puse_normalization) {};

    void FCWT_LIBRARY_API create_FFT_optimization_plan(int pmaxsize, int poptimizationflags);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, complex<float>* poutput, Scales *scales);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, Scales *scales, complex<float>** poutput, int* pnoutput);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pnoutput);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2, float *pfreqs, int pnf);

    Wavelet *wavelet;
    
private:
    void cwt_static(float *pinput, int psize, float* poutput, float* scales);
    void cwt_dynamic(float *pinput, int psize, float* poutput, float* scales);
    void convolve(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, complex<float> *out, Wavelet *wav, int size, int newsize, float scale, bool lastscale);
    void convolve(float* in, complex<float> *out, Wavelet *wav, float scale);
    void fftbased(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, float *out, float* mother, int size, float scale, bool imaginary, bool doublesided);
    void firbased(float* in, float *out, Wavelet *wav, float scale);
    void fft_normalize(complex<float>* out, int size);
    void main(float *Rinput,float *Routput);
    void load_FFT_optimization_plan();
    void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float const *mother, float scale, int isize, bool imaginary, bool doublesided);
    
    void calculate_logscale_array(float base, float four_wavl, float *scales);
    void calculate_linscale_array(float four_wavl, float *scales);
    
    int threads;
    int size;
    float fs, f0, f1, fn;
    bool use_optimalization_schemes;
    bool use_normalization;
};

inline int find2power(int n)
{
    int m, m2;
        
    m = 0;
    m2 = 1<<m; /* 2 to the power of m */
    while (m2-n < 0) {
        m++;
        m2 <<= 1; /* m2 = m2*2 */
    }
    return(m);
}

inline double factorial(int N) {
    static const double fact[41] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000, 51090942171709440000.0, 1124000727777607680000.0, 25852016738884976640000.0, 620448401733239439360000.0, 15511210043330985984000000.0, 403291461126605635584000000.0, 10888869450418352160768000000.0, 304888344611713860501504000000.0, 8841761993739701954543616000000.0, 265252859812191058636308480000000.0, 8222838654177922817725562880000000.0, 263130836933693530167218012160000000.0, 8683317618811886495518194401280000000.0, 295232799039604140847618609643520000000.0, 10333147966386144929666651337523200000000.0, 371993326789901217467999448150835200000000.0, 13763753091226345046315979581580902400000000.0, 523022617466601111760007224100074291200000000.0, 20397882081197443358640281739902897356800000000.0, 815915283247897734345611269596115894272000000000.0 };
    return fact[N];
}

inline float gamma_dog(int N) {
    static const float gamma[11] = {0.751126, 1.0623, 0.8673, 0.5485, 0.2932, 0.1382, 0.0589, 0.0231, 0.0084, 0.0029, 0.0009};
    return gamma[N];
}

#endif