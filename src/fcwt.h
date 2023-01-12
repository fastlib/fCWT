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
#pragma once

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
#include <math.h>

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
<<<<<<< HEAD
//check if avx is supported and include the header
#if defined(__AVX__)
    #include <immintrin.h>
    #define AVX

    union U256f {
        __m256 v;
        float a[8];
    };
#endif

#define PI                    3.14159265358979323846264338f
#define sqrt2PI               2.50662827463100050241576528f
#define IPI4                  0.75112554446f

using namespace std;

enum SCALETYPE {FCWT_LINSCALES,FCWT_LOGSCALES,FCWT_LINFREQS};

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

class Scales {
public:
    FCWT_LIBRARY_API Scales(Wavelet *pwav, SCALETYPE st, int fs, float f0, float f1, int fn);
=======
#include <immintrin.h>

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f
>>>>>>> parent of fe242c7 (fCWT 2.0 initial commit)

namespace fcwt {

    void precalculate_morlet(float* mother, float cf, int isize);
    void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float *mother, float scale, int isize);
    void FCWT_LIBRARY_API create_optimization_schemes(int maxsize, int threads, int optimizationflags);
    void load_optimization_schemes(bool use_optimalization_schemes, int size, int nthread);
    void main(float *Rinput,float *Routput,int *stboctave, int *endoctave, int *pnbvoice, int *pinputsize, float *pcenterfrequency, int nthreads, bool use_optimalization_schemes);
    void FCWT_LIBRARY_API cwt(float *input, int inputsize, float* output, int stboctave, int endoctave, int pnbvoice, float c0, int threads, bool use_optimalization_schemes);

<<<<<<< HEAD
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

    //Used in Python package
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
=======
>>>>>>> parent of fe242c7 (fCWT 2.0 initial commit)
}
