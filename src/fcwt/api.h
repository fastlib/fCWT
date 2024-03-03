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



#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>
#include <complex>

#include <iostream>

#ifndef SINGLE_THREAD
    #include <omp.h>
#endif
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif
#include "fftw3.h"
#include <memory>
//check if avx is supported and include the header
#if defined(__AVX__)
    #include <immintrin.h>
    #define AVX
    union U256f {
        __m256 v;
        float a[8];
    };
#endif

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f
#define sqrt2PI               2.50662827463100050241576528f
#define IPI4                  0.75112554446f


#include "wavelet.h"
#include "scales.h"

namespace fcwt {
    class API {
        public:
        FCWT_LIBRARY_API API(Wavelet *pwav, int pthreads=1, bool puse_optimalization_schemes=false, bool puse_normalization=false):
            wavelet(pwav),
            threads(pthreads),
            use_optimalization_schemes(puse_optimalization_schemes),
            use_normalization(puse_normalization) {};

        void FCWT_LIBRARY_API create_FFT_optimization_plan(int maxsize, int flags) const;
        void FCWT_LIBRARY_API create_FFT_optimization_plan(int pmaxsize, std::string poptimizationflags);
        void FCWT_LIBRARY_API cwt(float *pinput, int psize, std::complex<float>* poutput, Scales *scales);
        void FCWT_LIBRARY_API cwt(std::complex<float> *pinput, int psize, std::complex<float>* poutput, Scales *scales);
        void FCWT_LIBRARY_API cwt(float *pinput, int psize, Scales *scales, std::complex<float>* poutput, int pn1, int pn2);
        void FCWT_LIBRARY_API cwt(float *pinput, int psize, std::complex<float>* poutput, Scales *scales, bool complexinput);
        void FCWT_LIBRARY_API cwt(std::complex<float> *pinput, int psize, Scales* scales, std::complex<float>* poutput, int pn1, int pn2);
        Wavelet *wavelet;

        private:

        void convolve(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, std::complex<float> *out, Wavelet *wav, int size, int newsize, float scale, bool lastscale);

        void fftbased(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, float *out, float* mother, int size, float scale, bool imaginary, bool doublesided);

        void fft_normalize(std::complex<float>* out, int size);

        void load_FFT_optimization_plan();

        void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float const *mother,
                                             float scale, int isize, bool imaginary, bool doublesided) const;

        static int find2power(const int n)
        {
            int m = 0;
            int m2 = 1 << m; /* 2 to the power of m */
            while (m2 - n < 0) {
                m++;
                m2 <<= 1; /* m2 = m2*2 */
            }
            return(m);
        }

        int threads;
        int size;
        float fs, f0, f1, fn;
        bool use_optimalization_schemes;
        bool use_normalization;
    };

}

#endif