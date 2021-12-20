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

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <unistd.h>
#include <stdbool.h>

#include <iostream>
#include <omp.h>
#include <fftw3.h>
#include <complex.h>
#include <memory>
#include <immintrin.h>

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f

namespace fcwt {

    void precalculate_morlet(float* mother, float cf, int isize);
    void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float *mother, float scale, int isize);
    void create_optimization_schemes(int maxsize, int threads, int optimizationflags);
    void load_optimization_schemes(bool use_optimalization_schemes, int size, int nthread);
    void main(float *Rinput,float *Routput,int *stboctave, int *endoctave, int *pnbvoice, int *pinputsize, float *pcenterfrequency, int nthreads, bool use_optimalization_schemes);
    void cwt(float *input, int inputsize, float* output, int stboctave, int endoctave, int pnbvoice, float c0, int threads, bool use_optimalization_schemes);
    void cwt_test(float *input, int inputsize, int nscales, float* output);

}
