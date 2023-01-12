//
//  fcwt.cpp
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

//
//  fcwt.cpp
//  fCWT-testing
//
//  Created by Lukas Arts on 21/12/2020.
//  Copyright © 2020 Lukas Arts.
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

#include "fcwt.h"

Morlet::Morlet(float bandwidth) {
    four_wavelen = 0.9876f;
    fb = bandwidth;
    fb2 = 2.0f*fb*fb;
    ifb = 1.0f/fb;
    imag_frequency = false;
    doublesided = false;
    mother = NULL;
}

void Morlet::generate(int size) {
    //Frequency domain, because we only need size. Default scale is always 2;
    width = size;
    
    float tmp1;
    float toradians = (2*PI)/(float)size;
    float norm = sqrt(2*PI)*IPI4;
    
    mother = (float*)malloc(sizeof(float)*width);
    
    //calculate array
    for(int w = 0; w < width; w++) {
        tmp1 = (2.0f * ((float)w * toradians) * fb - 2.0f*PI*fb);
        tmp1 = -(tmp1 * tmp1)/2;
        mother[w] = (norm*exp(tmp1));
    }
}
void Morlet::generate(float* real, float* imag, int size, float scale) {
    //Time domain because we know size from scale
    float tmp1, tmp2;
    width = getSupport(scale);
    float norm = (float)size * ifb * IPI4;
    
    //cout << scale << " [";
    for(int t=0; t < width*2+1; t++) {
        tmp1 = (float)(t - width)/scale;
        tmp2 = exp(-(tmp1*tmp1)/(fb2));
        
        real[t] = norm*tmp2*cos(tmp1*2.0f*PI)/scale;
        imag[t] = norm*tmp2*sin(tmp1*2.0f*PI)/scale;
        //cout << real[t]*real[t]+imag[t]*imag[t] << ",";
    }
    //cout << "]" << endl;
}

// //Frequency domain
// DOG::DOG(int porder) {
//     four_wavelen = (2*PI)/(sqrt((float)porder + 0.5));
//     fb = 1.0f;
//     fb2 = 2.0f*fb*fb;
//     ifb = 1.0f/fb;
//     order = porder;
//     doublesided = true;
    
//     if(order % 2 == 0) {
//         imag_frequency = false;
//     } else {
//         imag_frequency = true;
//     }
// }

// void DOG::generate(int size) {
    
//     width = size;
//     float toradians = (2*PI)/(float)size;
//     float norm = sqrt(2*PI);
    
//     mother = (float*)malloc(sizeof(float)*width);
    
//     float sign;
//     if (order % 4 == 0 || order % 4 == 1) {
//         sign = -1;
//     }
//     else {
//         sign = 1;
//     }
    
//     //float sum = 0;
//     for(int w=0; w < width; w++) {
//         float exponent = -pow((2.0f*((float)w * toradians)),2.0)*0.5f;
//         mother[w] = norm * gamma_dog(order) * sign * pow(2.0f*((float)w * toradians),order)*exp(exponent);
//         //sum += complex[w]*complex[w];
//     }
//     //cout << "sum: " << sum << endl;
    
// }

// void DOG::generate(float* real, float* imag, int size, float scale) {
    
//     width = getSupport(scale);
    
//     float tmp1, tmp2, tmp4, exponent;
//     float sign = -1.0f;
//     if((order+1) % 2 ==0) {
//         sign = 1.0f;
//     }
//     float constant = ((float)size*(sign*gamma_dog(order)))/scale;
    
//     if (order == 1) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * -tmp1 * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 2) {
//         float sum = 0;
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * (-1.0 + tmp2) * exp(exponent));
//             imag[t] = 0.0f;
//             sum += real[t]*real[t];
//         }
//     } else if (order == 3) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * -tmp1 * (3.0 - tmp2) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 4) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * (tmp2*tmp2 - 6.0*tmp2 + 3.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 5) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * -tmp1 * (tmp2*tmp2 - 10.0*tmp2 + 15.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 6) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * (tmp2*tmp2*tmp2 - 15.0*tmp2*tmp2 + 45.0*tmp2 - 15.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 7) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * -tmp1 * (tmp2*tmp2*tmp2 - 21.0*tmp2*tmp2 + 105.0*tmp2 - 105.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 8) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             tmp4 = tmp2*tmp2;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * (tmp4*tmp4 - 28.0 * tmp4*tmp2 + 210.0 * tmp4 - 420.0 * tmp2 + 105.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 9) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             tmp4 = tmp2*tmp2;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * -tmp1 * (tmp4*tmp4 - 36.0 * tmp4*tmp2 + 378.0 * tmp4 - 1260.0 * tmp2 + 945.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else if (order == 10) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = tmp1*tmp1;
//             tmp4 = tmp2*tmp2;
//             exponent = -(tmp2)/fb2;
//             real[t] = (float)(constant * (tmp4*tmp4*tmp2 - 45.0 * tmp4*tmp4 + 630*tmp4*tmp2 - 3150.0*tmp4 + 4725.0*tmp2 - 945.0) * exp(exponent));
//             imag[t] = 0.0f;
//         }
//     } else {
//         cerr << "ERROR: DOG's are only defined for orders <= 10\n" << endl;
//     }
// }

// //Frequency domain
// Paul::Paul(int porder) {
//     order = porder;
//     four_wavelen = (4*PI)/(2*order+1);
//     fb = 1.0f;
//     imag_frequency = false;
//     doublesided = false;
// }

// void Paul::generate(int size) {
    
//     width = size;
//     float norm = sqrt(2*PI)*(pow(2.0,(double)order) / sqrt((double)(order*factorial(2 * order - 1))));
//     float toradians = (2*PI)/(float)width;
    
//     mother = (float*)malloc(sizeof(float)*width);
    
//     for (int w = 0; w < width; ++w) {
//         float temp = 2.0f * w * toradians;
//         mother[w] = norm * pow(temp,(double)order) * exp(-temp);
//     }
// }
// void Paul::generate(float* real, float* imag, int size, float scale) {
//     width = getSupport(scale);
    
//     float norm = ((float)size * (pow(2.0,(float)order)*(float)factorial(order)) / sqrt((float)(PI*(float)factorial(2 * order))));
//     float tmp1, tmp2, tmp22, tmp24;
    
//     if(order==1) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),2.0);
//             real[t] = norm*(2.0*tmp1)/(tmp2*scale);
//             imag[t] = -norm*(1.0 - tmp1*tmp1)/(tmp2*scale);
//         }
//     } else if (order==2) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),3.0);
//             tmp22 = tmp1*tmp1;
//             real[t] = norm*(-1.0 + 3.0*tmp22)/(tmp2*scale);
//             imag[t] = -norm*(-3.0*tmp1 + tmp22*tmp1)/(tmp2*scale);
//         }
//     } else if (order==3) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),4.0);
//             tmp22 = tmp1*tmp1;
//             real[t] = norm*(4.0 - 4.0*tmp22*tmp1)/(tmp2*scale);
//             imag[t] = -norm*(-1.0 + 6.0*tmp22 - tmp22*tmp22)/(tmp2*scale);
//         }
//     } else if (order==4) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),5.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(1.0 - 10.0*tmp22 + 5.0*tmp24)/(tmp2*scale);
//             imag[t] = norm*(5.0*tmp1 -10.0*tmp22*tmp1 + tmp24*tmp1)/(tmp2*scale);
//         }
//     } else if (order==5) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),6.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(-6.0*tmp1 + 20.0*tmp22*tmp1 - 6.0*tmp24*tmp1)/(tmp2*scale);
//             imag[t] = norm*(1.0 - 15.0*tmp22 + 15.0*tmp24 - tmp24*tmp22)/(tmp2*scale);
//         }
//     } else if (order==6) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),7.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(-1.0 + 21.0*tmp22 - 35.0*tmp24 + 7.0*tmp24*tmp22)/(tmp2*scale);
//             imag[t] = norm*(-7.0*tmp1 + 35.0*tmp22*tmp1 - 21.0*tmp24*tmp1 + tmp24*tmp22*tmp1)/(tmp2*scale);
//         }
//     } else if (order==7) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),8.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(8.0*tmp1 - 56.0*tmp22*tmp1 + 56.0*tmp24*tmp1 - 8.0*tmp24*tmp22*tmp1)/(tmp2*scale);
//             imag[t] = norm*(-1.0 + 28.0*tmp22 - 70.0*tmp24 + 28.0*tmp24*tmp22 - tmp24*tmp24)/(tmp2*scale);
//         }
//     } else if (order==8) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),9.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(1.0 - 36.0*tmp22 + 126.0*tmp24 - 84.0*tmp24*tmp22 + 9.0*tmp24*tmp24)/(tmp2*scale);
//             imag[t] = norm*(-9.0*tmp1 - 84.0*tmp22*tmp1 + 126.0*tmp24*tmp1 - 36.0*tmp24*tmp22*tmp1 + tmp24*tmp24*tmp1)/(tmp2*scale);
//         }
//     } else if (order==9) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),10.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(-10.0*tmp1 + 120.0*tmp22*tmp1 - 252.0*tmp24*tmp1 + 120.0*tmp24*tmp22*tmp1 - 10.0*tmp24*tmp24*tmp1)/(tmp2*scale);
//             imag[t] = norm*(1.0 - 45.0*tmp22 + 210.0*tmp24 - 210.0*tmp24*tmp22 + 45.0*tmp24*tmp24 - tmp24*tmp24*tmp22)/(tmp2*scale);
//         }
//     } else if (order==10) {
//         for(int t=0; t < width*2+1; t++) {
//             tmp1 = (float)(t - width)/scale;
//             tmp2 = pow((1.0 + tmp1*tmp1),11.0);
//             tmp22 = tmp1*tmp1;
//             tmp24 = tmp22*tmp22;
//             real[t] = norm*(-1.0 + 55.0*tmp22 - 330*tmp24 + 462*tmp24*tmp22 - 165.0*tmp24*tmp24 + 11.0*tmp24*tmp24*tmp22)/(tmp2*scale);
//             imag[t] = norm*(-11.0*tmp1 + 165.0*tmp22*tmp1 - 462.0*tmp24*tmp1 + 330.0*tmp24*tmp22*tmp1 - 55.0*tmp24*tmp24*tmp1 + tmp24*tmp24*tmp22*tmp1)/(tmp2*scale);
//         }
//     } else {
//         cerr << "ERROR: PAUL is only defined for orders <= 10\n" << endl;
//     }
// }


//==============================================================//
//================== Scales =====================================//
//==============================================================//

Scales::Scales(Wavelet *wav, SCALETYPE st, int afs, float af0, float af1, int afn) {
    
    fs = afs;
    scales = (float*)malloc(afn*sizeof(float));
    fourwavl = wav->four_wavelen;
    nscales = afn;

    if(st==SCALETYPE::FCWT_LOGSCALES)
        calculate_logscale_array(2.0f, wav->four_wavelen, afs, af0, af1, afn);
    else if(st==SCALETYPE::FCWT_LINSCALES)
        calculate_linscale_array(wav->four_wavelen, afs, af0, af1, afn);
    else 
        calculate_linfreq_array(wav->four_wavelen, afs, af0, af1, afn);

}

void Scales::getScales(float *pfreqs, int pnf) { 
    for(int i=0;i<pnf;i++) { 
        pfreqs[i]=scales[i]; 
    }; 
};

void Scales::getFrequencies(float *pfreqs, int pnf) { 
    for(int i=0;i<pnf;i++) { 
        pfreqs[i]=((float)fs)/scales[i]; 
    }; 
};

void Scales::calculate_logscale_array(float base, float four_wavl, int fs, float f0, float f1, int fn) {
    
    //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    float nf0 = f0;
    float nf1 = f1;
    float s0 = (fs/nf1);
    float s1 = (fs/nf0);
    
    //Cannot pass the nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs/2));

    float power0 = log(s0)/log(base);
    float power1 = log(s1)/log(base);
    float dpower = power1-power0;

    for(int i=0; i<fn; i++) {
        float power = power0 + (dpower/(fn-1))*i;
        scales[i] = pow(base,power);
    }
}

void Scales::calculate_linfreq_array(float four_wavl, int fs, float f0, float f1, int fn) {
    
    float nf0 = f0;
    float nf1 = f1;
    //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    
    //Cannot pass the nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs/2));
    float df = nf1-nf0;

    for(int i=0; i<fn; i++) {
        scales[fn-i-1] = (((float)fs)/(nf0 + (df/fn)*(float)i));
    }
}

void Scales::calculate_linscale_array(float four_wavl, int fs, float f0, float f1, int fn) {
    
    float nf0 = f0;
    float nf1 = f1;
    //If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    float s0 = fs/nf1;
    float s1 = fs/nf0;
    
    //Cannot pass the nyquist frequency
    assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs/2));
    float ds = s1-s0;

    for(int i=0; i<fn; i++) {
        scales[i] = (s0 + (ds/fn)*i);
    }
}


//==============================================================//
//================== FCWT =====================================//
//==============================================================//

void FCWT::daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float const *mother, float scale, int isize, bool imaginary, bool doublesided)
{
    float isizef = ((float)(isize));
    float endpointf = fmin(isizef/2.0,((isizef*2.0)/scale));
    float step = (scale/2.0);
    int endpoint = ((int)endpointf);
    int endpoint4 = endpoint>>2;

    #ifdef AVX
        //has avx instructions
        __m256* O8 = (__m256*)output;
        __m256* I8 = (__m256*)input;
        __m256 step4 = _mm256_set1_ps(step);
        __m256 offset = _mm256_set_ps(3,3,2,2,1,1,0,0);
        __m256 maximum = _mm256_set1_ps(isizef-1);
        
        int athreads = min(threads,max(1,endpoint4/16));
        int batchsize = (endpoint4/athreads);
        int s4 = (isize>>2)-1;

        #pragma omp parallel for
        for(int i=0; i<athreads; i++) {
            int start = batchsize*i;
            int end = batchsize*(i+1);
            
            for(int q4=start; q4<end; q4++) {
                float q = (float)q4*4;

                __m256 qq = _mm256_set1_ps(q);
                
                U256f tmp = {_mm256_min_ps(maximum,_mm256_mul_ps(step4,_mm256_add_ps(qq,offset)))};
                //U256f tmp = {_mm256_mul_ps(step4,_mm256_add_ps(qq,offset))};

                __m256 wav = _mm256_set_ps(
                                        mother[(int)tmp.a[7]]*(1-2*imaginary),
                                        mother[(int)tmp.a[6]],
                                        mother[(int)tmp.a[5]]*(1-2*imaginary),
                                        mother[(int)tmp.a[4]],
                                        mother[(int)tmp.a[3]]*(1-2*imaginary),
                                        mother[(int)tmp.a[2]],
                                        mother[(int)tmp.a[1]]*(1-2*imaginary),
                                        mother[(int)tmp.a[0]]);
                
                if(imaginary) {
                    __m256 tmp2 = _mm256_mul_ps(I8[q4],wav);
                    O8[q4] = _mm256_shuffle_ps(tmp2, tmp2, 177);
                } else {
                    O8[q4] = _mm256_mul_ps(I8[q4],wav);
                }
            }

            if(doublesided) {
                for(int q4=start; q4<end; q4++) {
                    float q = (float)(q4*4);
                    
                    __m256 qq = _mm256_set1_ps(q);
                    U256f tmp = {_mm256_mul_ps(step4,_mm256_add_ps(qq,offset))};
                    
                    __m256 wav = _mm256_set_ps(
                                            mother[(int)tmp.a[0]]*(1-2*imaginary),
                                            mother[(int)tmp.a[1]],
                                            mother[(int)tmp.a[2]]*(1-2*imaginary),
                                            mother[(int)tmp.a[3]],
                                            mother[(int)tmp.a[4]]*(1-2*imaginary),
                                            mother[(int)tmp.a[5]],
                                            mother[(int)tmp.a[6]]*(1-2*imaginary),
                                            mother[(int)tmp.a[7]]);
                    
                    if(imaginary) {
                        __m256 tmp2 = _mm256_mul_ps(I8[s4-q4],wav);
                        O8[s4-q4] = _mm256_shuffle_ps(tmp2, tmp2, 177);
                    } else {
                        O8[s4-q4] = _mm256_mul_ps(I8[s4-q4],wav);
                    }
                }
            }
        }
    #else
        int athreads = min(threads,max(1,endpoint/16));
        int batchsize = (endpoint/athreads);
        float maximum = isizef-1;
        int s1 = isize-1;

        #pragma omp parallel for
        for(int i=0; i<athreads; i++) {
            int start = batchsize*i;
            int end = batchsize*(i+1);
            
            for(int q1=start; q1<end; q1++) {
                float q = (float)q1;
                float tmp = min(maximum,step*q);
                
                output[q1][0] = input[q1][0]*mother[(int)tmp];
                output[q1][1] = input[q1][1]*mother[(int)tmp]*(1-2*imaginary);
            }

            if(doublesided) {
                for(int q1=start; q1<end; q1++) {
                    float q = (float)q1;
                    float tmp = min(maximum,step*q);
                    
                    output[s1-q1][0] = input[s1-q1][0]*mother[(int)tmp]*(1-2*imaginary);
                    output[s1-q1][1] = input[s1-q1][1]*mother[(int)tmp];
                }
            }
        }

    #endif
    
    return;
}

void FCWT::create_FFT_optimization_plan(int maxsize, int flags) {
    
    int nt = find2power(maxsize);
    
    if(nt <= 10) {
        std::cerr << "Maxsize is too small (<=1024)... please use a larger number" << std::endl;
    }
    
    for(int i=11; i<=nt; i++) {
        int n = 1 << i;
        
        float *dat = (float*)malloc(sizeof(float)*n);
        fftwf_complex *O1 = fftwf_alloc_complex(n);
        fftwf_complex *out = fftwf_alloc_complex(n);
        
        omp_set_num_threads(threads);
        std::cout << "Threads:" << omp_get_max_threads() << std::endl;
        
        fftwf_init_threads();
        fftwf_plan_with_nthreads(omp_get_max_threads());
        fftwf_plan p_for;
        fftwf_plan p_back;
        
        char file_for[50];
        sprintf(file_for, "n%d_t%d.wis", n, threads);
        
        std::cout << "Calculating optimal scheme for forward FFT with N:" << n << std::endl;
        p_for = fftwf_plan_dft_r2c_1d(n, dat, O1, flags);
        
        std::cout << "Calculating optimal scheme for backward FFT with N:" << n << std::endl;
        p_back = fftwf_plan_dft_1d(n, O1, out, FFTW_BACKWARD, flags);
        
        fftwf_export_wisdom_to_filename(file_for);
        
        free(dat);
        fftwf_free(O1);
        fftwf_free(out);
        
        std::cout << "Optimization schemes for N: " << n << " have been calculated. Next time you use fCWT it will automatically choose the right optimization scheme based on number of threads and signal length." << std::endl;
    }
}

void FCWT::load_FFT_optimization_plan() {
    const int nt = find2power(size);
    const int newsize = 1 << nt;
    
    if(use_optimalization_schemes) {
        if(newsize <= 1024) {
            std::cout << "Inputsize is too small (N <= 1024) to use optimization." << std::endl;
            return;
        }
        
        char file_for[50];
        sprintf(file_for, "n%d_t%d.wis", newsize, threads);
        
        if(!fftwf_import_wisdom_from_filename(file_for)) {
            std::cout << "WARNING: Optimization scheme '" << file_for << "' was not found, fallback to calculation without optimization." << std::endl;
        }
    }
}

//Convolve in time domain using a single wavelet
void FCWT::convolve(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, complex<float> *out, Wavelet *wav, int size, int newsize, float scale, bool lastscale) {
    
    if(lastscale) {
        #ifdef _WIN32
            fftwf_complex *lastscalemem = (fftwf_complex*)_aligned_malloc(newsize*sizeof(fftwf_complex), 32);
        #else
            fftwf_complex *lastscalemem = (fftwf_complex*)aligned_alloc(32, newsize*sizeof(fftwf_complex));
        #endif
        memset(lastscalemem,0,sizeof(fftwf_complex)*newsize);
        
        fftbased(p, Ihat, O1, (float*)lastscalemem, wav->mother, newsize, scale, wav->imag_frequency, wav->doublesided);
        if(use_normalization) fft_normalize((complex<float>*)lastscalemem, newsize);
        memcpy(out, (complex<float>*)lastscalemem, sizeof(complex<float>)*size);
    } else {
        if(!out) {
            std::cout << "OUT NOT A POINTER" << std::endl;
        }
        fftbased(p, Ihat, O1, (float*)out, wav->mother, newsize, scale, wav->imag_frequency, wav->doublesided);
        if(use_normalization) fft_normalize(out, newsize);
    }
}

void FCWT::fftbased(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, float *out, float* mother, int size, float scale, bool imaginary, bool doublesided) {
    
    void *pt = out;
    
    //Perform daughter wavelet generation and multiplication with the Fourier transformed input signal
    daughter_wavelet_multiplication(Ihat,O1,mother,scale,size,imaginary,doublesided);

    std::size_t space = 16;
    std::align(16,sizeof(fftwf_complex),pt,space);
    
    fftwf_execute_dft(p,O1,(fftwf_complex*)pt);
}

void FCWT::fft_normalize(complex<float>* out, int size) {

    int nbatch = threads;
    int batchsize = (int)ceil((float)size/((float)threads));
    
    #pragma omp parallel for
    for(int i=0; i<nbatch; i++) {
        int start = batchsize*i;
        int end = min(size,batchsize*(i+1));
        
        for(int i8=start; i8<end; i8++) {
            out[i8] = out[i8] / (float)size;
        }
    }
}

void FCWT::cwt(float *pinput, int psize, complex<float>* poutput, Scales *scales) {
    
    fftwf_complex *Ihat, *O1;
    size = psize;
    
    //Find nearest power of 2
    const int nt = find2power(size);
    const int newsize = 1 << nt;
    
    //Initialize input and mother wavelet memory
    float* input = (float*)calloc(newsize,sizeof(float));

    //Initialize intermediate result
    #ifdef _WIN32
        Ihat = (fftwf_complex*)_aligned_malloc(newsize*sizeof(fftwf_complex), 32);
        O1 = (fftwf_complex*)_aligned_malloc(newsize*sizeof(fftwf_complex), 32);
    #else
        Ihat = (fftwf_complex*)aligned_alloc(32, newsize*sizeof(fftwf_complex));
        O1 = (fftwf_complex*)aligned_alloc(32, newsize*sizeof(fftwf_complex));
    #endif
    
    //Copy input to new input buffer
    memcpy(input,pinput,sizeof(float)*size);
    memset(Ihat,0,sizeof(fftwf_complex)*newsize);
    memset(O1,0,sizeof(fftwf_complex)*newsize);
    
    // //Initialize FFTW plans
    omp_set_num_threads(threads);
    
    // //Initialize FFTW plans
    fftwf_init_threads();
    
    fftwf_plan_with_nthreads(omp_get_max_threads());
    fftwf_plan pinv;
    fftwf_plan p;
    
    // //Load optimization schemes if necessary
    load_FFT_optimization_plan();

    // //Perform forward FFT on input signal
    p = fftwf_plan_dft_r2c_1d(newsize, &input[0], Ihat, FFTW_ESTIMATE);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    
    pinv = fftwf_plan_dft_1d(newsize, O1, (fftwf_complex*)poutput, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    //Generate mother wavelet function
    wavelet->generate(newsize);
    
    for(int i=1; i<(newsize>>1); i++) {
        Ihat[newsize-i][0] = Ihat[i][0];
        Ihat[newsize-i][1] = -Ihat[i][1];
    }
    
    complex<float> *out = poutput;
    
    for(int i = 0; i < scales->nscales; i++) {
        //FFT-base convolution in the frequency domain
        convolve(pinv, Ihat, O1, out, wavelet, size, newsize, scales->scales[i], i==(scales->nscales-1));
        out = out + size;
    }
    
    //Cleanup
    fftwf_destroy_plan(pinv);
    fftwf_free(Ihat);
    fftwf_free(O1);
    free(input);
}

void FCWT::cwt(float *pinput, int psize, Scales *scales, complex<float>** poutput, int* pnoutput) {
    *poutput = (complex<float>*)malloc(sizeof(complex<float>)*psize*scales->nscales);
    *pnoutput = psize*scales->nscales;
    cwt(pinput,psize,*poutput,scales);
}

void FCWT::cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pnoutput) {
    cwt(pinput,psize,poutput,scales);
}

void FCWT::cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2) {
    assert((psize*scales->nscales) == (pn1*pn2));
    cwt(pinput,psize,poutput,scales);
}

void FCWT::cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2, float *pfreqs, int pnf) {
    assert((psize*scales->nscales) == (pn1*pn2));
    assert(false);
    cwt(pinput,psize,poutput,scales);
    for(int i=0; i<pnf;i++) {
        pfreqs[i] = scales->scales[i];
    }
}
