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

#include "fcwt.h"

namespace fcwt {
    
    int find2power(int n)
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

    //We used Morlet in this benchmark, but as we state in the article, performance is independent of wavelet type. We plan to include more wavelet types at publication.
    void precalculate_morlet(float *wavelet, float cf, int isize) {
        float tmp;
        int i, nt, newsize;
        
        nt = find2power(isize);
        newsize = 1 << nt;
        
        float newsizef = cf/newsize;
        
        //calculate array
        for(i = 0; i < newsize; i++) {
            tmp = (float)(2.0 * i * newsizef - cf);
            tmp = -(tmp * tmp)/2;
            wavelet[i] = exp(tmp);
        }
    }

    void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float *mother, float scale, int isize)
    {
        float isizef = ((float)isize-1);
        float endpointf = (isizef*2.0)/scale;
        float step = (scale/2.0);
        int endpoint = (int)endpointf;
        int endpoint4 = endpoint>>2;
        
        float wav[8];
        __m256* O8 = (__m256*)output;
        __m256* I8 = (__m256*)input;
        __m256* tmp = (__m256*)&wav;
        __m128 ind4 = _mm_set_ps1(0.0);
        __m128 step4 = _mm_set_ps1(step);
        __m128 offset = _mm_set_ps(3,2,1,0);
        
        int ind1;
        
        for(int q4=0; q4<endpoint4; q4++) {
            float q = (float)q4*4;
            
            __m128 qq = _mm_set_ps1(q);
            ind4 = _mm_mul_ps(step4,_mm_add_ps(qq,offset));
            float find4[4]; _mm_store_ps(find4, ind4);
            
            wav[0] = wav[1] = mother[(int)find4[0]];
            wav[2] = wav[3] = mother[(int)find4[1]];
            wav[4] = wav[5] = mother[(int)find4[2]];
            wav[6] = wav[7] = mother[(int)find4[3]];
            
            O8[q4] = _mm256_mul_ps(I8[q4],tmp[0]);
        }
        
        for(float q=endpoint4*4; q<endpoint; q+=1.0) {
            ind1 = step*q;
            output[(int)q][0] = input[(int)q][0]*mother[ind1];
            output[(int)q][1] = input[(int)q][1]*mother[ind1];
        }
        
        return;
    }

    void create_optimization_schemes(int maxsize, int threads, int flags) {
        
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

    void load_optimization_schemes(bool use_optimalization_schemes, int size, int nthread) {
        if(use_optimalization_schemes) {
            if(size <= 1024) {
                std::cout << "Inputsize is too small (N <= 1024) to use optimization." << std::endl;
                return;
            }
            
            char file_for[50];
            sprintf(file_for, "n%d_t%d.wis", size, nthread);
            
            if(!fftwf_import_wisdom_from_filename(file_for)) {
                std::cout << "WARNING: Optimization scheme '" << file_for << "' was not found, fallback to calculation without optimization." << std::endl;
            }
        }
    }

    void main(float *Rinput,float *Routput,int *pstoctave, int *pendoctave, int *pnbvoice,
                     int *pinputsize, float *pcenterfrequency, int nthreads, bool use_optimalization_schemes)
    {
        int abase;
        float centerfrequency, a1;
        fftwf_complex *Ihat, *O1, *O2;
        
        centerfrequency = *pcenterfrequency;
        const int stoctave = *pstoctave;
        const int endoctave = *pendoctave;
        const int nboctave = endoctave-stoctave+1;
        const int nbvoice = *pnbvoice;
        const int inputsize = *pinputsize;
        
        //Find nearest power of 2
        const int nt = find2power(inputsize);
        const int newsize = 1 << nt;
        
        //Initialize input and mother wavelet memory
        float* input = (float*)calloc(newsize,sizeof(float));
        float* mother = (float*)calloc(newsize,sizeof(float));

        //Initialize intermediate result
        Ihat = fftwf_alloc_complex(newsize);
        O1 = fftwf_alloc_complex(newsize);
        
        //Copy input to new input buffer
        memcpy(input,Rinput,sizeof(float)*inputsize);
        memset(Ihat,0,sizeof(fftwf_complex)*newsize);
        memset(O1,0,sizeof(fftwf_complex)*newsize);
        
        omp_set_num_threads(nthreads);
        
        //Initialize FFTW plans
        fftwf_init_threads();
        fftwf_plan pinv;
        
        fftwf_plan_with_nthreads(omp_get_max_threads());
        fftwf_plan p;
        
        //Load optimization schemes if necessary
        load_optimization_schemes(use_optimalization_schemes,newsize,nthreads);

        //Perform forward FFT on input signal
        p = fftwf_plan_dft_r2c_1d(newsize, input, Ihat, FFTW_ESTIMATE);
        fftwf_execute(p);
        fftwf_destroy_plan(p);
        
        pinv = fftwf_plan_dft_1d(newsize, O1, (fftwf_complex*)Routput, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        //Generate mother wavelet function (in this case a Morlet wavelet, but other types are possible without performance drop)
        precalculate_morlet(mother,centerfrequency,newsize);
        
        float nbvoicef = (float)nbvoice;
        void *pt;
        
        //Loop over all scales and perform scale-dependent operations
        long ind=0;
        for(int i = stoctave; i <= endoctave; i++) {
            abase = 1 << i;
            for(int j=0; j < nbvoice; j++) {
                
                //Calculate scale
                a1 = abase*(pow(2.0,(float)(j)/(nbvoicef)));
                
                //Perform daughter wavelet generation and multiplication with the Fourier transformed input signal
                daughter_wavelet_multiplication(Ihat,O1,mother,a1,newsize);
                
                if((inputsize*nboctave*nbvoice*2 - (ind+(newsize*2))) > 0) {
                    pt = Routput;
                    
                    std::size_t space = 16;
                    std::align(16,sizeof(fftwf_complex),pt,space);
                    
                    fftwf_execute_dft(pinv,O1,(fftwf_complex*)pt);
                    
                    Routput = Routput + inputsize*2;
                    
                }
                ind += inputsize*2;
            }
        }
        
        //Cleanup
        fftwf_destroy_plan(pinv);
        fftwf_free(Ihat);
        fftwf_free(O1);
        fftwf_cleanup_threads();
        free(input);
        free(mother);
    }

    void cwt(float *input, int inputsize, float *output, int pstoctave, int pendoctave, int pnbvoice, float c0, int nthreads, bool use_optimalization_schemes) {
            main(input, output, &pstoctave, &pendoctave, &pnbvoice, &inputsize, &c0, nthreads, use_optimalization_schemes);
        }
}
