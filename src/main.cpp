//
//  main.cpp
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

#include "main.h"


int main2(int argc, char * argv[]) {
    
    int n = 1000; //signal length
    const int fs = 1000; //sampling frequency
    float twopi = 2.0*3.1415;
    
    //3000 frequencies spread logartihmically between 1 and 32 Hz
    const float f0 = 0.1;
    const float f1 = 20;
    const int fn = 200;

    //Define number of threads for multithreaded use
    const int nthreads = 8;

    //input: n real numbers
    std::vector<float> sig(n);

    //input: n complex numbers
    std::vector<std::complex<float>> sigc(n);
    
    //output: n x scales x 2 (complex numbers consist of two parts)
    std::vector<std::complex<float>> tfm(n * fn * 2);
    
    //initialize with 1 Hz cosine wave
    for(auto& el : sig) {
        el = std::cos(twopi* (static_cast<float>(&el - &sig[0]) / static_cast<float>(fs)));
    }

    //initialize with 1 Hz cosine wave
    for(auto& el : sigc) {
        el = std::complex<float>(cos(twopi*((float)(&el - &sigc[0])/(float)fs)), 0.0f);
    }

    //Start timing
    const auto start = std::chrono::high_resolution_clock::now();

    // Initialize a Morlet wavelet having sigma=1.0;
    fcwt::Morlet morl(1.0f);
    fcwt::Wavelet *wavelet = &morl;

    //Create the continuous wavelet transform object
    //constructor(wavelet, nthreads, optplan)
    //
    //Arguments
    //wavelet   - pointer to wavelet object
    //nthreads  - number of threads to use
    //optplan   - use FFTW optimization plans if true
    fcwt::API fcwt(wavelet, nthreads, true, false);

    //Generate frequencies
    //constructor(wavelet, dist, fs, f0, f1, fn)
    //
    //Arguments
    //dist      - FCWT_LOGSCALES | FCWT_LINSCALES for logarithmic or linear distribution of scales across frequency range
    //fs        - sample frequency
    //f0        - beginning of frequency range
    //f1        - end of frequency range
    //fn        - number of wavelets to generate across frequency range
    fcwt::Scales scs(fcwt::ScaleType::FCWT_LINFREQS, fs, f0, f1, fn);

    //Perform a CWT
    //cwt(input, length, output, scales)
    //
    //Arguments:
    //input     - floating pointer to input array
    //length    - integer signal length
    //output    - floating pointer to output array
    //scales    - pointer to scales object
    fcwt.cwt(sigc.data(), n, tfm.data(), &scs);

    //End timing
    const auto finish = std::chrono::high_resolution_clock::now();

    //Calculate total duration
    const std::chrono::duration<double> elapsed = finish - start;
    
    std::cout << "=== fCWT example ===\n";
    std::cout << "Calculate CWT of a 100k sample sinusodial signal using a [" << f0 << "-" << f1 << "] Hz linear frequency range and " << fn << " wavelets.\n";
    std::cout << "====================\n";
    std::cout << "fCWT finished in " << elapsed.count() << " s\n";

    return 0;
}
