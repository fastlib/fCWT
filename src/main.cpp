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

using namespace std;

int main(int argc, char * argv[]) {
    
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
    std::vector<complex<float>> sigc(n);
    
    //output: n x scales x 2 (complex numbers consist of two parts)
    std::vector<complex<float>> tfm(n*fn);
    
    //initialize with 1 Hz cosine wave
    for(auto& el : sig) {
        el = cos(twopi*((float)(&el - &sig[0])/(float)fs));
    }

    //initialize with 1 Hz cosine wave
    for(auto& el : sigc) {
        el = complex<float>(cos(twopi*((float)(&el - &sigc[0])/(float)fs)), 0.0f);
    }
    
    //Start timing
    auto start = chrono::high_resolution_clock::now();
    
    //Create a wavelet object
    Wavelet *wavelet;
    
    //Initialize a Morlet wavelet having sigma=1.0;
    Morlet morl(1.0f);
    wavelet = &morl;

    //Other wavelets are also possible
    //DOG dog(int order); 
    //Paul paul(int order);

    //Create the continuous wavelet transform object
    //constructor(wavelet, nthreads, optplan)
    //
    //Arguments
    //wavelet   - pointer to wavelet object
    //nthreads  - number of threads to use
    //optplan   - use FFTW optimization plans if true
    FCWT fcwt(wavelet, nthreads, true, false);

    //Generate frequencies
    //constructor(wavelet, dist, fs, f0, f1, fn)
    //
    //Arguments
    //wavelet   - pointer to wavelet object
    //dist      - FCWT_LOGSCALES | FCWT_LINSCALES for logarithmic or linear distribution of scales across frequency range
    //fs        - sample frequency
    //f0        - beginning of frequency range
    //f1        - end of frequency range
    //fn        - number of wavelets to generate across frequency range
    Scales scs(wavelet, FCWT_LINFREQS, fs, f0, f1, fn);

    //Perform a CWT
    //cwt(input, length, output, scales)
    //
    //Arguments:
    //input     - floating pointer to input array
    //length    - integer signal length
    //output    - floating pointer to output array
    //scales    - pointer to scales object
    fcwt.cwt(&sigc[0], n, &tfm[0], &scs);
        
    //End timing
    auto finish = chrono::high_resolution_clock::now();

    //Calculate total duration
    chrono::duration<double> elapsed = finish - start;
    
    cout << "=== fCWT example ===" << endl;
    cout << "Calculate CWT of a 100k sample sinusodial signal using a [" << f0 << "-" << f1 << "] Hz linear frequency range and " << fn << " wavelets." << endl;
    cout << "====================" << endl;
    cout << "fCWT finished in " << elapsed.count() << "s" << endl;

    return 0;
}
