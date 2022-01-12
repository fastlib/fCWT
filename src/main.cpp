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

//Print CLI help text
static void show_usage(std::string name)
{
    std::cerr << "Usage: ./fcwt <method> <number of samples> <number of threads> <OPTIONAL:optimization method>\n"
              << "Method:\n"
              << "\t-h,--help\t Show this help message\n"
              << "\t-fcwt\t\t Run fCWT benchmark for N samples and K threads\n"
              << "\t-fcwtoptimize\t Calculate fCWT optimization plan for N samples, K threads using specified optimization method\n"
              << "\t-rwave\t\t Run Rwave benchmark for N samples\n"
              << "\t-wavelib\t Run Wavelib benchmark for N samples\n"
    << "Optimization methods (when method= -fcwtoptimize):\n"
    << "\t-estimate\t Fast but low performance\n"
    << "\t-measure\t Medium fast and medium performance performance\n"
    << "\t-patient\t Slow but higher performance\n"
    << "\t-exhaustive\t Slowest but (sometimes) highest performance\n"
              << std::endl;
}

//Calculate and print mean and variance of times array
static void show_stats(chrono::duration<double> *times, int runs)
{
    double total = 0.0;
    double mean = 0.0;
    double std = 0.0;
    for(int i=0; i<runs; i++) {
        total += times[i].count();
    }
    mean = total/runs;
    total = 0.0;

    for(int i=0; i<runs; i++) {
        total += (times[i].count() - mean)*(times[i].count() - mean);
    }
    std = sqrt(total/(runs-1));

    cout << " | elapsed avg time: " << mean << "s (sd: " << std << "s) on " << runs << " runs\n";
    
    cout << "[";
    for(int i=0; i<runs; i++) {
        cout << times[i].count() << ",";
    }
    cout << "]\n";
}


int main(int argc, char * argv[]) {
    
    //Initialize variables with default values
    string algorithm = "";
    int size = 1000;
    int nthreads = 8;
    string optimization = "";
    
    //If user requests help
    if ((argv[1] == "-h") || (argv[1] == "--help")) {
        show_usage(argv[0]);
        return 0;
    }
    
    //Return help text if the user inputs the wrong number of arguments
    if(argc<4) {
        show_usage(argv[0]);
        return 0;
    }
    
    regex numberregex("^[0-9]+$");
    regex letterregex("^-[a-zA-Z]+$");
    
    //Check algorithm, number of samples and number of threads
    if(regex_match(argv[1], letterregex) && regex_match(argv[2], numberregex) && regex_match(argv[3], numberregex)) {
        algorithm = argv[1];
        size = stoi(argv[2]);
        nthreads = stoi(argv[3]);
        cout << "Algorithm: " << algorithm << endl;
        cout << "Number of samples: " << size << endl;
        cout << "Number of threads: " << nthreads << endl;
    } else {
        cerr << "ERROR: Please define the method, number of samples and number of threads correctly." <<  '\n';
        return 0;
    }
    
    //If optimisation is requested, one additional argument explaining the method is needed
    if(algorithm=="-fcwtoptimize") {
        if(argc==5) {
            if(regex_match(argv[4], letterregex)) {
                optimization = argv[4];
            } else {
                optimization = "measure";
                cout << "ERROR: No optimization method recognized, use -h to see which methods are available" <<  '\n';
                return 0;
            }
        } else {
            optimization = "measure";
            cout << "ERROR: No optimization method defined, use -h to see which methods are available" <<  '\n';
            return 0;
        }
    }

    //Set CWT parameters and initialize demo signals
    const int noct = 6;
    const int nvoi = 500;
    const int sigoutsize = size*noct*nvoi*2;
    float c0 = 2*PI;
    float hz = 1;
    int runs = 5;
    
    float *sig1 = (float*)malloc(sizeof(float)*size);
    float *sig2 = (float*)malloc(sizeof(float)*size);
    float *sig3 = (float*)malloc(sizeof(float)*size);
    float *sigout = (float*)malloc(sizeof(float)*sigoutsize);
    
    double *sig1d = (double*)malloc(sizeof(double)*size);
    double *sig2d = (double*)malloc(sizeof(double)*size);
    double *sig3d = (double*)malloc(sizeof(double)*size);
    double *sigoutd = (double*)malloc(sizeof(double)*sigoutsize);
    
    for(int i=0; i<size; i++) {
        //Sig1: dynamic sine wave with varying frequency from 1Hz-7Hz, sampling rate of 64Hz.
        sig1[i] = cos((2.0*3.1415*(hz+((float)(7*i)/size)))*((float)i/64.0));
        sig1d[i] = cos((2.0*3.1415*(hz+((double)(7*i)/size)))*((double)i/64.0));
        
        //Sig2: random numbers between 0-10.
        sig2[i] = ((float)(rand() % 1000))/100.0;
        sig2d[i] = ((double)(rand() % 1000))/100.0;
        
        //Sig3: repeating non-smooth function.
        sig3[i] = (i%10==0);
        sig3d[i] = (i%10==0);
    }
    
    auto start = chrono::high_resolution_clock::now();
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed;

    //Initialize timing array
    chrono::duration<double> *times = (chrono::duration<double>*)malloc(sizeof(chrono::duration<double>)*runs);
    

    if(algorithm=="-fcwtoptimize") {
        
        //FCWT optimize
        cout << "=========== OPTIMIZING FCWT ============" << endl;
        
        //Use FFTW_MEASURE for low-quality fast optimization
        //Use FFTW_PATIENT for mid-quality optimization
        //Use FFTW_EXHAUSTIVE for high-quality slow optimization
        
        int opt = 0;
        
        if(optimization=="-estimate") {
            opt = FFTW_ESTIMATE;
            cout << "Using very fast but low-quality optimization: FFTW_ESTIMATE" << endl;
        }
        if(optimization=="-measure") {
            opt = FFTW_MEASURE;
            cout << "Using fast but low-quality optimization: FFTW_MEASURE" << endl;
        }
        if(optimization=="-patient") {
            opt = FFTW_PATIENT;
            cout << "Using slow but high-quality optimization: FFTW_PATIENT" << endl;
        }
        if(optimization=="-exhaustive") {
            opt = FFTW_EXHAUSTIVE;
            cout << "Using very slow but very high-quality optimization: FFTW_EXHAUSTIVE" << endl;
        }
        
        fcwt::create_optimization_schemes(size,nthreads,opt);
        
        cout << "=========== OPTIMIZING END ============" << endl;
    }
    if(algorithm=="-fcwt") {
        
        cout << "=========== BENCHMARKING FCWT ============" << endl;
        //FCWT sig1
        cout << "----- First test -----" << endl;
        cout << "Testing with sig1: dynamic sine wave (1Hz-7Hz)" << endl;
        cout << "Sample sig1: [";
        for(int n=0; n< 10; n++) {
            cout << sig1[n] << ",";
        }
        cout << "...]" << endl;
        
        for(int k=0; k<runs; k++) {
            cout << ".";
            
            start = chrono::high_resolution_clock::now();
            
            fcwt::cwt(sig1, size, sigout, 1, noct, nvoi, c0, nthreads, true);
            
            finish = chrono::high_resolution_clock::now();
            times[k] = finish - start;
            
            this_thread::sleep_for(chrono::microseconds(10000000));

        }
        cout << endl;
        cout << algorithm << " on sig1 with length N: " << size;
        show_stats(times,runs);
        
        
        //FCWT sig2
        cout << "----- Second test -----" << endl;
        cout << "Testing with sig2: random floats between 0-10" << endl;
        cout << "Sample sig2: [";
        for(int n=0; n< 10; n++) {
            cout << sig2[n] << ",";
        }
        cout << "...]" << endl;
        
        for(int k=0; k<runs; k++) {
            cout << ".";
            start = chrono::high_resolution_clock::now();

            fcwt::cwt(sig2, size, sigout, 1, noct, nvoi, c0, nthreads, true);

            finish = chrono::high_resolution_clock::now();
            times[k] = finish - start;
            
            this_thread::sleep_for(chrono::microseconds(10000000));
        }
        cout << endl;
        cout << algorithm << " on sig2 with length N: " << size;
        show_stats(times,runs);
        
        
        //FCWT sig3
        cout << "----- Third test -----" << endl;
        cout << "Testing with sig3: Repeating non-smooth function (x%10==0)" << endl;
        cout << "Sample sig3: [";
        for(int n=0; n< 15; n++) {
            cout << sig3[n] << ",";
        }
        cout << "...]" << endl;
        
        for(int k=0; k<runs; k++) {
            cout << ".";
            start = chrono::high_resolution_clock::now();

            fcwt::cwt(sig3, size, sigout, 1, noct, nvoi, c0, nthreads, true);

            finish = chrono::high_resolution_clock::now();
            times[k] = finish - start;
            
            this_thread::sleep_for(chrono::microseconds(10000000));
        }
        cout << endl;
        cout << algorithm << " on sig3 with length N: " << size;
        show_stats(times,runs);
        
        cout << "=========== BENCHMARKING END ============" << endl;
    }
    if(algorithm=="-rwave") {
        
        //RWAVE
        cout << "=========== BENCHMARKING RWAVE ============" << endl;
        for(int k=0; k<runs; k++) {
            cout << ".";
            start = chrono::high_resolution_clock::now();

            rwave::cwt(sig1d, size, sigoutd, noct, nvoi, size);

            finish = chrono::high_resolution_clock::now();
            times[k] = finish - start;
            
            this_thread::sleep_for(chrono::microseconds(10000000));
        }
        cout << endl;
        cout << algorithm << " on sig1 with length N: " << size;
        show_stats(times,runs);
    }
    if(algorithm=="-wavelib") {
        
        //WAVELIB
        cout << "=========== BENCHMARKING WAVELIB ============" << endl;
        for(int k=0; k<runs; k++) {
            cout << ".";
            start = chrono::high_resolution_clock::now();

            wavelib::cwt(sig1d, nvoi, noct, size);

            finish = chrono::high_resolution_clock::now();
            times[k] = finish - start;
            cout << times[k].count();
            
            this_thread::sleep_for(chrono::microseconds(10000000));
        }
        cout << endl;
        cout << algorithm << " on sig1 with length N: " << size;
        show_stats(times,runs);
    }
    
    delete sig1;
    delete sig2;
    delete sig3;
    delete sig1d;
    delete sig2d;
    delete sig3d;
    delete sigout;
    delete sigoutd;
    delete times;
    
    return 0;
}
