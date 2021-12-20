//
//  fcwtplan.cpp
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

#include <stdio.h>
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include "../fcwt.h"

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f

//The gateway function
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("fCWT:nrhs","Three inputs required. (size, num. of threads, optimization type)");
    }

    if( !mxIsDouble(prhs[0]) ||
        mxIsComplex(prhs[0]) ||
        mxGetNumberOfElements(prhs[0]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:size:notScalar",
                        "The maximal size must be a scalar.");
    }
    if( !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1]) ||
        mxGetNumberOfElements(prhs[1]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:nthreads:notScalar",
                        "Number of threads must be a scalar.");
    }
    if ( mxIsChar(prhs[2]) != 1)
        mexErrMsgIdAndTxt( "fCWT:type:notString",
                        "Optimization type must be a string.");
    
    int buflen;
    buflen = (mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
    char *opt = (char*)mxCalloc(buflen,sizeof(char));
    mxGetString(prhs[2],opt,buflen);
    
    int size = mxGetScalar(prhs[0]);
    int nthreads = mxGetScalar(prhs[1]);
    
    int method = 0;
    if(!strcmp(opt,"estimate")) {
        method = FFTW_ESTIMATE;
        mexPrintf("Using 'estimate' to calculate plan ");
    }
    if(!strcmp(opt,"measure")) {
        method = FFTW_MEASURE;
        mexPrintf("Using 'measure' to calculate plan ");
    }
    if(!strcmp(opt,"patient")) {
        method = FFTW_PATIENT;
        mexPrintf("Using 'patient' to calculate plan ");
    }
    if(!strcmp(opt,"exhaustive")) {
        method = FFTW_EXHAUSTIVE;
        mexPrintf("Using 'exhaustive' to calculate plan ");
    }
    
    mexPrintf("using N:%d and %d threads.",size,nthreads);
    fcwt::create_optimization_schemes(size,nthreads,method);
}
