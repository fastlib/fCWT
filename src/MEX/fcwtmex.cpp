//
//  fcwtmex.cpp
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

#include "mex.h"
#include "matrix.h"
#include "../fcwt/fcwt.h"

//The gateway function
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs != 7) {
        mexErrMsgIdAndTxt("fCWT:nrhs","Six inputs required. (inputArray, c0, fs, f0, f1, fn, nthreads)");
    }
    if(nlhs != 2) {
        mexErrMsgIdAndTxt("fCWT:nlhs","Two outputs are required. (time-frequency matrix, frequencies)");
    }

    if( !mxIsSingle(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("fCWT:notSingle","Input matrix must be type float (single).");
    }
    if(mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("fCWT:notRowVector",
                          "Input must be a row vector.");
    }

    if( !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1]) ||
        mxGetNumberOfElements(prhs[1]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Wavelet coefficient must be scalar.");
    }
    
    if( !mxIsDouble(prhs[2]) ||
        mxIsComplex(prhs[2]) ||
        mxGetNumberOfElements(prhs[2]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Sampling frequency must be a scalar.");
    }
    if( !mxIsDouble(prhs[3]) ||
        mxIsComplex(prhs[3]) ||
        mxGetNumberOfElements(prhs[3]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Starting frequency must be a scalar.");
    }
    if( !mxIsDouble(prhs[4]) ||
        mxIsComplex(prhs[4]) ||
        mxGetNumberOfElements(prhs[4]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Ending frequency must be a scalar.");
    }
    if( !mxIsDouble(prhs[5]) ||
        mxIsComplex(prhs[5]) ||
        mxGetNumberOfElements(prhs[5]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Number of frequencies must be a scalar.");
    }
    if( !mxIsDouble(prhs[6]) ||
        mxIsComplex(prhs[6]) ||
        mxGetNumberOfElements(prhs[6]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Number of threads must be a scalar.");
    }
    
    float *inMatrix;       /* 1xN input matrix */
    float f0, f1, c0;
    mwSize fs, fn, nthreads;
    mwSize ncols;           /* size of matrix */
    
    inMatrix = mxGetSingles(prhs[0]);
    ncols = mxGetN(prhs[0]);

    c0 = mxGetScalar(prhs[1]);
    fs = mxGetScalar(prhs[2]);
    f0 = mxGetScalar(prhs[3]);
    f1 = mxGetScalar(prhs[4]);
    fn = mxGetScalar(prhs[5]);
    nthreads = mxGetScalar(prhs[6]);
    
    mwSize dims[2] = {ncols,fn};
    
    plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxCOMPLEX);
    complex<float>* outMatrixf = (complex<float>*)mxGetComplexSingles(plhs[0]);
    
    mwSize dimsscale[2] = {1,fn};
    plhs[1] = mxCreateNumericArray(2,dimsscale,mxSINGLE_CLASS, mxREAL);
    float* outScales = (float*)mxGetSingles(plhs[1]);
    
    Wavelet *wavelet;
    Morlet morl(c0);
    wavelet = &morl;

    FCWT fcwt(wavelet, 1, true, true); //threads = 1 because of OMP bug
    Scales scs(wavelet, FCWT_LOGSCALES, fs, f0, f1, fn);

    scs.getScales(outScales, fn);

    for (int i = 0; i < fn; i++) {
        outScales[i] = ((float)fs)/outScales[i];
    }

    mexWarnMsgIdAndTxt("fcwt:nothreads","Threads are currently not supported in Matlab. Using nthreads=1. See Issue #17 on Github.");

    fcwt.cwt(inMatrix, ncols, outMatrixf, &scs);
}

