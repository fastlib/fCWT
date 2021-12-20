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
#include "../fcwt.h"

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f

//The gateway function
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs != 6) {
        mexErrMsgIdAndTxt("fCWT:nrhs","Six inputs required. (inputArray, c0, fs, s0, s1, Octave resolution)");
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
                        "Starting octave must be a scalar.");
    }
    if( !mxIsDouble(prhs[4]) ||
        mxIsComplex(prhs[4]) ||
        mxGetNumberOfElements(prhs[4]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Ending octave must be a scalar.");
    }
    if( !mxIsDouble(prhs[5]) ||
        mxIsComplex(prhs[5]) ||
        mxGetNumberOfElements(prhs[5]) != 1 ) {
        mexErrMsgIdAndTxt("fCWT:notScalar",
                        "Octave resolution must be a scalar.");
    }
    
    float *inMatrix;       /* 1xN input matrix */
    mwSize c0, fs, s0, s1, res;
    mwSize ncols, nvoi;           /* size of matrix */
    
    inMatrix = mxGetSingles(prhs[0]);
    ncols = mxGetN(prhs[0]);

    c0 = mxGetScalar(prhs[1]);
    fs = mxGetScalar(prhs[2]);
    s0 = mxGetScalar(prhs[3]);
    s1 = mxGetScalar(prhs[4]);
    nvoi = mxGetScalar(prhs[5]);
    
    if(pow(2,s1) > ncols) {
        mexErrMsgIdAndTxt("fCWT:length",
                        "Ending octave measures wavelenghts longer than the window size...");
    }
    
    mwSize noct = s1-s0+1;
    mwSize dims[2] = {ncols,noct*nvoi};
    
    plhs[0] = mxCreateNumericArray(2,dims, mxSINGLE_CLASS, mxCOMPLEX);
    float* outMatrixf = (float*)mxGetComplexSingles(plhs[0]);
    
    mwSize dimsscale[2] = {1,noct*nvoi};
    plhs[1] = mxCreateNumericArray(2,dimsscale,mxSINGLE_CLASS, mxREAL);
    float* outScales = (float*)mxGetSingles(plhs[1]);
    
    float abase;
    int ind = 0;
    for(int i = s0; i <= s1; i++) {
        abase = 1 << i;
        for(int j=0; j < nvoi; j++) {
            outScales[ind] = fs/(abase*(pow(2.0,(float)(j)/((float)nvoi))));
            ind++;
        }
    }

    fcwt::cwt(inMatrix, ncols, outMatrixf, (int)s0, (int)s1, (int)nvoi, (float)c0, 8, true);
}

