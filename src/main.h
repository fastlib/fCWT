//
//  main.h
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
#include <unistd.h>
#include <string>
#include <fstream>
#include <regex>

#include <iostream>

#include <thread>
#include "libwavelib/wavelib.h"
#include <fftw3.h>
#ifdef _WIN32
    #include <windows.h>
#else
    #include <pthread.h>
#endif
#include <complex.h>
#include <immintrin.h>
//
//#include <TAxis.h>
//#include <TGraph.h>
//#include <TGraph2D.h>
//#include <TH2F.h>
//#include <TMultiGraph.h>
//#include <TCanvas.h>
//#include <TApplication.h>
//#include <TStyle.h>
//#include <TPad.h>
//#include <TROOT.h>
//#include <TColor.h>
//#include <TFrame.h>
//#include <TVirtualPad.h>

using namespace std;

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f

//GLOBAL FUNCTIONS
inline string GetCurrentWorkingDir( void ) {
    char buff[FILENAME_MAX];
    getcwd( buff, FILENAME_MAX );
    std::string current_working_dir(buff);
    return current_working_dir;
}

#include "fcwt.h"
#include "rwave.h"
#include "wavelib.h"
