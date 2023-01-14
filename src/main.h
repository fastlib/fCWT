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
#include <string>
#include <fstream>
#include <regex>

#include <iostream>

#include <chrono>
#include <thread>
#include "../libs/fftw3.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
    #include <pthread.h>
#endif
#if defined(__AVX__)
    #include <immintrin.h>
    #define AVX
#endif

using namespace std;

#define PI                    3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f

#include "fcwt/fcwt.h"
#include "rwave-bench.h"
#include "wavelib-bench.h"
