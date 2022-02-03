![](https://github.com/fastlib/fCWT/blob/main/githubart.png)

The fast Continuous Wavelet Transform (fCWT)
====================================

The fast Continuous Wavelet Transform (fCWT) is a highly optimized C++ library for very fast calculation of the CWT. 


---
**fCWT has been featured on the January 2022 cover of NATURE Computational Science**. In this article, fCWT is compared against eight competitor algorithms, tested on noise resistance and validated on synthetic electroencephalography and in vivo extracellular local field potential data.

> Arts, L.P.A., van den Broek, E.L. The fast continuous wavelet transformation (fCWT) for real-time, high-quality, noise-resistant time–frequency analysis. _Nat Comput Sci_ **2**, 47–58 (2022). https://doi.org/10.1038/s43588-021-00183-z

Features
========

- Calculating CWT 34-120x faster than all competitors
- Very high time-frequency resolution (i.e., it does not rely on wavelet estimation)
- Real-time CWT for signals having sample frequencies of up to 200kHz
- Easy MATLAB integration via compiled MEX-files
- Easy extendable to other wavelet types

Example
=======

```cpp
#include <iostream>
#include <vector>
#include <math.h>
#include <fcwt.h>

int main(int argc, const char * argv[]) {
    
    int n = 100000; //signal length
    float fs = 1000; //sampling frequency
    float twopi = 2.0*3.1415;
    
    //3000 scales spread over 5 octaves
    const int noctaves = 5;
    const int nsuboctaves = 600;

    //input: n real numbers
    std::vector<float> sig(n);
    
    //output: n x scales complex numbers
    std::vector<float> tfm(n*noctaves*nsuboctaves*2);
    
    //initialize with 1 Hz cosine wave
    for(auto& el : sig) {
        el = cos(twopi*((float)(&el - &sig[0])/fs));
    }
    
    //Arguments:
    //input     - floating pointer to input array
    //length    - integer signal length
    //output    - floating pointer to output array
    //startoct  - scale range begin (2^startoct)
    //endoct    - scale range end (2^endoct)
    //suboct    - exponential subdivisions of each octave
    //sigma     - parameter to control time-frequency precision
    //nthreads  - number of threads to use
    //optplans  - use FFTW optimization plans if true
    
    fcwt::cwt(&sig[0], n, &tfm[0], 1, noctaves, nsuboctaves, twopi, 8, false);
    
    return 0;
}
```
_Note: Compile with AVX and OpenMP._

Installation
==========================

Dependencies
------------

- [Cmake] >=3.10
- C++17 compiler ([GCC] 10+, [Clang] or [Microsoft Visual C++][Visual_Studio] 15.7+);
- [FFTW] >=3.3  (is included in the zip-file)
- [OpenMP] >=5

fCWT has been tested on Mac OSX Mojave 10.14.5, Big Sur 11.6 and Windows 10. The benchmark has been performed on a MacBook Pro 2019 having a 2,3 GHz Intel Core i9, 16 GB 2400 MHz DDR4. 

Build time settings
-------------------

Settings that may be specified at build time by using [CMake] variables are:
  1. the flag to build a shared library instead of static (default is on);
  2. the flag to build Matlab MEX-files (default is off);
  3. installation directories.

Details:

|CMake variable|Possible values|Default on Unix|Default on Windows|
|:-------------|:--------------|:--------------|:-----------------|
|**The flag to build the shared library**||||
|BUILD_SHARED_LIBS|On \| Off|On|On|
|**The flag to build Matlab MEX-files**||||
|BUILD_MATLAB|On \| Off|Off|Off|
|**Installation directories**||||
|FCWT_MATLAB_DIR|*a path relative to `build`*|"../MATLAB"|"../MATLAB"|
|CMAKE_INSTALL_PREFIX|*an absolute path*|"/usr/local"|"%ProgramFiles (x86)%\fCWT"|
|FCWT_CMAKE_INSTALL_DIR|*a path relative to CMAKE_INSTALL_PREFIX*|"share/fcwt/cmake"|"cmake"|
|FCWT_LIB_INSTALL_DIR|*a path relative to CMAKE_INSTALL_PREFIX*|"lib"|"lib"|
|FCWT_INCLUDE_INSTALL_DIR|*a path relative to CMAKE_INSTALL_PREFIX*|"include"|"include"|

Installation on Unix
---------------------

    $ git clone https://github.com/fastlib/fCWT.git
    $ cd fCWT
    $ mkdir -p build
    $ cd build
    $ cmake ../
    $ make 
    $ sudo make install

Installation on Microsoft Windows
---------------------------------

Run the Developer Command Prompt for Visual Studio and type:

    > git clone https://github.com/fastlib/fCWT.git
    > cd fCWT
    > mkdir -p build
    > cd build
    > cmake -G "Visual Studio 17 2022" ..
    > cmake --build .

A Visual Studio .SLN file has now been created in the build-folder. This project includes several build targets. To build a shared/static library of fCWT, build the 'fCWT' target. To run the example code, set the 'fCWT_example' target as the start-up project and run the code as 'release'.

To install fCWT, run the Elevated Command Prompt (i.e. the command prompt with administrator privileges) in the build-folder and type:

    > cmake -P cmake_install.cmake

To make the installed DLL available for *any* application that depends on it, the symbolic link to the
fcwt.dll should be created:

  - in %SYSTEMROOT%\System32 for the 64-bit DLL on 64-bit host (or for 32-bit DLL on 32-bit host);
  - in %SYSTEMROOT%\SysWOW64 for the 32-bit DLL on 64-bit host.


Usage
=====

Example code
------------

One can run several commands to test fCWT:

`$ ./fCWT_example -h`
Get information about fCWT usage

`$ ./fCWT_example -fcwtoptimize [Length=100000] [Number of threads=8] [Optimization=estimate|measure|patient|exhaustive]`
Calculate all optimization schemes for all powers of two up to the maximal power of two that is calculated based on the [Length] argument, using [Number of threads] threads and using the optimization method specified. If no optimization method is defined it will use the `measure` method. Estimate is the fastest but lowest quality optimization, exhaustive is the slowest but highest quality optimization. Warning: Exhaustive will take some time!

`$ ./fCWT_example -fcwt [Length=100000] [Number of threads=8]`
Calculate the CWT using fCWT for three signals with length=[Length]. The signals are generated by a dynamic sine wave function, random function and non-smooth function generator, respectively. Depending on optimization usage and system hardware, a run-time in the order of seconds is expected for lengths up to 100.000. Furthermore, one will see that signal content is independent of performance. See the source code for more information.

`$ ./fCWT_example -wavelib [Length=100000] [Number of threads=8]`
Calculate the CWT using Wavelib`s C implementation [1] for a signal with length=[Length]. Depending on system hardware, a run-time in the order of several minutes is expected for lengths up to 100.000.

`$ ./fCWT_example -rwave [Length=100000] [Number of threads=8]`
Calculate the CWT using Rwave's C implementation [2] for a signal with length=[Length]. Depending on system hardware, a run-time in the order of several minutes is expected for lengths up to 100.000.

Benchmark reproduction
----------------------

If one wants to reproduce the benchmark results as stated in the article, one has to use a signal length of 100000, 8 threads and calculate optimization plans with the exhaustive method: `$ ./fCWT_example -fcwtoptimize 100000 8 exhaustive` After that:

1. Run fCWT with: `$ ./fCWT_example -fcwt 100000 8`,
2. Run RWave with: `$ ./fCWT_example -rwave 100000 8`,
3. Run Wavelib with: `$ ./fCWT_example -wavelib 100000 8`,
4. Run additional Python, Matlab and Mathematica scripts found in `/src/benchmark`

By default, the source code performs 10 runs for demonstration purposes. To match the number of runs in the paper, adjust the `runs` variable in `main.cpp:128`. It is recommended to close any background processes. 

MATLAB
---------
You can also enjoy fCWT's extremely fast computation in Matlab using CMake's MEX build option. Unfortunately, MEX-files have to be generated as these files are system and Matlab-version dependent. In the MATLAB-folder you can find MEX-files generated for Mac OSX and Matlab r2020b. I will upload additional MEX-files for other systems and version soon.

When you generated your MEX-files, you can run the example live-script titled `example.mlx`. The example includes basic MATLAB-implementation on how to generate fCWT optimization plans and calculate time-frequency matrices using fCWT. 

Note: Expect a decrease in performance when using fCWT via MATLAB. The official benchmark tests MATLAB's CWT implementation via the MATLAB interface and fCWT via the command line. We advice to use fCWT's MATLAB interface solely for plotting purposes. 


License
=======

fCWT is distributed under Apache 2.0 license. For conditions of distribution and use, see file `LICENSE.TXT`


References
==========

[1]     Rafat. Wavelib. https://github.com/rafat/wavelib. Accessed: 2020-04-22.\
[2]     R. Carmona, W.-L. Hwang, and B. Torresani.Practical Time-Frequency Analysis:  Gabor and wavelettransforms, with an implementation in S. Academic Press, 1998. [p1, 7, 10]

[CMake]: https://cmake.org/
[GCC]: https://gcc.gnu.org/
[Clang]: https://clang.llvm.org/
[OpenMP]: https://www.openmp.org/
[FFTW]: https://www.fftw.org/
[Visual_Studio]: https://www.visualstudio.com/
