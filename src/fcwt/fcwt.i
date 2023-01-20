
%module fcwt
%{
    #define SWIG_FILE_WITH_INIT
    #include "fcwt.h"
%}

%include "numpy.i"
%include "exception.i"
%include "std_string.i"

%init %{
import_array();
%}

%rename(cwt) cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2);
%rename(ccwt) cwt(complex<float>*pcinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2);

%numpy_typemaps(complex<double> , NPY_CDOUBLE, int)
%numpy_typemaps(complex<float> , NPY_CFLOAT , int)

%apply (float* IN_ARRAY1, int DIM1) {(float* pinput, int psize)};
%apply (complex<float>* IN_ARRAY1, int DIM1) {(complex<float>*pcinput, int psize)};
%apply (float* INPLACE_ARRAY1, int DIM1) {(float *pfreqs, int pnf)};
%apply (complex<float>* INPLACE_ARRAY2, int DIM1, int DIM2) {(complex<float>* poutput, int pn1, int pn2)};
%apply (complex<float>* INPLACE_ARRAY1, int DIM1) {(complex<float>* pwav, int pn)};

%typemap(check) float f0, float f1, int fn {
  if ($1 <= 0) {
    SWIG_exception(SWIG_ValueError, "Expected a non-zero, positive value.");
  }
}

enum SCALETYPE {FCWT_LINSCALES,FCWT_LOGSCALES,FCWT_LINFREQS};

class Wavelet {
public:
    Wavelet() {};
    virtual void generate(float* real, float* imag, int size, float scale) { printf("ERROR [generate time complex]: Override this virtual class"); };
    virtual void generate(int size) { printf("ERROR [generate freq]: Override this virtual class"); };
    virtual int getSupport(float scale) { printf("ERROR [getsupport]: Override this virtual class"); return 0; };
    virtual void getWavelet(float scale, complex<float>* pwav, int pn) { printf("ERROR [getsupport]: Override this virtual class"); };
    
    int width;
    float four_wavelen;
    bool imag_frequency, doublesided;
    float *mother;
}; 


class Morlet : public Wavelet {
public:
    Morlet(float bandwidth); //frequency domain
    
    void generate(int size); //frequency domain
    void generate(float* real, float* imag, int size, float scale); //time domain
    int getSupport(float scale) { return (int)(fb*scale*3.0f); };
    void getWavelet(float scale, complex<float>* pwav, int pn);
    float fb;
};

class Scales {
public:
    Scales(Wavelet *pwav, SCALETYPE st, int fs, float f0, float f1, int fn);
    void getScales(float *pfreqs, int pnf);
    void getFrequencies(float *pfreqs, int pnf);

    float* scales;
    float fs;
    int nscales;
};


class FCWT {
public:
    FCWT(Wavelet *pwav, int pthreads, bool puse_optimalization_schemes, bool puse_normalization);
    
    void create_FFT_optimization_plan(int pmaxsize, std::string poptimizationflags);
    void cwt(float *pinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2);
    void cwt(complex<float>*pcinput, int psize, Scales *scales, complex<float>* poutput, int pn1, int pn2);

    Wavelet *wavelet;
};