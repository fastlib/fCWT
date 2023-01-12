#pragma once

namespace rwave {
    
    void four1(double data[],int nn,int isign);
    void double_fft(double *Or,double *Oi,double *Ir,double *Ii,int isize,int isign);
    void multi(double *Ri1, double *Ii1, double *Ri2, double *Or,double *Oi, int isize);
    void morlet_frequency(double cf,double scale,double *w,int isize);
    void Scwt_morlet(double *Rinput,double *Iinput,double *Oreal,double *Oimage,int *pnboctave,int *pnbvoice,int *pinputsize, int *pdefaultsize, double *pcenterfrequency);
    void cwt(double *input, int inputsize, double* output, int pnboctave,int pnbvoice, int pinputsize);
    
}
