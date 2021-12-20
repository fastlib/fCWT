
/***************************************************************
 *              (c) Copyright  1997                             *
 *                         by                                   *
 *     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
 *                 Princeton University                         *
 *                 All right reserved                           *
 ***************************************************************/

#include "main.h"

namespace rwave {
    
    int find2power(int n)
    {
        long m, m2;
            
        m = 0;
        m2 = 1<<m; /* 2 to the power of m */
        while (m2-n < 0) {
            m++;
            m2 <<= 1; /* m2 = m2*2 */
        }
        return(m);
    }

    #define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
        
    void four1(double data[], int nn, int isign)
    {
        int n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        double tempr,tempi;
        
        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
            if (j > i) {
                SWAP(data[j],data[i]);
                SWAP(data[j+1],data[i+1]);
            }
            m=n >> 1;
            while (m >= 2 && j > m) {
                j -= m;
                m >>= 1;
            }
            j += m;
        }
        mmax=2;
        while (n > mmax) {
            istep=2*mmax;
            theta=6.28318530717959/(isign*mmax);
            wtemp=sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=sin(theta);
            wr=1.0;
            wi=0.0;
            for (m=1;m<mmax;m+=2) {
                for (i=m;i<=n;i+=istep) {
                    j=i+mmax;
                    tempr=wr*data[j]-wi*data[j+1];
                    tempi=wr*data[j+1]+wi*data[j];
                    data[j]=data[i]-tempr;
                    data[j+1]=data[i+1]-tempi;
                    data[i] += tempr;
                    data[i+1] += tempi;
                }
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            mmax=istep;
        }
    }
        
    void double_fft(double *Or,double *Oi,double *Ir,double *Ii,
                    int isize,int isign)
    {
        double *tmp;
        int nt, newsize, i;
        
        nt = find2power(isize);
        newsize = 1 << nt;
        
        tmp = (double *)malloc(newsize * 2 * sizeof(double));
        memset(tmp,0,newsize * 2 * sizeof(double));
        
        for(i = 0; i < isize; i++) {
            tmp[2 * i] = Ir[i];
            tmp[2 * i + 1] = Ii[i];
        }
        four1(tmp-1,newsize,isign);
        
        
        for(i = 0; i < isize; i++) {
            if(isign == -1) {
                Or[i] = tmp[2 * i]/newsize;
                Oi[i] = tmp[2 * i + 1]/newsize;
            }
            else {
                Or[i] = tmp[2 * i];
                Oi[i] = tmp[2 * i + 1];
            }
        }
        delete tmp;
    }

    void multi(double *Ri1, double *Ii1, double *Ri2, double *Or,
               double *Oi, int isize)
    {
        int i;
        
        for(i = 0; i < isize; i++) {
            Or[i] = Ri1[i] * Ri2[i];
            Oi[i] = Ii1[i] * Ri2[i];
        }
        return;
    }

    void morlet_frequency(double cf,double scale,double *w,int isize)
    {
        double tmp;
        int i, nt, newsize;
        double twopi;
        
        twopi = 6.28318530717959;
        
        nt = find2power(isize);
        newsize = 1 << nt;
        
        for(i = 0; i < isize; i++) {
            tmp = (double)(scale * i * twopi/newsize - cf);
            tmp = -(tmp * tmp)/2;
            w[i] = exp(tmp);
        }
        return;
    }

    void Scwt_morlet(double *Rinput,double *Iinput, double *Oreal, double *Oimage, int *pnboctave, int *pnbvoice,
                     int *pinputsize, int *pdefaultsize, double *pcenterfrequency)
    {
        int nboctave, nbvoice, i, j, inputsize, defaultsize;
        double centerfrequency, a;
        double *Ri2, *Ri, *Ii, *Ri1, *Ii1;
        
        centerfrequency = *pcenterfrequency;
        nboctave = *pnboctave;
        nbvoice = *pnbvoice;
        inputsize = *pinputsize;
        defaultsize = *pdefaultsize;
        
        Ri2 = (double *)malloc(sizeof(double)*inputsize);
        Ri1 = (double *)malloc(sizeof(double)*inputsize);
        Ii1 = (double *)malloc(sizeof(double)*inputsize);
        Ri = (double *)malloc(sizeof(double)*inputsize);
        Ii = (double *)malloc(sizeof(double)*inputsize);
        
        for(i = 0; i < inputsize; i++) {
            Ri[i] = (double)Rinput[i];
            Ii[i] = (double)Iinput[i];
        }
        
        double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);
        
        for(i = 1; i <= nboctave; i++) {
            for(j=0; j < nbvoice; j++) {
                a = (double)(pow(2.0,(double)(i+j/(double)nbvoice)));
                morlet_frequency(centerfrequency,a,Ri2,inputsize);
                multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
                
                double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1);

                Oreal = Oreal + inputsize;
                Oimage = Oimage + inputsize;
            }
        }
        
        delete Ri;
        delete Ii;
        delete Ri2;
        delete Ri1;
        delete Ii1;
    }

    void cwt(double *input, int defaultwindow, double *output, int pnboctave, int pnbvoice, int pinputsize) {
        double w0 = 2*PI;
        int pp = pnbvoice * pnboctave;
        int size = pinputsize;
        
        double *Iinput = (double*)malloc(sizeof(double)*size);
        memset(Iinput,0,sizeof(double)*size);
        
        Scwt_morlet(input, Iinput, &output[0], &output[size*pp], &pnboctave, &pnbvoice, &pinputsize, &defaultwindow, &w0);
        
        delete Iinput;
    }
}
