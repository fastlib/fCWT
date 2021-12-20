#include "main.h"

namespace wavelib {
    void cwt(double *input, int nvoi, int noct, int size) {
        //Wavelib Rafat
        char *wave = (char*)"morlet";
        char *type = (char*)"pow";

        double dj = 1.0 / (double)nvoi;
        int J = noct * nvoi;
        int a0 = 2;
        double s0 = 2.0;

        cwt_object wt = cwt_init(wave, 2*PI, size, 1.0, J);
        setCWTScales(wt, s0, dj, type, a0);
        
        cwt(wt,input);
        cwt_free(wt);
    }
}
