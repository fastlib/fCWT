#ifndef SCALES_H
#define SCALES_H

#include "wavelet.h"

namespace fcwt {
    enum class ScaleType {
        FCWT_LINSCALES,
        FCWT_LOGSCALES,
        FCWT_LINFREQS
    };

    class Scales {
        public:
        FCWT_LIBRARY_API Scales(ScaleType st, int fs, float f0, float f1, int fn);

        void FCWT_LIBRARY_API getScales(const std::vector<float>& pfreqs) noexcept;

        void FCWT_LIBRARY_API getFrequencies(std::vector<float>& pfreqs) const noexcept;

        std::vector<float> scales;

        int fs;

        int nscales;

        private:
        static bool check_nyquist_satisfied(float f, int fs) noexcept;

        void calculate_logscale_array(float base, int fs, float f0, float f1, int fn) noexcept;

        void calculate_linscale_array(int fs, float f0, float f1, int fn) noexcept;

        void calculate_linfreq_array(int fs, float f0, float f1, int fn) noexcept;
    };
} // fcwt

#endif //SCALES_H
