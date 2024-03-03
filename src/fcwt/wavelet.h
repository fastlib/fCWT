#ifndef WAVELET_H
#define WAVELET_H

#include <vector>
#include <complex>
#ifdef _WIN32
  #ifdef FCWT_LIBRARY_DLL_BUILDING
    #define FCWT_LIBRARY_API __declspec(dllexport)
  #else
    #if FCWT_LIBRARY_DLL
      #define FCWT_LIBRARY_API __declspec(dllimport)
    #else /* static or header-only library on Windows */
      #define FCWT_LIBRARY_API
    #endif
  #endif
#else /* Unix */
  #define FCWT_LIBRARY_API
#endif

namespace fcwt {
    class Wavelet {
        public:
        virtual ~Wavelet() = default;

        virtual void generate(std::vector<std::complex<float>> &pwav, int size, float scale) = 0;

        virtual void generate(int size) = 0;

        [[nodiscard]] virtual int getSupport(float scale) const noexcept = 0;

        virtual void getWavelet(float scale, std::vector<std::complex<float>>& pwav, int pn) = 0;

        int width = 0;

        bool imag_frequency = false;

        bool doublesided = false;

        std::vector<float> mother;
    };
};

#endif //WAVELET_H
