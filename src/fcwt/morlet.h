//
// Created by jonas on 29.02.24.
//

#ifndef MORLET_H
#define MORLET_H

#include <complex>
#include <vector>

#include "wavelet.h"
#include "api.h"

namespace fcwt {
    class Morlet: public Wavelet {
        public:
        FCWT_LIBRARY_API explicit Morlet(float bandwidth); // frequency domain

        ~Morlet() override = default;

        // Frequency domain
        void generate(int size) noexcept override;

        // Time domain
        void generate(std::vector<std::complex<float>> &pwav, int size, float scale) noexcept override;

        [[nodiscard]] int getSupport(float scale) const noexcept override;

        void getWavelet(float scale, std::vector<std::complex<float>>& pwav, int pn) noexcept override;

        float fb;

        private:
        float ifb, fb2;
    };
}

#endif //MORLET_H
