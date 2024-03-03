#include "scales.h"

#include <algorithm>
#include <iostream>

fcwt::Scales::Scales(const ScaleType st, const int fs, const float f0, const float f1, const int fn): fs(fs), nscales(fn)  {
    scales.resize(fn);
    switch(st) {
        case ScaleType::FCWT_LINSCALES:
            calculate_logscale_array(2.0f, fs, f0, f1, fn);
            break;
        case ScaleType::FCWT_LOGSCALES:
            calculate_linscale_array(fs, f0, f1, fn);
            break;
        default:
            calculate_linfreq_array(fs, f0, f1, fn);
    }
}

bool fcwt::Scales::check_nyquist_satisfied(const float f, const int fs) noexcept {
    if (f > static_cast<float>(fs) / 2) { [[unlikely]]
        std::cerr << "Max frequency cannot be higher than the Nyquist frequency fs/2\n";
        return false;
    }
    return true;
}

void fcwt::Scales::getScales(const std::vector<float>& pfreqs) noexcept {
    scales = pfreqs;
};

void fcwt::Scales::getFrequencies(std::vector<float>& pfreqs) const noexcept {
    std::ranges::transform(scales, pfreqs.begin(),
                           [&](const float scale) {return static_cast<float>(fs) / scale;});
};

void fcwt::Scales::calculate_logscale_array(const float base, const int fs, const float f0, const float f1,
                                            const int fn) noexcept {

    // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    const float nf0 = f0;
    const float nf1 = f1;
    const float s0 = static_cast<float>(fs) / nf1;
    const float s1 = static_cast<float>(fs) / nf0;

    if (!check_nyquist_satisfied(f1, fs)) {
        // Cannot pass the nyquist frequency
        return;
    }

    const float power[2] = {std::log(s0) / std::log(base), std::log(s1) / std::log(base)};
    const float dpower = power[1] - power[0];

    for(int i = 0; i < fn; i++) {
        const float log_power = power[0] + (dpower / static_cast<float>(fn - 1)) * static_cast<float>(i);
        scales[i] = std::pow(base, log_power);
    }
}

void fcwt::Scales::calculate_linfreq_array(const int fs, const float f0, const float f1, const int fn) noexcept {

    const float nf0 = f0;
    const float nf1 = f1;
    // If a signal has fs=100hz and you want to measure [0.1-50] Hz, you need scales 2 to 1000;

    if (!check_nyquist_satisfied(f1, fs)) {
        // Cannot pass the nyquist frequency
        return;
    }

    const float df = nf1 - nf0;

    for(int i=0; i < fn; i++) {
        scales[fn - i- 1] = static_cast<float>(fs) / (nf0 + df / static_cast<float>(fn) * static_cast<float>(i));
    }
}

void fcwt::Scales::calculate_linscale_array(const int fs, const float f0, const float f1, const int fn) noexcept {
    // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
    const float s0 = static_cast<float>(fs) / f1;
    const float s1 = static_cast<float>(fs) / f0;

    if (!check_nyquist_satisfied(f1, fs)) {
        // Cannot pass the nyquist frequency
        return;
    }

    const float ds = s1 - s0;

    for (int i = 0; i < fn; i++) {
        scales[i] = s0 + ds/ static_cast<float>(fn) * static_cast<float>(i);
    }
}
