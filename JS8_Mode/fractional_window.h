/** Shared, validity-aware extraction of complex JS8 symbol windows. */
#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <numbers>

namespace js8 {
inline bool fractionalWindowValid(double start, int count, int available) {
    if (!std::isfinite(start) || count <= 0 || available < count ||
        start < 0.0 || start > available - count)
        return false;
    double const first = std::floor(start);
    if (start == first) return true; // Integer extraction needs no filter halo.
    return first >= 7.0 && first + count - 1 + 8 < available;
}

inline bool extractSymbolWindow(std::complex<float> const *samples, int available,
                                double start, int count,
                                std::complex<float> *out) {
    if (!samples || !out || !fractionalWindowValid(start, count, available))
        return false;
    int const first = static_cast<int>(std::floor(start));
    double const fraction = start - first;
    if (fraction == 0.0) {
        for (int n = 0; n < count; ++n) {
            auto const v = samples[first + n];
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) return false;
            out[n] = v;
        }
        return true;
    }
    // The 12-sample modes alias tones 6/7 above Nyquist. Interpolate around
    // the middle of the actual eight-tone passband (3.5 tones), not around DC,
    // so fractional displacement retains their POSITIVE-tone phase convention.
    constexpr double pi = std::numbers::pi;
    std::array<std::complex<float>, 16> weights{};
    double normalization = 0.0;
    for (int tap = 0; tap < 16; ++tap) {
        double const delta = fraction - (tap - 7);
        double const sinc = std::abs(delta) < 1e-12 ? 1.0 :
                            std::sin(pi * delta) / (pi * delta);
        double const lanczos = sinc * std::sin(pi * delta / 8.0) / (pi * delta / 8.0);
        normalization += lanczos;
        weights[tap] = static_cast<float>(lanczos) * std::polar(
            1.0f, static_cast<float>(2.0 * pi * 3.5 * delta / count));
    }
    if (!(normalization > 0.0) || !std::isfinite(normalization)) return false;
    for (auto &w : weights) w /= static_cast<float>(normalization);
    for (int n = 0; n < count; ++n) {
        std::complex<float> value{};
        for (int tap = 0; tap < 16; ++tap) {
            auto const sample = samples[first + n + tap - 7];
            if (!std::isfinite(sample.real()) || !std::isfinite(sample.imag())) return false;
            value += sample * weights[tap];
        }
        out[n] = value;
    }
    return true;
}
} // namespace js8
