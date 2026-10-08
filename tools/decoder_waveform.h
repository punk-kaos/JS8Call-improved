// Deterministic, continuous-phase real FSK with AWGN for decoder experiments.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace js8::test {
struct Waveform {
    std::vector<std::int16_t> samples;
    double pcmScale = 0.0;
    double signalPower = 0.0; // During the actual transmission, before PCM scaling.
    double noisePower = 0.0;
    double sampleSnrDb = 0.0;
    double snr2500Db = 0.0;
};

inline Waveform synthesize(std::span<int const> tones, int symbolSamples,
                           int count, double snrSettingDb, double baseHz,
                           double phase, double offsetHz, double driftHzPerSec,
                           double startSamples, unsigned seed,
                           double clockPpm = 0.0) {
    constexpr double fs = 12000.0;
    constexpr double tau = 2.0 * std::numbers::pi;
    if (tones.empty() || symbolSamples <= 0 || count <= 0 ||
        !std::isfinite(snrSettingDb) || !std::isfinite(startSamples) ||
        !std::isfinite(baseHz) || !std::isfinite(phase) ||
        !std::isfinite(offsetHz) || !std::isfinite(driftHzPerSec) ||
        !std::isfinite(clockPpm) || !(1.0 + clockPpm * 1e-6 > 0.0))
        throw std::invalid_argument("invalid waveform configuration");
    double const variance = std::pow(10.0, -snrSettingDb / 10.0);
    double const sigma = std::sqrt(variance);
    if (!(sigma > 0.0) || !std::isfinite(sigma))
        throw std::invalid_argument("invalid noise variance");
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, sigma);
    std::vector<double> raw(static_cast<std::size_t>(count));
    // A complete symbol advances integer tone phase by 2*pi*tone, so only
    // the current fractional symbol is needed for the continuous FSK phase.
    double maxAbs = 1.0, signalPower = 0.0, noisePower = 0.0;
    int observed = 0;
    for (int i = 0; i < count; ++i) {
        double const s = (i - startSamples) / (1.0 + clockPpm * 1e-6);
        double signal = 0.0;
        if (s >= 0.0 && s < tones.size() * static_cast<double>(symbolSamples)) {
            auto const symbol = static_cast<std::size_t>(s / symbolSamples);
            double const t = s / fs;
            double const fraction = s / symbolSamples - symbol;
            signal = std::cos(phase + tau * ((baseHz + offsetHz) * t +
                                0.5 * driftHzPerSec * t * t +
                                tones[symbol] * fraction));
            signalPower += signal * signal;
            ++observed;
        }
        double const n = noise(rng);
        noisePower += n * n;
        raw[i] = signal + n;
        maxAbs = std::max(maxAbs, std::abs(raw[i]));
    }
    Waveform result;
    // Preserve SNR, leave conversion headroom, and never clip Gaussian noise.
    result.pcmScale = std::min(2000.0, 30000.0 / maxAbs);
    result.samples.reserve(raw.size());
    for (double value : raw)
        result.samples.push_back(static_cast<std::int16_t>(
            std::lround(value * result.pcmScale)));
    result.signalPower = observed ? signalPower / observed : 0.0;
    result.noisePower = noisePower / count;
    result.sampleSnrDb = result.signalPower > 0.0
        ? 10.0 * std::log10(result.signalPower / result.noisePower)
        : -std::numeric_limits<double>::infinity();
    // White real noise occupies 0..fs/2. Report the conventional 2500-Hz
    // reference separately from both the setting and full-band sample SNR.
    result.snr2500Db = result.sampleSnrDb + 10.0 * std::log10((fs / 2.0) / 2500.0);
    return result;
}
} // namespace js8::test
