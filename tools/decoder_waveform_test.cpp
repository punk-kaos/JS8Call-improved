#include "decoder_waveform.h"
#include <array>
#include <cstdio>

int main() {
    std::array<int, 79> tones{};
    for (int i = 0; i < 79; ++i) tones[i] = i % 8;
    int failures = 0;
    for (double snr : {-12.0, -24.0, -32.0, -40.0}) {
        auto const a = js8::test::synthesize(tones, 1920, 180000, snr,
                                            1003.7, 0.3, 0.7, 0.02, 6000.5, 42);
        auto const b = js8::test::synthesize(tones, 1920, 180000, snr,
                                            1003.7, 0.3, 0.7, 0.02, 6000.5, 42);
        bool const ok = a.samples == b.samples &&
            std::abs(a.sampleSnrDb - (snr - 10.0 * std::log10(2.0))) < 0.05 &&
            std::abs(a.snr2500Db - a.sampleSnrDb -
                     10.0 * std::log10(6000.0 / 2500.0)) < 1e-10 &&
            std::all_of(a.samples.begin(), a.samples.end(),
                        [](auto v) { return std::abs(static_cast<int>(v)) <= 30000; });
        std::printf("%s setting=%.1f sampleSNR=%.3f SNR2500=%.3f scale=%.3f\n",
                    ok ? "PASS" : "FAIL", snr, a.sampleSnrDb, a.snr2500Db, a.pcmScale);
        failures += !ok;
    }
    // A frame entirely outside the window contributes no invented edge tone.
    auto const absent = js8::test::synthesize(tones, 1920, 1000, -20,
                                             1003.7, 0, 0, 0, 2000, 42);
    failures += absent.signalPower != 0.0;
    std::printf("%s out-of-window signal is absent\n",
                absent.signalPower == 0.0 ? "PASS" : "FAIL");
    return failures ? 1 : 0;
}
