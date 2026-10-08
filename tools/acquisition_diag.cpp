// Production decoder acquisition/tracking replay. See tools/README.txt.
#include <QCoreApplication>
#include <QLoggingCategory>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <mutex>
#include "JS8_Include/commons.h"
#include "JS8_Mode/JS8.h"
#include "decoder_waveform.h"

struct dec_data dec_data;
struct specData specData;
std::mutex fftw_mutex;
Q_LOGGING_CATEGORY(decoder_js8, "decoder.js8", QtWarningMsg)
#include "../JS8_Mode/JS8.cpp"

template <typename Mode> void sweep(char const *name, int trials, double snr,
                                    unsigned firstSeed, bool impairments) {
    auto decoder = std::make_unique<DecodeMode<Mode>>();
    constexpr char message[] = "TESTTEST1234";
    int tones[NN]{};
    JS8::encode(0, JS8::Costas::array(Mode::NCOSTAS), message, tones);
    struct Scenario { char const *name; double bin, timing, drift, coverage;
                      double clockPpm = 0.0; };
    std::vector<Scenario> scenarios{
             Scenario{"center", 0.0, 0.0, 0.0, 1.0},
             Scenario{"bin_plus", 0.49, 0.0, 0.0, 1.0},
             Scenario{"bin_minus", -0.49, 0.0, 0.0, 1.0},
             Scenario{"late", 0.0, 0.2, 0.0, 1.0},
             Scenario{"early", 0.0, -0.2, 0.0, 1.0},
             Scenario{"drift", 0.2, 0.0, 0.08, 1.0},
             Scenario{"partial", 0.0, 0.0, 0.0, 0.92}};
    if (impairments) {
        scenarios.push_back({"drift_minus", 0.2, 0.0, -0.08, 1.0});
        scenarios.push_back({"fractional", 0.0, 0.375 / Mode::NDOWNSPS, 0.0, 1.0});
        scenarios.push_back({"clock_plus", 0.0, 0.0, 0.0, 1.0, 300.0});
        scenarios.push_back({"clock_minus", 0.0, 0.0, 0.0, 1.0, -300.0});
    }
    for (auto const &scenario : scenarios) {
        int exact = 0, falseOutputs = 0, detected = 0;
        double frequencyError = 0.0, timingError = 0.0;
        double sampleSnr = 0.0, snr2500 = 0.0;
        long long micros = 0;
        double const base = (std::round(1000.0 / Mode::DF) + scenario.bin) * Mode::DF;
        double const start = Mode::ASTART * 12000.0 + scenario.timing * Mode::NSPS;
        for (int trial = 0; trial < trials; ++trial) {
            auto const w = js8::test::synthesize(tones, Mode::NSPS, Mode::NMAX,
                snr, base, 0.31, 0.0, scenario.drift, start, firstSeed + trial,
                scenario.clockPpm);
            sampleSnr += w.sampleSnrDb;
            snr2500 += w.snr2500Db;
            std::copy(w.samples.begin(), w.samples.end(), dec_data.d2);
            dec_data.params.nfa = 800;
            dec_data.params.nfb = 1400;
            dec_data.params.nfqso = static_cast<int>(base);
            dec_data.params.syncStats = true;
            dec_data.params.nutc = code_time(0, 0, 0);
            bool recovered = false, candidate = false;
            auto const *traceCase = std::getenv("JS8_ACQUISITION_TRACE");
            bool const trace = traceCase &&
                std::string(traceCase) == std::string(name) + ":" + scenario.name;
            if (trace) {
                QLoggingCategory::setFilterRules("decoder.js8.debug=true");
                std::fprintf(stderr, "TRACE %s:%s seed=%u begin\n", name,
                             scenario.name, firstSeed + trial);
            }
            auto const before = std::chrono::steady_clock::now();
            (*decoder)(dec_data, 0, static_cast<int>(Mode::NMAX * scenario.coverage),
                [&](JS8::Event::Variant const &event) {
                    if (auto const *s = std::get_if<JS8::Event::SyncState>(&event))
                        candidate |= s->type == JS8::Event::SyncState::Type::CANDIDATE &&
                                     std::abs(s->frequency - base) < Mode::DF;
                    if (auto const *d = std::get_if<JS8::Event::Decoded>(&event)) {
                        if (d->data == message) {
                            if (!recovered) {
                                // Synthesis drift is measured from frame start,
                                // and clock error stretches both phase and time.
                                double const midpoint = 0.5 * NN * Mode::NSPS / 12000.0;
                                frequencyError += std::abs(d->frequency -
                                    (base + scenario.drift * midpoint) /
                                    (1.0 + scenario.clockPpm * 1e-6));
                                timingError += std::abs(d->xdt -
                                    (start / 12000.0 - Mode::ASTART));
                            }
                            recovered = true;
                        } else ++falseOutputs;
                    }
                });
            micros += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - before).count();
            exact += recovered;
            detected += candidate;
            if (trace) {
                std::fprintf(stderr, "TRACE %s:%s seed=%u exact=%d end\n", name,
                             scenario.name, firstSeed + trial, recovered);
                QLoggingCategory::setFilterRules("decoder.js8.debug=false");
            }
        }
        std::printf("%s,%s,%.1f,%d,%d,%d,%d,%.6f,%.6f,%.3f,%.3f,%.3f\n", name,
            scenario.name, snr, trials, detected, exact, falseOutputs,
            exact ? frequencyError / exact : NAN,
            exact ? timingError / exact : NAN, micros / (1000.0 * trials),
            sampleSnr / trials, snr2500 / trials);
        std::fflush(stdout);
    }
}

int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    int const trials = argc > 1 ? std::atoi(argv[1]) : 2;
    double const snr = argc > 2 ? std::atof(argv[2]) : -24.0;
    unsigned const firstSeed = argc > 3 ? static_cast<unsigned>(std::strtoul(argv[3], nullptr, 10)) : 1234;
    bool const impairments = argc > 4 && std::string(argv[4]) == "impairments";
    if (trials <= 0 || !std::isfinite(snr)) return 2;
    // Replays must not accumulate the same deterministic reception across trials.
    ::setenv("JS8_SOFT_COMBINING", "0", 1);
    std::puts("mode,scenario,settingDb,trials,candidate,exact,falseOutputs,meanFreqErrorHz,meanTimingErrorSec,meanMs,sampleSnrDb,snr2500Db");
    sweep<ModeA>("A", trials, snr, firstSeed, impairments);
    sweep<ModeB>("B", trials, snr, firstSeed, impairments);
    sweep<ModeC>("C", trials, snr, firstSeed, impairments);
    sweep<ModeE>("E", trials, snr, firstSeed, impairments);
    sweep<ModeI>("I", trials, snr, firstSeed, impairments);
}
