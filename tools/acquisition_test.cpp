// Deterministic tests against the production acquisition/extraction code.
#include <QCoreApplication>
#include <QLoggingCategory>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <mutex>
#include "JS8_Include/commons.h"
#include "JS8_Mode/JS8.h"
#include "decoder_waveform.h"

struct dec_data dec_data;
struct specData specData;
std::mutex fftw_mutex;
Q_LOGGING_CATEGORY(decoder_js8, "decoder.js8", QtWarningMsg)
#ifdef JS8_BENCHMARK_VERIFY_COARSE_SYNC
std::atomic<std::uint64_t> benchmarkSyncMismatches{};
#endif
#ifdef JS8_BENCHMARK_VERIFY_RANK_SCAN
std::atomic<std::uint64_t> benchmarkRankMismatches{};
std::atomic<std::uint64_t> benchmarkRankComparedScores{};
#endif
#define JS8_DECODER_TEST_ACCESS
#include "../JS8_Mode/JS8.cpp"

int failures = 0;
int captureBoundaryReports = 0;
int acceptedPilotModels = 0;
void check(bool ok, char const *name) {
    std::printf("%s %s\n", ok ? "PASS" : "FAIL", name);
    failures += !ok;
}

template <typename Mode> void testMode(char const *name) {
    std::printf("[%s]\n", name);
    auto d = std::make_unique<DecodeMode<Mode>>();
    constexpr float fs = 12000.0f / Mode::NDOWN;
    for (int tone = 0; tone < 8; ++tone) {
        std::array<std::complex<float>, Mode::NDOWNSPS> raw{};
        for (int n = 0; n < Mode::NDOWNSPS; ++n)
            raw[n] = std::polar(1.0f, TAU * tone * n / Mode::NDOWNSPS);
        float const energy = d->matchedToneEnergy(raw, tone);
        check(std::abs(energy / (Mode::NDOWNSPS * Mode::NDOWNSPS) - 1.0f) < 1e-5f,
              "timing correlator matches the positive-frequency pilot");
        for (float residual : {-0.2f, 0.2f}) {
            for (int n = 0; n < Mode::NDOWNSPS; ++n)
                raw[n] = std::polar(1.0f, TAU * (tone + residual) * n /
                                                  Mode::NDOWNSPS);
            d->csymb = raw;
            fftwf_execute(d->plans[DecodeMode<Mode>::Plan::CS]);
            auto const measured = d->pilotResidualHz(tone, fs);
            check(measured && ((*measured > 0.0f) == (residual > 0.0f)),
                  "residual estimator sign is correct, including wrapped DC");
            js8::FrequencyTracker tracker;
            tracker.reset(residual * fs / Mode::NDOWNSPS, fs);
            auto replay = raw;
            tracker.apply(raw.data(), Mode::NDOWNSPS);
            js8::aided::replayTrackerCorrection(replay.data(), Mode::NDOWNSPS,
                                               tracker.currentHz(), fs);
            float const corrected = d->matchedToneEnergy(raw, tone);
            check(raw == replay &&
                  corrected / (Mode::NDOWNSPS * Mode::NDOWNSPS) > 0.9999f,
                  "tracker removes residual and aided replay is identical");
        }
    }
    int const start = Mode::NDOWNSPS;
    d->m_receivedSamples = Mode::NMAX;
    d->m_basebandSamples = Mode::NDFFT2;
    d->cd0.fill(ZERO);
    for (int block = 0; block < 3; ++block)
        for (int k = 0; k < 7; ++k)
            for (int n = 0; n < Mode::NDOWNSPS; ++n)
                d->cd0[start + (36 * block + k) * Mode::NDOWNSPS + n] =
                    d->csyncs[block][k][n];
    float const expected = 3.0f * std::pow(7.0f * Mode::NDOWNSPS, 2);
    float const all = d->syncjs8d(start, 0.0f);
    check(std::abs(all - expected) / expected < 1e-5f,
          "positive start uses all three available Costas blocks");
    // Sweep the full coarse-bin uncertainty with real complex pilot waveforms.
    auto const pilotLine = [&](float residual) {
        d->cd0.fill(ZERO);
        for (int block = 0; block < 3; ++block)
            for (int k = 0; k < 7; ++k)
                for (int n = 0; n < Mode::NDOWNSPS; ++n) {
                    int const sample = start + (36 * block + k) * Mode::NDOWNSPS + n;
                    d->cd0[sample] = d->csyncs[block][k][n] *
                        std::polar(1.0f, TAU * residual * sample / fs);
                }
    };
    for (float fraction : {-0.5f, -0.37f, -0.1f, 0.0f, 0.1f, 0.37f, 0.5f}) {
        float const residual = fraction * Mode::DF;
        pilotLine(residual);
        auto const found = d->findFineSync(start + 1);
        check(found && found->start == start &&
              std::abs(found->residualHz - residual) <= d->kFineHz * 0.51f,
              "joint fine sync covers the full coarse FFT bin");
    }
    pilotLine((d->kCoarseRadius + 1) * d->kCoarseHz);
    auto const extended = d->findFineSync(start);
    int const limit = (2 * Mode::NQSYMBOL + 1) *
        (2 * d->kCoarseRadius + 1 + d->kBoundaryExtension) +
        2 * 3 * (2 * d->kFineDivisions + 1) + 2;
    check(extended && extended->extended && !extended->atBoundary &&
          extended->evaluations <= limit &&
          std::abs(extended->residualHz -
              (d->kCoarseRadius + 1) * d->kCoarseHz) < d->kFineHz * 0.51f,
          "frequency-edge optimum gets one bounded extension");
    pilotLine(d->kFineRadius * d->kFineHz);
    auto const boundary = d->findFineSync(start);
    // Outside the nominal half-bin capture, sidelobes can win the first stage.
    // Require bounded work; do not promise acquisition beyond that envelope.
    check(boundary && boundary->evaluations <= limit &&
          std::abs(boundary->residualHz) <= d->kFineRadius * d->kFineHz,
          "out-of-nominal-capture signal cannot force recursive widening");
    if (boundary && boundary->atBoundary) ++captureBoundaryReports;
    bool sameScores = true;
    for (int grid = -d->kFineRadius; grid <= d->kFineRadius; ++grid) {
        float const hz = grid * d->kFineHz;
        sameScores &= d->syncjs8d(start, hz) == d->template syncjs8d<true>(
            start, hz, &d->syncReferences[grid + d->kFineRadius]);
    }
    check(sameScores, "cached references match direct scores over the entire grid");
    pilotLine(0.0f);
    int const outsideInitial = start - Mode::NQSYMBOL - 1;
    auto const constrained = d->findFineSync(outsideInitial);
    check(constrained && constrained->start >= outsideInitial - Mode::NQSYMBOL &&
          constrained->start <= outsideInitial + Mode::NQSYMBOL,
          "local frequency refinement cannot widen the coarse timing envelope");
    d->m_basebandSamples = start + 72 * Mode::NDOWNSPS;
    float const partial = d->syncjs8d(start, 0.0f);
    check(std::abs(partial - expected) / expected < 1e-5f,
          "two-block sync has the same expected score scale");
    d->m_basebandSamples = 36 * Mode::NDOWNSPS;
    check(d->syncjs8d(start, 0.0f) == 0.0f,
          "insufficient pilot coverage rejects acquisition");
    check(!d->findFineSync(start), "fine search rejects missing pilot coverage");
    check(!d->symbolWindowValid(-1) &&
          !d->symbolWindowValid(d->m_basebandSamples) &&
          d->symbolWindowValid(d->m_basebandSamples - Mode::NDOWNSPS),
          "symbol bounds reject missing windows instead of clamping");

    std::array<std::array<float, NN>, NROWS> mags{};
    std::array<std::array<std::complex<float>, NN>, NROWS> bins{};
    std::array<double, NN> starts{}, shifts{};
    std::array<float, NN> hz{};
    std::array<bool, NN> valid{};
    valid.fill(true);
    for (int k = 0; k < NN; ++k) {
        starts[k] = k * Mode::NDOWNSPS;
        for (int tone = 0; tone < NROWS; ++tone)
            mags[tone][k] = 1.0f + 0.1f * tone;
        mags[k % NROWS][k] = 3.0f;
    }
    d->m_enableCoherentData = false;
    valid[7] = false;
    auto const a = d->processToneBins(mags, bins, starts, shifts, hz, valid);
    for (auto &row : mags) row[7] = 1e6f;
    auto const b = d->processToneBins(mags, bins, starts, shifts, hz, valid);
    check(a.llr0 == b.llr0 && a.llr1 == b.llr1 &&
          a.llr0[0] == 0.0f && a.llr0[1] == 0.0f && a.llr0[2] == 0.0f,
          "missing symbols neither bias noise estimates nor supply LLRs");
    auto const savedMags=mags; auto const savedBins=bins;
    auto const savedStarts=starts, savedShifts=shifts;
    auto const savedHz=hz; auto const savedValid=valid;
    check(!d->refinePilots(mags,bins,starts,shifts,hz,valid) &&
          mags==savedMags && bins==savedBins && starts==savedStarts &&
          shifts==savedShifts && hz==savedHz && valid==savedValid,
          "rejected pilot model preserves every first-pass observation");

    dec_data.params.nfa = 800;
    dec_data.params.nfb = 1400;
    dec_data.params.nfqso = 1000;
    dec_data.params.syncStats = false;
    std::fill(std::begin(dec_data.d2), std::end(dec_data.d2), 0);
    check((*d)(dec_data, 0, Mode::NMAX, [](auto const &) {}) == 0,
          "zero-filled input produces no candidates/outputs");
    check((*d)(dec_data, 0, 0, [](auto const &) {}) == 0,
          "empty input produces no outputs");
    // Both unused FFT tail and previous candidate scratch must be overwritten.
    d->cd0.fill({123.0f, -456.0f});
    d->computeBasebandFFT();
    d->js8_downsample(1000.0f);
    check(std::all_of(d->cd0.begin(), d->cd0.end(),
                      [](auto value) { return value == ZERO; }),
          "downsampling cannot retain a stale candidate tail");

    int tones[NN]{};
    JS8::encode(0, JS8::Costas::array(Mode::NCOSTAS), "TESTTEST1234", tones);
    auto const wave = js8::test::synthesize(tones, Mode::NSPS, Mode::NMAX,
        -6, 1000.0, 0.31, 0.0, 0.0, Mode::ASTART * 12000.0, 1234);
    std::copy(wave.samples.begin(), wave.samples.end(), dec_data.d2);
    auto collect = [&](int pos, int count) {
        std::vector<std::string> outputs;
        (*d)(dec_data, pos, count, [&](JS8::Event::Variant const &event) {
            if (auto const *decoded = std::get_if<JS8::Event::Decoded>(&event))
                outputs.push_back(decoded->data);
        });
        return outputs;
    };
    auto const unwrapped = collect(0, Mode::NMAX);
    check(std::find(unwrapped.begin(), unwrapped.end(), "TESTTEST1234") !=
              unwrapped.end(),
          "coherent-disabled/aided-enabled strong frame decodes");
    // Rotate the same ring into a wrapped representation of the reception.
    constexpr int tail = 300;
    std::rotate(std::begin(dec_data.d2), std::begin(dec_data.d2) + tail,
                std::end(dec_data.d2));
    auto const wrapped = collect(JS8_RX_SAMPLE_SIZE - tail, Mode::NMAX);
    check(wrapped == unwrapped, "wrapped and unwrapped reception outputs match");
    // Scratch/input beyond ksz must not affect a partial reception.
    std::copy(wave.samples.begin(), wave.samples.end(), dec_data.d2);
    int const partialSize = Mode::NMAX * 9 / 10;
    auto const partialOutputs = collect(0, partialSize);
    std::fill(std::begin(dec_data.d2) + partialSize, std::end(dec_data.d2),
              std::numeric_limits<std::int16_t>::max());
    check(collect(0, partialSize) == partialOutputs,
          "unreceived input tail cannot change partial-frame outputs");

    typename DecodeMode<Mode>::TimingScores surface{};
    surface.fill(1.0f);
    for (int offset : {-24, -16, -8, 0, 8, 16, 24})
        surface[Mode::JZ + offset] = 2.0f + 0.01f * offset;
    surface[Mode::JZ + 17] = surface[Mode::JZ + 16]; // Equal-score plateau.
    auto const peaks = d->timingPeaks(1000.0f, surface);
    check(peaks.size() == d->kTimingPeaksPerBin &&
          peaks[0].sync > peaks[1].sync && peaks[1].sync > peaks[2].sync,
          "timing peak list is strength-ordered and capped per frequency");
    surface.fill(0.0f);
    surface[Mode::JZ] = 3.0f;
    surface[Mode::JZ + 1] = 3.0f;
    surface[Mode::JZ + 3] = 2.5f; // Separate maximum inside duplicate radius.
    surface[Mode::JZ + 8] = 1.6f; // Weaker independent timing hypothesis.
    auto const separate = d->timingPeaks(1000.0f, surface);
    check(separate.size() == 2 && separate[0].sync == 3.0f &&
          separate[1].sync == 1.6f &&
          separate[0].step == Mode::TSTEP * 0.5f,
          "plateaus/neighbors collapse but a weaker separated peak survives");
    surface.fill(std::numeric_limits<float>::quiet_NaN());
    check(d->timingPeaks(1000.0f, surface).empty(),
          "nonfinite timing surfaces supply no candidates");

    d->sync.clear();
    std::vector<float> maxima{1, 2, 3, 4, 5};
    d->sync.emplace(1000.0f, 0.0f, 3.0f);
    for (int i = 0; i < 10; ++i)
        d->sync.emplace(1001.0f, i * Mode::TSTEP, 0.1f);
    bool const normalized = d->normalizeSync(maxima);
    check(normalized && d->sync.template get<Tag::Freq>().begin()->sync == 1.0f,
          "extra timing peaks do not alter the per-frequency normalization reference");
    std::vector<float> empty;
    std::vector<float> zeros(5, 0.0f);
    check(!d->normalizeSync(empty) && !d->normalizeSync(zeros),
          "empty/zero normalization populations reject safely");

    d->sync.clear();
    float const symbolSeconds = Mode::NSPS / 12000.0f;
    d->sync.emplace(1000.0f, 0.0f, 3.0f);
    d->sync.emplace(1000.0f + Mode::AZ * 0.5f, symbolSeconds * 0.5f, 2.8f);
    d->sync.emplace(1000.0f, symbolSeconds * 2.0f, 2.0f);
    d->sync.emplace(1000.0f + Mode::AZ * 2.0f, 0.0f, 1.8f);
    d->sync.emplace(1000.0f, symbolSeconds * 4.0f, 1.3f);
    d->m_enableDeepSearch = true;
    auto const retained = d->extractSyncCandidates();
    check(retained.normal.size() == 3 && retained.deep.size() == 1,
          "2D suppression preserves time-separated and frequency-separated signals");
    d->sync.clear();
    for (std::size_t i = 0; i < NMAXCAND + 10; ++i)
        d->sync.emplace(1000.0f + i * Mode::AZ * 2.0f, 0.0f, 2.0f);
    check(d->extractSyncCandidates().normal.size() == NMAXCAND,
          "normal candidate budget remains capped");
    d->sync.clear();
    for (std::size_t i = 0; i < NMAXDEEP + 10; ++i)
        d->sync.emplace(1000.0f + i * Mode::AZ * 2.0f, 0.0f, 1.3f);
    check(d->extractSyncCandidates().deep.size() == NMAXDEEP,
          "deep candidate budget remains capped");
    d->sync.clear();
    for (std::size_t i = 0; i < NMAXSECONDARY + 10; ++i)
        d->sync.emplace(1000.0f + i * Mode::AZ * 2.0f, 0.0f, 2.0f, true);
    d->sync.emplace(500.0f, 0.0f, 1.8f);
    auto const limited = d->extractSyncCandidates();
    check(limited.normal.size() == NMAXSECONDARY + 1 &&
          std::count_if(limited.normal.begin(), limited.normal.end(),
                        [](Sync const &s) { return s.secondary; }) == NMAXSECONDARY,
          "additional timing budget cannot starve primary frequency hypotheses");

    for (float fraction : {-0.49f, 0.49f}) {
        double const carrier = (std::round(1000.0 / Mode::DF) + fraction) * Mode::DF;
        auto const offsetWave = js8::test::synthesize(tones, Mode::NSPS,
            Mode::NMAX, -6, carrier, 0.31, 0.0, 0.0, Mode::ASTART * 12000.0, 1234);
        std::copy(offsetWave.samples.begin(), offsetWave.samples.end(), dec_data.d2);
        auto const decoded = collect(0, Mode::NMAX);
        check(std::find(decoded.begin(), decoded.end(), "TESTTEST1234") != decoded.end(),
              "production decoder recovers both sides of a coarse-bin midpoint");
    }

    int otherTones[NN]{};
    JS8::encode(0, JS8::Costas::array(Mode::NCOSTAS), "OTHERMSG1234", otherTones);
    double const carrier = std::round(1000.0 / Mode::DF) * Mode::DF;
    constexpr int separationSymbols = 7;
    auto const first = js8::test::synthesize(tones, Mode::NSPS, Mode::NMAX,
        20, carrier, 0.31, 0.0, 0.0, Mode::ASTART * 12000.0, 1234);
    auto const second = js8::test::synthesize(otherTones, Mode::NSPS, Mode::NMAX,
        20, carrier, 1.1, 0.0, 0.0,
        Mode::ASTART * 12000.0 + separationSymbols * Mode::NSPS, 1235);
    for (int i = 0; i < Mode::NMAX; ++i)
        d->dd[i] = first.samples[i] + second.samples[i];
    d->m_receivedSamples = Mode::NMAX;
    d->m_basebandSamples = std::min(Mode::NDFFT2, Mode::NMAX / Mode::NDOWN);
    auto const overlap = d->syncjs8(800, 1400);
    bool firstFound = false, secondFound = false;
    for (auto const *list : {&overlap.normal, &overlap.deep}) {
        for (auto const &p : *list) {
            if (std::abs(p.freq - carrier) <= Mode::DF) {
                firstFound |= std::abs(p.step) < symbolSeconds * 0.5f;
                secondFound |= std::abs(p.step - separationSymbols * symbolSeconds) <
                               symbolSeconds * 0.5f;
            }
        }
    }
    check(firstFound && secondFound,
          "pre-SIC acquisition retains both overlapping same-frequency frames");
    if (!firstFound || !secondFound)
        for (auto const *list : {&overlap.normal, &overlap.deep})
            for (auto const &p : *list)
                if (std::abs(p.freq - carrier) <= Mode::DF * 2)
                    std::printf("  candidate f=%.3f time=%.6f sync=%.3f\n",
                                p.freq, p.step, p.sync);

    // A real complex FSK frame, not pre-manufactured pilot bins: exercise
    // accepted fractional extraction and its exact aided replay contract.
    d->m_enableFreqTracking=d->m_enableTimingTracking=true;
    d->m_enablePilotSmoothing=d->m_enableFractionalTiming=true;
    d->cd0.fill(ZERO);
    constexpr int nominal=64;
    double const actual=nominal+0.375;
    for (int n=0; n<Mode::NDFFT2; ++n) {
        double const s=(n-actual)/Mode::NDOWNSPS;
        if (s<0 || s>=NN) continue;
        int const symbol=static_cast<int>(s);
        double const t=n/static_cast<double>(fs);
        double const phase=TAU*(0.15*t+0.01*t*t+tones[symbol]*(s-symbol));
        d->cd0[n]=std::polar(1000.0f,static_cast<float>(phase));
    }
    mags={}; bins={}; hz.fill(0); shifts.fill(0); valid.fill(true);
    for (int k=0; k<NN; ++k) {
        starts[k]=nominal+k*Mode::NDOWNSPS;
        js8::extractSymbolWindow(d->cd0.data(),d->m_basebandSamples,starts[k],
                                Mode::NDOWNSPS,d->csymb.data());
        fftwf_execute(d->plans[DecodeMode<Mode>::Plan::CS]);
        for (int tone=0; tone<8; ++tone) {
            bins[tone][k]=d->csymb[tone]/1000.0f;
            mags[tone][k]=std::abs(bins[tone][k]);
        }
    }
    auto const accepted=d->refinePilots(mags,bins,starts,shifts,hz,valid);
    if (accepted) {
        ++acceptedPilotModels;
        check(std::abs(starts[0]-actual)<0.3,"accepted pilot model recovers physical fractional frame start");
        js8::extractSymbolWindow(d->cd0.data(),d->m_basebandSamples,starts[7]+shifts[7],
                                Mode::NDOWNSPS,d->csymb.data());
        js8::aided::replayTrackerCorrection(d->csymb.data(),Mode::NDOWNSPS,hz[7],fs);
        fftwf_execute(d->plans[DecodeMode<Mode>::Plan::CS]);
        bool same=true;
        for (int tone=0; tone<8; ++tone) same &= bins[tone][7]==d->csymb[tone]/1000.0f;
        check(same,"aided replay reproduces accepted fractional tone bins bit-exactly");
    }
}

int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    ::setenv("JS8_SOFT_COMBINING", "0", 1);
    testMode<ModeA>("A"); testMode<ModeB>("B"); testMode<ModeC>("C");
    testMode<ModeE>("E"); testMode<ModeI>("I");
    check(captureBoundaryReports > 0,
          "hard frequency-boundary telemetry is exercised by physical signals");
    check(acceptedPilotModels>0,"production pilot refinement activates on real FSK samples");
#ifdef JS8_BENCHMARK_VERIFY_COARSE_SYNC
    check(benchmarkSyncMismatches.load() == 0,
          "all instrumented cached/direct fine-sync scores are bit-identical");
#endif
#ifdef JS8_BENCHMARK_VERIFY_RANK_SCAN
    check(benchmarkRankMismatches.load() == 0 &&
          benchmarkRankComparedScores.load() > 0,
          "all optimized/original coarse timing scores are bit-identical");
#endif
    return failures ? 1 : 0;
}
