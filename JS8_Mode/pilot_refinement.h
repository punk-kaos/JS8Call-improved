/** Bounded pilot-only frequency, drift and symbol-clock inference. */
#pragma once
#include "coherent_likelihood.h"
#include <utility>

namespace js8::pilot {
// ~7.5 dB pilot phase SNR. Steering is stricter than the likelihood blend's
// 0.60-rad fallback, but does not require the 0.20-rad full-coherence regime.
inline constexpr double kMaximumPhaseRms = 0.30;
inline constexpr double kMinimumCoherence = 0.95;
inline constexpr double kMaximumDriftHzPerSecond = 0.15;
struct Observation {
    CoherentPilot pilot;
    double noisePower = 0.0; // Mean power in the other seven valid tones.
    bool winner = false;
};

struct Model {
    bool trusted = false;
    CarrierPhaseFit carrier;
    double frequencyRefSeconds = 0.0;
    double timingRefSeconds = 0.0;
    double timingDeltaSamples = 0.0;
    double timingRateSamplesPerSecond = 0.0;
    int fits = 0;
    double frequencyAt(double seconds) const {
        return carrier.deltaF + carrier.fdot * (seconds - frequencyRefSeconds);
    }
    double timingAt(double seconds) const {
        return timingDeltaSamples + timingRateSamplesPerSecond * (seconds - timingRefSeconds);
    }
};

inline Model fit(std::vector<Observation> const &observations, int window,
                 double rate, bool timingEnabled) {
    Model out;
    if (window < 12 || window > 32 || !(rate > 0.0) || !std::isfinite(rate))
        return out;
    double const baud = rate / window;
    double const radius = std::min(4.0, window * 0.125);
    std::array<int, 3> coverage{}, winners{};
    std::array<bool, 79> seen{};
    std::vector<CoherentPilot> pilots;
    std::vector<double> ratios;
    pilots.reserve(21);
    ratios.reserve(21);
    for (auto const &o : observations) {
        int const block = o.pilot.symbolIndex / 36;
        double const power = std::norm(o.pilot.value);
        if (o.pilot.symbolIndex < 0 || o.pilot.symbolIndex > 78 ||
            o.pilot.tone < 0 || o.pilot.tone >= 8 ||
            !std::isfinite(o.pilot.baseTimeSeconds) ||
            !std::isfinite(o.pilot.timingShiftSamples) ||
            !std::isfinite(o.pilot.trackerHz) ||
            block < 0 || block > 2 || o.pilot.symbolIndex % 36 >= 7 ||
            !std::isfinite(power) || !(power > 0.0) ||
            !std::isfinite(o.noisePower) || o.noisePower < 0.0)
            continue;
        double const variance = std::max(o.noisePower, power * 1e-6);
        if (seen[o.pilot.symbolIndex]) continue;
        seen[o.pilot.symbolIndex] = true;
        auto p = o.pilot;
        p.value /= std::sqrt(variance); // Noise-aware reliability, not amplitude fitting.
        pilots.push_back(p);
        ratios.push_back(power / variance);
        ++coverage[block];
        winners[block] += o.winner;
    }
    if (pilots.size() < 15 ||
        *std::min_element(coverage.begin(), coverage.end()) < 4 ||
        *std::min_element(winners.begin(), winners.end()) < 4)
        return out;
    std::sort(ratios.begin(), ratios.end());
    if (ratios[ratios.size() / 2] < 4.0) return out;
    std::sort(pilots.begin(), pilots.end(), [](auto const &a, auto const &b) {
        return a.symbolIndex < b.symbolIndex;
    });
    auto const fitted = [&](std::vector<CoherentPilot> const &p, int minimum, double span) {
        ++out.fits;
        return fitCarrierPhase(p, minimum, span, window, rate, baud * 0.45,
                               kMaximumDriftHzPerSecond);
    };
    auto const shifted = [&](double delta, double slope, double ref) {
        auto p = pilots;
        for (auto &v : p) {
            double const t = effectiveSymbolTimeSeconds(v.baseTimeSeconds, v.timingShiftSamples, rate);
            v.value *= std::exp(std::complex<double>{0.0,
                2.0 * std::numbers::pi * v.tone * (delta + slope * (t - ref)) / window});
        }
        return p;
    };
    double bestRms = std::numeric_limits<double>::infinity();
    double bestDelta = 0.0;
    auto bestPilots = pilots;
    int const steps = timingEnabled ? static_cast<int>(radius * 4.0) : 0;
    // Baseline wins exact ties. No data tones/codeword guesses enter this fit.
    for (int n = 0; n <= 2 * steps; ++n) {
        int const index = n == 0 ? 0 : (n + 1) / 2 * (n % 2 ? -1 : 1);
        double const delta = index * 0.25;
        auto p = shifted(delta, 0.0, 0.0);
        auto const f = fitted(p, 15, 70.0 * window / rate);
        if (f.fitted && std::isfinite(f.rmsRad) && f.rmsRad < bestRms) {
            bestRms = f.rmsRad;
            out.carrier = f;
            bestDelta = delta;
            bestPilots = std::move(p);
        }
    }
    if (!out.carrier.fitted || bestRms > 0.60) return out;
    out.timingRefSeconds = out.carrier.refTimeSeconds;
    out.timingDeltaSamples = bestDelta;
    // Block-wise timing slopes remove an independent phase intercept in each
    // Costas block. A symbol-clock slope is used only when resolved >3 sigma.
    if (timingEnabled) {
        std::array<double, 3> times{}, deltas{}, variances{};
        for (int b = 0; b < 3; ++b) {
            double sw = 0, st = 0, sx = 0, sy = 0;
            struct TimingRow { double tone, phase, weight; };
            std::vector<TimingRow> rows;
            rows.reserve(7);
            for (auto const &p : bestPilots) {
                if (p.symbolIndex / 36 != b) continue;
                double const t = effectiveSymbolTimeSeconds(p.baseTimeSeconds, p.timingShiftSamples, rate);
                auto y = normalizeFrequencyTrackerPhase(p.value, p.trackerHz, rate, window);
                y = normalizeTimingPhase(y, p.tone, p.timingShiftSamples, window);
                double const residual = detail::wrapPhase(std::arg(y) - predictCarrierPhase(out.carrier, t));
                double const w = std::min(std::norm(p.value), 100.0);
                rows.push_back({static_cast<double>(p.tone), residual, w});
                sw += w;
                st += w * t;
                sx += w * p.tone;
                sy += w * residual;
            }
            double const mx = sx/sw, my = sy/sw;
            double xx = 0, xy = 0, rss = 0;
            for (auto const &r : rows) {
                xx += r.weight * (r.tone - mx) * (r.tone - mx);
                xy += r.weight * (r.tone - mx) * (r.phase - my);
            }
            if (!(xx > 0.0)) return out;
            double const slope = xy/xx;
            for (auto const &r : rows)
                rss += r.weight * std::pow(r.phase - my - slope * (r.tone - mx), 2);
            times[b] = st/sw - out.timingRefSeconds;
            deltas[b] = bestDelta - window/(2.0*std::numbers::pi)*slope;
            variances[b] = std::max(1e-6, std::pow(window/(2.0*std::numbers::pi),2)*rss /
                                         (std::max(1.0, static_cast<double>(rows.size())-2.0)*xx));
        }
        double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int b = 0; b < 3; ++b) {
            double const w = 1.0 / variances[b];
            sw += w;
            sx += w * times[b];
            sy += w * deltas[b];
            sxx += w * times[b] * times[b];
            sxy += w * times[b] * deltas[b];
        }
        double const det=sw*sxx-sx*sx;
        if (!(det>0.0)) return out;
        double const slope=(sw*sxy-sx*sy)/det;
        double const intercept=(sy-slope*sx)/sw;
        double residual=0;
        for (int b=0; b<3; ++b) residual+=std::pow(deltas[b]-intercept-slope*times[b],2);
        if (std::sqrt(residual/3.0) <= 0.10 && std::abs(intercept) <= radius) {
            out.timingDeltaSamples=intercept;
            if (std::abs(slope)>3.0*std::sqrt(sw/det) && std::abs(slope)<=rate*0.001)
                out.timingRateSamplesPerSecond=slope;
        }
        bestPilots=shifted(out.timingDeltaSamples,out.timingRateSamplesPerSecond,out.timingRefSeconds);
        out.carrier=fitted(bestPilots,15,70.0*window/rate);
    }
    if (!out.carrier.fitted || out.carrier.rmsRad > kMaximumPhaseRms ||
        out.carrier.coherence < kMinimumCoherence) return out;
    out.frequencyRefSeconds=out.carrier.refTimeSeconds+(window-1.0)/(2.0*rate);
    for (auto const &p : bestPilots) {
        double const t=effectiveSymbolTimeSeconds(p.baseTimeSeconds,p.timingShiftSamples,rate);
        if (std::abs(out.frequencyAt(t))>=baud*0.45 ||
            std::abs(out.timingAt(t))>radius) return out;
    }
    // Independent local slope support in each block. This does not unwrap
    // across the long data gaps: two-block absolute-phase fits can choose
    // different cycle aliases even on good recordings. Local frequency
    // evidence must instead agree with the global trajectory within its
    // empirically estimated uncertainty (three sigma).
    for (int block = 0; block < 3; ++block) {
        std::vector<double> times, phases, weights;
        double sw=0.0, st=0.0, previous=0.0;
        for (auto const &p : bestPilots) {
            if (p.symbolIndex/36 != block) continue;
            double const t=effectiveSymbolTimeSeconds(p.baseTimeSeconds,p.timingShiftSamples,rate);
            auto value=normalizeFrequencyTrackerPhase(p.value,p.trackerHz,rate,window);
            value=normalizeTimingPhase(value,p.tone,p.timingShiftSamples,window);
            double const phase=times.empty() ? std::arg(value) :
                previous+detail::wrapPhase(std::arg(value)-previous);
            previous=phase;
            double const w=std::min(std::norm(p.value),100.0);
            times.push_back(t);
            phases.push_back(phase);
            weights.push_back(w);
            sw += w;
            st += w * t;
        }
        if (times.size()<4 || !(sw>0.0)) return out;
        double const mean=st/sw;
        std::vector<double> ones(times.size(),1.0), centered(times.size());
        double xx=0.0;
        for (std::size_t i=0; i<times.size(); ++i) {
            centered[i]=times[i]-mean;
            phases[i]-=std::numbers::pi*out.carrier.fdot*centered[i]*centered[i];
            xx+=weights[i]*centered[i]*centered[i];
        }
        double intercept=0.0,slope=0.0;
        if (!detail::solveWeightedLine2(ones.data(),centered.data(),phases.data(),
                                       weights.data(),times.size(),intercept,slope)) return out;
        double rss=0.0;
        for (std::size_t i=0; i<times.size(); ++i)
            rss+=weights[i]*std::pow(phases[i]-intercept-slope*centered[i],2);
        double const standardError=std::sqrt(rss/((times.size()-2)*xx))/(2.0*std::numbers::pi);
        double const predicted=out.carrier.deltaF+out.carrier.fdot*(mean-out.carrier.refTimeSeconds);
        if (!std::isfinite(standardError) ||
            std::abs(slope/(2.0*std::numbers::pi)-predicted)>
                std::max(baud*0.003,3.0*standardError)) return out;
    }
    out.trusted = true;
    return out;
}
} // namespace js8::pilot
