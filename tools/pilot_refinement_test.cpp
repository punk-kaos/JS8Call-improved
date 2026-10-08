#include "JS8_Mode/fractional_window.h"
#include "JS8_Mode/pilot_refinement.h"
#include "JS8_Mode/decoder_aided_redemod.h"
#include <cstdio>
#include <random>

int failures=0;
void check(bool ok, char const *name) {
    std::printf("%s %s\n",ok?"PASS":"FAIL",name); failures+=!ok;
}

void interpolation(int window) {
    std::vector<std::complex<float>> line(256), out(window);
    for (int tone=0; tone<8; ++tone) {
        for (int n=0; n<256; ++n)
            line[n]=std::polar(1.0f,static_cast<float>(2*std::numbers::pi*tone*n/window));
        for (double fraction : {0.0,0.125,0.375,0.75}) {
            double const start=64+fraction;
            bool ok=js8::extractSymbolWindow(line.data(),line.size(),start,window,out.data());
            double error=0;
            for (int n=0; n<window; ++n) {
                auto const expected=std::polar(1.0f,static_cast<float>(2*std::numbers::pi*tone*(start+n)/window));
                error=std::max(error,static_cast<double>(std::abs(out[n]-expected)));
            }
            check(ok && error<0.006,"fractional extraction preserves all eight positive-tone phases");
            if (fraction==0.0)
                check(std::equal(out.begin(),out.end(),line.begin()+64),"integer extraction is bit-exact");
        }
    }
    check(!js8::fractionalWindowValid(0.5,window,256) &&
          js8::fractionalWindowValid(0,window,256) &&
          !js8::fractionalWindowValid(256-window+0.5,window,256) &&
          !js8::fractionalWindowValid(NAN,window,256),"fractional filter halo and invalid starts are guarded");
    line[64]={NAN,0};
    check(!js8::extractSymbolWindow(line.data(),line.size(),64,window,out.data()),"corrupt sample refuses extraction");
}

std::vector<js8::pilot::Observation> observations(int window,double rate,
                                                double delta,double drift,double clock) {
    constexpr int costas[3][7]={{0,6,2,3,5,4,1},{1,5,0,2,3,6,4},{2,5,0,6,4,1,3}};
    double const mid=(64+79*window*0.5)/rate;
    std::vector<js8::pilot::Observation> result;
    for (int b=0; b<3; ++b) for (int c=0; c<7; ++c) {
        int const k=36*b+c,tone=costas[b][c];
        double const base=(64+k*window)/rate;
        double const shift=(k%3-1)*0.375, tracker=(k%4-2)*0.07;
        double const t=base+shift/rate, center=t+(window-1)*0.5/rate;
        double const timing=delta+clock*(t-mid);
        double const phase=0.7+2*std::numbers::pi*0.1*center+
            std::numbers::pi*drift*(center-mid)*(center-mid)+
            2*std::numbers::pi*tone*(shift-timing)/window-
            2*std::numbers::pi*tracker*(window+1)*0.5/rate;
        result.push_back({js8::CoherentPilot{base,std::polar(1.0,phase),tone,k,shift,tracker},0.01,true});
    }
    return result;
}

void fitting(int window,double rate) {
    double const mid=(64+79*window*0.5)/rate;
    for (double drift : {-0.08,0.0,0.08}) {
        auto const obs=observations(window,rate,0.75,drift,0.0);
        auto const model=js8::pilot::fit(obs,window,rate,true);
        if (!model.trusted)
            std::printf("  fit failed N=%d rate=%.3f drift=%.3f rms=%.6f delta=%.6f\n",window,rate,drift,model.carrier.rmsRad,model.timingDeltaSamples);
        check(model.trusted && std::abs(model.carrier.fdot-drift)<0.003 &&
              std::abs(model.frequencyAt(mid)-0.1)<0.02 &&
              std::abs(model.timingAt(mid)-0.75)<0.08,
              "pilot fit recovers fractional timing and both drift signs");
        check(model.fits<=40,"pilot fitting work is bounded");
    }
    for (double ppm : {-300.0,300.0}) {
        double const slope=rate*ppm*1e-6;
        auto const model=js8::pilot::fit(observations(window,rate,-0.75,0.0,slope),window,rate,true);
        if (!model.trusted || std::abs(model.timingRateSamplesPerSecond-slope)>0.2*std::abs(slope))
            std::printf("  clock N=%d expected=%.6f fitted=%.6f rms=%.6f trusted=%d\n",window,slope,model.timingRateSamplesPerSecond,model.carrier.rmsRad,model.trusted);
        check(model.trusted && std::abs(model.timingRateSamplesPerSecond-slope)<0.2*std::abs(slope),
              "block-supported clock slope is recovered in both directions");
    }
    auto bad=observations(window,rate,0.0,0.0,0.0);
    std::mt19937 rng(123);
    std::uniform_real_distribution<double> phase(-std::numbers::pi,std::numbers::pi);
    for (auto &o:bad) o.pilot.value=std::polar(1.0,phase(rng));
    check(!js8::pilot::fit(bad,window,rate,true).trusted,"incoherent pilots cannot steer extraction");
    bad=observations(window,rate,0.0,0.0,0.0);
    for (auto &o:bad) o.noisePower=1.0;
    check(!js8::pilot::fit(bad,window,rate,true).trusted,"noise-level pilot evidence falls back");
    bad.resize(14);
    check(!js8::pilot::fit(bad,window,rate,true).trusted,"missing pilot coverage falls back");
    bad=observations(window,rate,0.0,0.0,0.0);
    for (auto &o:bad) if (o.pilot.symbolIndex/36==1) o.pilot.value=-o.pilot.value;
    check(!js8::pilot::fit(bad,window,rate,true).trusted,
          "a discontinuous block phase cannot masquerade as supported drift");
    bad=observations(window,rate,0.0,0.0,0.0);
    for (auto &o:bad) o.pilot.baseTimeSeconds=NAN;
    check(!js8::pilot::fit(bad,window,rate,true).trusted,"nonfinite pilot metadata is excluded safely");
    bad=observations(window,rate,0.0,0.0,0.0);
    bad.assign(21,bad.front());
    check(!js8::pilot::fit(bad,window,rate,true).trusted,"duplicate pilots cannot create coverage");
}

int main() {
    interpolation(12);interpolation(20);interpolation(32);
    fitting(32,200);fitting(20,200);fitting(12,240);fitting(32,100);fitting(12,375);
    return failures?1:0;
}
