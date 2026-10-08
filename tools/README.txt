In the tools directory there are the following shell scripts for linux:

- Debian_build-deb.sh is an experimental script that builds JS8Call and produces a Debian-style .deb installation file

- Linux-TestBuild-Qt6.sh is primarily for building Qt6 with the required modules for JS8Call from source code. By setting the QT_VERSION variable in the script you can build whatever version of Qt you want to test or experiment with

- Linux-User-Build-JS8Call.sh is an end-user script that will build and install JS8Call for Debian, Redhat/Fedora and Arch Linux systems for either x86_64 or arm64 architectures. This script should also work for Raspberry Pi 5

- tcp.py is a python script that connects to a server running on your own computer (127.0.0.1) at port 2442 to diagnose the JS8Call API for automated digital logging. Once connected, it automatically asks for the station's status, listens continuously for incoming messages from the server, prints those messages out, and ignores minor status updates it doesn't care about

- udp.py is the same thing as above, except using the UDP protocol instead of TCP

- tracking_diag.cpp is a source file that builds a diagnostic harness for JS8 frequency/timing tracking

- whitening_diag.cpp is a source file that builds a commandline harness for JS8 whitening/noise estimation

- llr_frame_benchmark.cpp runs paired-seed complex-bin AWGN trials through the
  actual 174-bit encoder, likelihood, BP, feedback and rescue code (sync and
  waveform extraction are assumed perfect). It reports complete correct-frame
  decodes as CSV. See its build command and arguments at the top of the file.
  Use its "sweep" argument for the 6x6 fixed-scale/erasure search, "final" for
  confirmation, and a final "coherent" argument to enable pilot-based blending.
  Its SNR is complex matched-bin signal power / noise power in dB.

- whitening_diag.cpp --llr-calibration-final TRIALS FIRST LAST STEP sweeps
  the waveform decoder around its 50% point with paired AWGN seeds; the
  --llr-calibration-gate variant compares the original phase RMS gate with
  0.20/0.60 and 0.25/0.75 rad. Its SNR is the real-sample synthesizer's dB
  setting (not the matched-bin SNR used by llr_frame_benchmark). Both report
  only decoded frames whose payload matches the transmitted message.

Acquisition/tracking correctness and replay (run from the repository root):

- decoder_waveform_test.cpp validates deterministic AWGN, measured power,
  full-band/2500-Hz SNR conversion, missing frame edges and PCM headroom.
  Build: clang++ -std=c++20 -O2 tools/decoder_waveform_test.cpp -o /tmp/waveform_test
  Run: /tmp/waveform_test

- acquisition_test.cpp exercises the real decoder's pilot coverage, sample
  bounds, missing-symbol noise/LLR neutrality, wrapped/partial input, stale
  scratch, tracker sign, timing correlator and exact aided replay in all modes.
- acquisition_diag.cpp runs all modes with bin-center/half-bin, early/late,
  drift and partial-window scenarios through the production decoder. It emits
  CSV with exact recovery, wrong-payload outputs, post-Costas-gate candidate
  counts, frequency/timing errors, runtime and measured SNR. Combining is
  disabled so repeated trials cannot reuse evidence. Use identical seeds and
  trial counts for baseline/candidate comparisons. Its drift frequency error
  is referenced to the frame midpoint; existing ordinary decoder output still
  reports its bulk acquisition frequency, not a smoothed midpoint estimate.

  Build either tool (replace TOOL and adjust the generated moc path):
    clang++ -std=c++20 -O2 -I. -Ibuild/JS8Call_autogen/include \
      $(pkg-config --cflags Qt6Core fftw3f) tools/TOOL.cpp \
      JS8_Mode/FrequencyTracker.cpp \
      build/JS8Call_autogen/33ZB6LRKMI/moc_JS8.cpp \
      $(pkg-config --libs Qt6Core fftw3f) -lpthread -o /tmp/TOOL
  Run: /tmp/acquisition_test
       /tmp/acquisition_diag 100 -18
       /tmp/acquisition_diag 32 -18 8192

  acquisition_diag accepts [TRIALS [SETTING_DB [FIRST_SEED [impairments]]]], defaulting to
  2, -24, 1234. The optional seed start supports separate validation sets.
  JS8_ACQUISITION_TRACE=A:late enables per-seed decoder logging for that
  mode/scenario only; traced timings must not be used for CPU comparisons.
  Appending "impairments" adds negative drift, fractional arrival and ±300-ppm
  sample-clock cases. Clock error stretches the transmitted waveform itself.
  Its frequency-error truth uses drift relative to frame start, including the
  clock-induced frequency scaling (the older start-time offset was incorrect).

  Checkpoint D controls (presence disables a feature, including a value of 0):
    JS8_DISABLE_PILOT_SMOOTHING=1   disables model-based frequency/drift steering
    JS8_DISABLE_FRACTIONAL_TIMING=1 disables pilot timing/clock steering
  The legacy trackers remain available as the conservative fallback. Disable
  both for a C observation-path comparison; the aided-SIC context fix is always
  active. New models are only accepted after pilot confidence/coverage checks
  and fresh-extraction validation. See acquisition_progress.md for limits.

  pilot_refinement_test.cpp is a standalone test of the pilot model and shared
  fractional extractor, including both drift/clock signs and aliased tone phase.
  Build: clang++ -std=c++20 -O2 -I. tools/pilot_refinement_test.cpp -o /tmp/pilot_test
  Run: /tmp/pilot_test

  Checkpoint A/B reference: build acquisition_diag with all three switches:
    -DJS8_BENCHMARK_LEGACY_FINE_SYNC
    -DJS8_BENCHMARK_SINGLE_TIMING_PEAK
    -DJS8_BENCHMARK_FREQUENCY_ONLY_SUPPRESSION
  Use the first switch alone to isolate timing-peak retention, or the last
  two to isolate fine-search changes. These are diagnostic build switches.
  They preserve the current sample-validity/tracker corrections.

  For exact score verification, build acquisition_test with:
    -DJS8_BENCHMARK_VERIFY_COARSE_SYNC -DJS8_BENCHMARK_VERIFY_RANK_SCAN
  Add -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer for sanitizer
  coverage. Test-local verification counters are provided by this harness;
  JS8_CPU_BENCHMARK is supplied only when building decoder_cpu_benchmark.

  decoder_waveform.h is shared by whitening_diag and acquisition_diag. The
  historical synthesis setting assigns noise variance 10^(-setting/10) and a
  unit-amplitude cosine (power approximately 0.5), so full-band sample SNR is
  approximately setting-3.0103 dB. SNR2500 adds 10*log10(6000/2500) dB. The
  generator now scales the complete signal+noise uniformly into PCM headroom
  without clipping and transmits only within the actual 79-symbol frame.
  Older overflowing waveform benchmark results are not directly comparable.
  Replay timings exclude construction but include waveform decoding/SIC;
  use decoder_cpu_benchmark for steady-state five-mode workload measurements.

- decoder_cpu_benchmark.cpp replays a deterministic 90-second 12-kHz ring
  containing tiled A/E WAV fixtures from media/tests, weak competing carriers,
  and fixed-seed noise through one persistent five-mode JS8 decoder worker.
  Its "noise" input tests ghost candidates without known transmissions; the
  optional fourth argument "autosync" exercises all five modes once per
  second (stress case). Build once normally and once with the five reference
  switches listed at the top of decoder_cpu_benchmark.cpp to compare all
  optimizations without changing any other decoder features. Both builds
  should include -DJS8_CPU_BENCHMARK for
  per-hotspot timing; the printed stage times include nested calls (e.g.,
  coherent scoring is part of tone-bin formation), while bpParts samples about
  one in 64 BP calls. Subtract nested stages rather than summing them. Its
  "diag" third argument counts attempted aided searches
  but enables debug formatting and must not be used for the CPU comparison.
  Output CPU duty is process CPU during decode calls divided by the virtual
  input duration, NOT a live receive-thread sample. The input is generated in
  advance; its virtual lag estimate does not simulate the UI busy queue.
  To isolate the BP edge lookup alone, compare two otherwise identical builds,
  one with -DJS8_BENCHMARK_LINEAR_BP_EDGES and one without.
  To isolate the check-centric BP loop, use
  -DJS8_BENCHMARK_ORIGINAL_CHECK_TO_BIT_LOOP in the reference build.
  To isolate timing-batched candidate scoring, use
  -DJS8_BENCHMARK_ORIGINAL_RANK_SCAN in the reference build;
  -DJS8_BENCHMARK_VERIFY_RANK_SCAN verifies every score against the original
  bit-for-bit and must not be used for CPU timings.

Decoder calibration overrides: JS8_LLR_SCALE (positive fixed multiplier,
default 2), JS8_LLR_ERASURE_THRESH (nonnegative, default 0),
JS8_COHERENT_GOOD_RMS_RAD / JS8_COHERENT_POOR_RMS_RAD (default 0.20/0.60).
For the previous normalized-likelihood baseline set JS8_LLR_SCALE=1,
JS8_LLR_ERASURE_THRESH=0.25, JS8_COHERENT_GOOD_RMS_RAD=0.12,
JS8_COHERENT_POOR_RMS_RAD=0.35, and JS8_LLR_FRAME_NORMALIZATION=1.
