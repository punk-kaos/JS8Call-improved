# Acquisition and tracking — checkpoints A through D

## Implemented in checkpoints A/B

- Deterministic continuous-phase waveform generation with measured signal/noise
  power, explicit 2500-Hz reference conversion, and uniform PCM headroom scaling.
- Separate received input length, generated baseband length, and frame length.
- Symbol windows reject missing observations instead of clamping onto unrelated
  samples. Missing data supply zero LLRs and are excluded from noise medians.
- Fine synchronization uses the observed baseband span and requires two complete
  Costas blocks. Its score is scaled to the equivalent three-block expectation;
  this equalizes the mean noise scale, not the full statistical distribution.
- Coarse spectra exclude unreceived input windows, clear unused columns, and
  reject degenerate normalization. Candidate scratch tails are reset.
- Positive physical frequency residuals receive negative-exponent correction.
  Coherent compensation, cached rotations, aided replay and midpoint reporting
  use the same convention. Tone zero uses its wrapped negative-frequency bin.
- Timing correlation now matches the forward-FFT positive-frequency tone.
- Extraction metadata is initialized and recorded independently of coherent
  likelihood enablement.

The FFT narrowband filter remains circular. Observation bounds exclude padding,
but do not constitute an edge-transient guard or a partial-symbol likelihood
model. The current conservative policy erases incomplete symbol windows;
edge-transient modeling needs further evaluation before changing that policy.

## Initial paired replay

Eight deterministic seeds, seven scenarios, five modes: 280 receptions per SNR.
The same validated generator was compiled against both the original production
decoder and this milestone. Soft combining was disabled. One payload and one
carrier phase were used; this is a regression experiment, not a sensitivity
curve or a statistically powered dB-gain claim.

| Synthesis setting | Original exact recoveries | Milestone exact recoveries | Wrong payloads, either variant |
| --- | ---: | ---: | ---: |
| -12 dB | 218 / 280 | 234 / 280 | 0 |
| -18 dB | 94 / 280 | 99 / 280 | 0 |
| -24 dB | 6 / 280 | 7 / 280 | 0 |

The synthesis setting is approximately 3.0103 dB above full-band sample SNR and
approximately 0.7918 dB below the 2500-Hz-referenced injected SNR. See
`tools/README.txt` for build/run instructions and the diagnostic CSV fields.
Aggregated improvement does not establish improvement in every scenario.

Validation includes all-mode deterministic acquisition/tracking tests under
AddressSanitizer/UndefinedBehaviorSanitizer, existing coherent/aided/BP tests,
and 90-second five-mode CPU replays. The noise replay produced zero outputs;
that duration alone cannot establish a long-term false-output rate.

## Checkpoint C — implemented

Fine synchronization now derives its grid from the seven-symbol pilot duration
and its capture width from the mode's coarse FFT spacing:

| Mode | Initial residual capture | Coarse step | Fine step |
| --- | ---: | ---: | ---: |
| A / Normal | ±2.5 Hz | 0.5 Hz | 0.1 Hz |
| B / Fast | ±3.0 Hz | 0.5 Hz | 0.1 Hz |
| C / JS8 40 | ±6.0 Hz | 1.0 Hz | 0.2 Hz |
| E / Slow | ±1.0 Hz | 0.2 Hz | 0.04 Hz |
| I / JS8 60 | ±10.0 Hz | 2.0 Hz | 0.4 Hz |

- At most two coarse neighborhoods receive joint local refinement. A winning
  frequency-edge hypothesis gets one extension of two coarse steps. All pilot
  references, including extension/refinement grids, are cached per mode.
- Timing remains within the original quarter-symbol search envelope. A local
  timing move needs ≥5% total pilot-power improvement and >1% improvement in
  every observed block, with at least two observed blocks. This avoids selecting
  a new symbol boundary based on one noisy block or silently widening timing
  capture. Earlier unrestricted refinement regressed the Normal late-start
  case; constraining the envelope corrected that regression.
- Up to three separated timing maxima survive per frequency bin. Duplicate
  suppression now requires proximity in both frequency and time (one symbol).
- The normalization reference remains the 40th percentile of exactly one
  maximum per frequency bin. Secondary timing peaks do not change that sample.
- Existing limits of 300 normal and 24 deep candidates remain. A separate cap of
  24 secondary timing hypotheses per mode/pass limits the additional FEC work.
  Exhausted secondary entries are discarded without suppressing primary peaks.

### Paired validation on additional seeds

32 seeds starting at 8192, seven scenarios and five modes: 1120 receptions per
setting and variant. The A/B reference switches reproduce the frozen A/B binary
on the original replay seeds (payload/candidate counts and synchronization
errors). One payload and one phase are still used; these are acquisition
regression measurements, not a measured sensitivity-curve shift.

| Synthesis setting | A/B exact recoveries | A/B/C exact recoveries | Wrong payloads, either variant |
| --- | ---: | ---: | ---: |
| -12 dB | 909 / 1120 | 979 / 1120 | 0 |
| -18 dB | 407 / 1120 | 413 / 1120 | 0 |
| -24 dB | 25 / 1120 | 26 / 1120 | 0 |

At -12 dB, JS8 40 improves from 174 to 200 correct recoveries and JS8 60 from
98 to 140 (224 receptions per mode). Normal and Fast totals remain 209 and 208.
Individual scenario results remain mixed: JS8 40 early-start recovery drops
from 16/32 to 11/32, despite its aggregate gain. At -18 dB, Normal loses two
recoveries overall. Further timing work must evaluate these cells explicitly;
aggregate gain does not establish uniform improvement.

### Runtime and correctness

- 378 checks pass with ASan/UBSan and both score-verification switches: full-bin
  capture, bounded extension, timing-envelope preservation, plateau handling,
  normalization stability, candidate budgets, overlapping equal-power frames
  separated by seven symbols, and original-versus-cached score equality.
- Seven-symbol separation validates independent pilot peaks during overlapping
  frames. It does not promise arbitrary same-frequency overlap: an exploratory
  three-symbol-separation case in JS8 60 remained below the acquisition gates.
- 90-second, five-mode autosync noise replay: zero outputs in both variants;
  CPU duty rises from 10.28% to 11.89%.
- Crowded replay: identical 72-output payload/mode fingerprint in both variants;
  CPU duty rises from 25.82% to 36.95%, with maximum pass time increasing from
  442 to 623 ms. Neither replay exceeds its virtual one-second deadline. These
  measurements do not simulate the live UI receive queue.
- Before the additional-hypothesis cap, crowded CPU duty reached about 48%;
  the cap reduces cost without changing the recorded sensitivity replay counts.

## Checkpoint D — implemented

### Pilot-supported tracking and extraction

- `pilot_refinement.h` fits residual frequency and linear drift from all three
  Costas blocks, using noise-scaled observations. Data/codeword guesses never
  enter this fit. Global timing is searched at quarter-sample steps within
  ±one-eighth of a symbol, capped at four downsampled samples.
- Block-wise timing/phase slopes refine the fractional arrival estimate. An
  affine symbol-clock slope is retained only when resolved above three sigma,
  with a 1000-ppm limit and a block-fit residual limit of 0.10 sample.
- At least 15 valid pilots, four winning tones in each block, median pilot
  power/noise ≥4, phase RMS ≤0.30 rad and coherence ≥0.95 are required. Each
  block's locally measured frequency must agree with the global trajectory
  within its measured uncertainty. Drift is bounded to 0.15 Hz/s and residual
  frequency to 0.45 tone spacing at the observed pilots.
- A proposed model re-extracts the existing observation rather than adding it
  as an independent reception. Re-extracted pilots must pass the gates again;
  observed data coverage must remain intact. Pilot contrast may fall by at
  most 5% in any block and 1% overall for a timing move. Frequency-only moves
  must improve total contrast by 0.5%. Rejected proposals leave C's bins and
  metadata untouched.
- `fractional_window.h` supplies a shared 16-tap Lanczos/sinc extractor. Its
  passband is centered at 3.5 tones so aliased tones 6/7 in 12-sample modes
  retain the positive-tone phase convention. Integer starts take an exact-copy
  path; fractional starts require their complete filter halo.
- Nominal starts and displacements remain double precision through coherent
  scoring and aided replay. Applied per-symbol frequency corrections are
  recorded exactly. Within-symbol drift curvature is negligible at the stated
  bound; the model supplies each window's midpoint frequency rather than
  applying a separate quadratic derotation inside that window.
- Reporting uses the carrier estimate from fresh pilots. Cancellation now
  receives both the accepted frequency and time explicitly; the previous
  hidden `xdt2` reference seed was incorrect after aided timing changes.

### Expanded paired replay

16 new seeds starting at 16384, 11 scenarios and five modes: 880 receptions
per setting/variant. The new cases include negative drift, fractional arrival,
and ±300-ppm clock error. The waveform generator stretches actual symbol
boundaries and carrier phase for clock error. Frequency-error truth now uses
drift measured from frame start, correcting the earlier diagnostic offset.

| Synthesis setting | Checkpoint C exact recoveries | Checkpoint D exact recoveries | Wrong payloads, either variant |
| --- | ---: | ---: | ---: |
| -12 dB | 751 / 880 | 751 / 880 | 0 |
| -18 dB | 310 / 880 | 316 / 880 | 0 |
| -24 dB | 15 / 880 | 16 / 880 | 0 |

At -12 dB, mean arrival error among successful late-start decodes decreases
from 9.18 to 2.09 ms; drift-case arrival error decreases from 16.97 to 9.26 ms.
Negative-clock-case frequency error decreases from 0.526 to 0.393 Hz. These
are aggregate successful-decode errors, not comparisons restricted to an
identical set of successful trials. Recovery gains are modest, not a measured
large sensitivity shift. Individual results remain mixed: JS8 60 clock-plus
recovery drops from 9/16 to 8/16 at -12 dB, offset by a gain elsewhere.

### Validation and runtime

- Pure pilot/window tests exercise all mode rates, both drift/clock signs,
  positive-tone phase through aliasing, integer-copy equivalence, filter-halo
  bounds, incoherent/weak/duplicate/nonfinite pilots and block phase jumps.
- Production acquisition tests exercise real FSK samples, accepted fractional
  starts, bit-exact aided replay and unchanged rejected proposals. Existing
  coherent/aided/SIC tests also pass.
- Disabling both new features reproduces the frozen C replay counts and
  synchronization errors on the checked subset. The independent aided-SIC
  context correction remains active.
- Paired 90-second five-mode noise replay: zero outputs in both variants;
  CPU duty is 12.71% for C and 13.74% for D. Crowded duty is 43.81% versus
  44.57%, with maximum pass times 755 versus 748 ms and no virtual deadline
  overruns. These timings are from the same paired run; compare neither with
  earlier sessions' absolute timings nor with a live UI deadline measurement.
- Crowded raw output count changes from 72 to 91. This workload's count and
  fingerprint are diagnostics, not an exact-payload accuracy oracle.

## Remaining work

1. Investigate remaining early-start/bin-offset and clock-error regression
   cells, without tuning against only one message/phase.
2. Broader payload/phase/seed coverage, sensitivity curves, recorded-audio and
   long-duration interference/noise validation.
