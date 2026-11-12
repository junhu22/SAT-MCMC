# Supplementary comparison package — MSSP26-3659 (revision)

*Amplitude-Dependent Bayesian Operational Modal Analysis for Tall Buildings
under Typhoon Excitation: Field Evidence and Identifiability Limits*

**Data version.** Everything here rests on the re-estimated per-window set,
`results/per_window_modes_v5.csv`, which is the set every field number in the
manuscript is computed from. `VERSIONS.md` states this file by file, including
the two frequency-side calibration curves that are measured at the
translational configuration and are reported as applying there.

This package is the material referred to in the response to Reviewer 2's first
comment and in the manuscript's Data availability statement. It is submitted
**with this revision** rather than promised on acceptance, so that the comparison
it supports can be started now. It becomes public with the code on acceptance.

The paper makes no claim that its per-window spectral fit is more informative
than FDD, SSI-COV or the Bayesian FFT method. What it offers sits above that
layer. A ranking of the per-window estimators would nonetheless be a reasonable
thing to want, and it requires each method to be calibrated under the same
synthetic conditions before its numbers are set beside anyone else's. That
requirement, and the materials for meeting it, are what this package contains.

## What is here

    README.md                    this file
    COMPARISON_PROTOCOL.md       metrics, mode-identity rule, admissible
                                 wording, and the fairness clauses F-1 to F-4
    VERSIONS.md                  which data version each file rests on and
                                 which script produced it

    results/
      per_window_modes_v5.csv        the re-estimated per-window frequency and
                                     damping: the set every field number in the
                                     manuscript is computed from
      per_window_modes_release_v5.csv  the same estimates with the damping
                                     de-biased and its interval, the measured
                                     bias ratio per mode, and window RMS
                                     amplitude, by anonymised window index.
                                     Damping cells the acceptance conditions
                                     withhold are left empty; frequencies are
                                     always present, with their verdict in
                                     frequency_valid_<mode> beside them. See
                                     VERSIONS.md for why the two differ

    baselines/
      bfft_gate.py               synthetic generator: a single-mode resonance
                                 realised at a stated damping ratio
      fdd.py                     FDD implementation
      ssi_cov.py                 SSI-COV implementation

    .audit/
      recipe_v3.py               the per-window spectral fit used for the field
                                 estimates
      fitguard.py                the three acceptance conditions deciding
                                 which windows contribute a damping estimate
      freq_calibration.py        bias curve: frequency against amplitude, at
                                 the translational configuration
      calibration_curve_v5.py    bias curve: damping at six true ratios,
                                 0.006 to 0.030, at the rotational
                                 configuration
      alpha_xi_calibration.py    recovery of the amplitude coefficient and the
                                 curvature exponent
      field_damping_table.py     how the reported field damping follows from
                                 the per-window file
      *.json                     our own measured curves, as produced by those
                                 scripts

## Layout

The directory names mirror the repository because the scripts derive their paths
from `__file__` and read and write relative to it. Run them from the package
root and they work unedited. They are the scripts that produced the numbers,
frozen as they were; we would rather ship an awkward directory name than an
edited script.

## How to enter a method into the comparison

1. Run your estimator on the synthetic set defined in `COMPARISON_PROTOCOL.md`
   §6, clause F-1: ξ ∈ {0.006, 0.010, 0.014, 0.018, 0.024, 0.030}, 1800 s at
   32 Hz.
   `baselines/bfft_gate.py` generates it.
2. Report your own ξ̂/ξ_true against the parameters the generator actually
   realised, not against the nominal input. Clause F-1 applies to us on the same
   terms, and our own ratios are in the JSON files here.
3. Then the per-window comparison of §2 is defined, and the sentence describing
   its result is taken from the ladder in §5 — which was fixed before any
   measurement, so that a conclusion cannot be chosen after the fact.

The expected outcome is the middle of that ladder rather than the top, and that
is not a failure: per-window precision is not what this paper claims to improve.

## What is not here

The field acceleration records. They cannot be shared under the confidentiality
agreement covering the monitored structure. Everything in this package is either
synthetic or a processed modal estimate carrying an anonymised window index, and
no file identifies the building, its location, or the instrumentation owner.

No ranking of methods is included, because at the time of writing no baseline
has reported its own F-1 ratio. That is a statement about what has been
measured, not about what the answer would be.
