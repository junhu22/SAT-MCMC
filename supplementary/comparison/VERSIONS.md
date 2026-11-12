# Data versions in this package

Every file below states which per-window estimate set it rests on and which
script produced it. This page exists because an earlier build of this package
did not say, and shipped the superseded set without anyone noticing.

**The reported set is `per_window_modes_v5.csv`.** Every amplitude-dependence
coefficient, every damping figure and every post-event comparison in the
manuscript is computed from it. Where a file below rests on something else, the
reason is given.

## What changed, and when

The rotational mode R1 was re-estimated on a wider band during the revision.
Its earlier band, 20 mHz, was narrower than the mode's own 30.4 mHz drift
across the storm, so at high amplitude the resonance left the band and the fit
returned damping values between 10⁻⁸ and 10⁻⁸². Those values passed a
positivity test and entered the amplitude--damping regression as near-zero
damping at the storm peak. Re-estimating R1 at 60 mHz changed its damping
coefficient from unresolved to the best determined of the three modes.

Nothing about TX1 or TY1 changed: their rows are identical between the two
sets. The change is confined to R1, and to every quantity computed from it.

## Per file

| File | Data version | Produced by | Note |
|---|---|---|---|
| `results/per_window_modes_v5.csv` | v5 | `.audit/refit_r1_v5.py` | the reported set; 145 windows, 141 after outlier rejection |
| `results/per_window_modes_release_v5.csv` | v5 | `.audit/build_release_v5.py` | de-biased damping, empty where the guard withholds it; frequencies are always present, with the verdict in `frequency_valid_*` beside them |
| `.audit/recipe_v3.py` | — | frozen | the per-window fit itself. It takes the band half-width as an argument and is used at 20 mHz for TX1 and TY1 and 60 mHz for R1; the name is historical and does not mean it is tied to the v3 set |
| `.audit/fitguard.py` | — | frozen | the three acceptance conditions: resolution floor, half-power points inside the band, and at or below the highest calibrated damping level |
| `.audit/calibration_curve_v5.py` | — | frozen | six damping levels, 0.006 to 0.030, at the **rotational configuration**: reference 0.447 Hz, half-width 60 mHz, 96 realisations |
| `.audit/calibration_curve_v5.json` | v5 configuration | the script above | the six measured levels: ratios 0.994, 0.947, 0.966, 0.969, 0.973, 0.970 |
| `.audit/calibration_curve_low.py` / `.json` | — | frozen | two further levels, 0.004 and 0.005, **sampled but not resolved**; see the note below |
| `.audit/field_damping_table.py` | takes its input as an argument | frozen | run as `field_damping_table.py results/per_window_modes_release_v5.csv .audit/calibration_curve_v5.json field_damping_table_v5.json` |
| `.audit/field_damping_table_v5.json` | v5 | the script above | the table the manuscript reports; records the source and calibration paths it actually read |
| `.audit/freq_calibration.py` / `.json` | — | frozen | frequency bias against amplitude, measured at the **translational configuration**: a single synthetic mode at 0.247 Hz with a 20 mHz half-width. See the note below |
| `.audit/alpha_xi_calibration.py` / `.json` | — | frozen | recovery of the amplitude coefficient and curvature exponent, same configuration and the same note |
| `baselines/bfft_gate.py` | — | frozen | the synthetic generator both calibration scripts draw from |
| `baselines/fdd.py`, `baselines/ssi_cov.py` | — | frozen | baseline implementations, no field data |

## A limitation of the two frequency-side calibrations, stated plainly

`freq_calibration.py` and `alpha_xi_calibration.py` simulate one synthetic mode
at 0.247 Hz and fit it at a 20 mHz half-width. That is the translational
configuration, and it matches TX1 and TY1. It does not match R1, which is
estimated at 60 mHz.

We have not re-measured those two curves at the rotational configuration, and
we do not claim they transfer to it. What is measured at R1's configuration is
the **damping** bias curve, `calibration_curve_v5`, which is the one the
reported damping table uses and the one R1's own bias band comes from. The
frequency-side curves are reported for the configuration they were measured at
and should be read as applying there.

This is a gap in coverage, not a defect in the numbers: each curve is correct
for the configuration it names. We prefer to name it than to widen the claim.

## The two lowest damping levels are sampled, not certified

The reported curve is measured at six levels from 0.006 to 0.030. The
resolution floor for a damping estimate is 0.004, so between 0.004 and 0.006
there are field windows that are resolved and uncalibrated.

Two further levels, 0.004 and 0.005, were run on the same generator, the same
fit, the same 60 mHz half-width about the same 0.447 Hz reference, 96
realisations each. They are shipped in `calibration_curve_low.json` and give
ratios of 1.036 and 0.968, with 95 per cent intervals of [0.987, 1.097] and
[0.926, 1.024].

**We report them as sampled rather than certified, and make no flatness claim
below 0.006.** On that segment there are three points, and a regression of the
ratio on the true damping gives a slope of -21.4 with a standard error of 26.9.
A departure large enough to matter -- the ratio moving 0.05 across 0.004 to
0.006 -- corresponds to a slope of 25, which that standard error cannot
resolve. The test would pass whether the segment were flat or not, so passing
it establishes nothing.

Two further cautions belong with the number. The ratio at 0.004 is 1.036, the
only level in the eight above unity, against a mean of 0.970 for the other
seven; its interval includes 1, so it is not significant on its own terms, but
it is the point nearest the floor and it behaves unlike the rest. And the
band-centre offset in this design is assigned by level index, 0 mHz at the
lowest rising to 30 mHz at the highest, so the two new levels carry no offset
while the upper five carry up to 30; a slope fitted across all eight is partly
a slope on offset rather than on damping.

None of this changes the correction applied to the reported table. That is
interpolated at each mode's median damping -- 0.0068, 0.0068 and 0.0105, all
above 0.006 -- and gives ratios identical to five decimals whether computed on
the six-level curve or the eight-level one.

## Frequency and damping are released under different conventions

A per-window **frequency** is always present in the released file, and the
acceptance verdict for it sits beside it in `frequency_valid_<mode>`. A
per-window **damping** value is absent where the acceptance conditions withhold
it.

The asymmetry is deliberate. The two frequency conditions are diagnostic: they
are what found the band misspecification for the rotational mode, and on the
widths now adopted exactly one of the 141 x 3 = 423 mode-window pairs fails
them. The paper's own Step 2 fits use every finite frequency rather than the
filtered subset, so blanking the flagged value would have made this file unable
to reproduce the paper. Keeping the value and the verdict together lets you
reproduce either: the published fits, by using every finite frequency, or the
filtered set, by dropping rows where the flag is False. That one window moves
the amplitude coefficients by at most 0.034 of a standard error.

Damping is different. A withheld damping value is one the bias calibration
cannot correct -- below the resolution floor, with half-power points outside
the band, or above the highest calibrated level -- so it is withheld from the
estimate rather than flagged and used.

## Reproducing the reported damping table

    python .audit/calibration_curve_v5.py
    python .audit/field_damping_table.py \
        results/per_window_modes_release_v5.csv \
        .audit/calibration_curve_v5.json \
        field_damping_table_v5.json

The table body it prints is the one in the manuscript: TX1 and TY1 at
alpha_xi 0.00126 and 0.00070 with z of 8.6 and 4.0, R1 at 0.00247 with z 13.7,
on 126, 129 and 139 contributing windows out of a 140-window base.
