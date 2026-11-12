# Comparison protocol and fairness clauses

*Supplementary material to MSSP26-3659 (revision). This document is the English
rendering of the internal pre-registration that was frozen before any comparison
was attempted; the internal copy is `.audit/h2h_preregistration.md` in the code
release. Nothing here was written after seeing a result.*

**Revision note, 2026-08-27.** The rules in this document are unchanged: the
method set M0-M3, the metrics P1-P3, the fairness clauses F-1 to F-4 and the
wording ladder T0-T3 are as they were frozen. What has changed is the data
version they point at and the band half-width they are applied at. The
rotational mode is now estimated on a 60 mHz half-width rather than 20 mHz,
which is a per-mode setting given in revised Table 3 of the manuscript, and
every file in this package rests on the re-estimated per-window set
(`VERSIONS.md`). Clause F-1's damping levels are also restated at the six that
were actually measured, in place of three that were named before the
measurement existed. None of these changes alters what a method must do to
enter the comparison or what wording a result would license.

The purpose of this document is to state, in advance, what a comparison between
this paper's per-window estimator and established modal identification methods
would have to satisfy in order to be fair, and what wording would then be
admissible. It is released with the revision so that the comparison can be
carried out by anyone, including by a reader who disagrees with us.

---

## 1 Methods

| ID | Method | Returns | Uncertainty |
|---|---|---|---|
| **M0** | FDD — frequency domain decomposition, peak picking on the first singular value of the PSD | f_i | none |
| **M1** | SSI-COV — covariance-driven stochastic subspace identification | f_i, ζ_i | stabilisation-diagram scatter, which is not a posterior |
| **M2** | Bayesian FFT / BAYOMA, per window | f_i, ζ_i | asymptotic posterior c.o.v. |
| **M3** | AD-BOMA two-step (this paper) | f_{i,0}, α_i^f, γ_i^f, and the damping pair | curve-fit covariance |

M0–M2 are to be run on the same 141 clean windows, the same 15 channels and the
same spectral settings; M3's existing results are used as they stand.

Two window counts appear in this work and they are not interchangeable.
**Frequency comparisons use the 141 clean windows.** **Damping comparisons use
140**: window 144 is a 10.7-minute window at n_avg = 1, outside the coverage of
the damping calibration, and is excluded from every damping estimate while
being retained for frequency. A method entering the comparison should report
each quantity on the count that applies to it and say which it used.

## 2 Metrics, frozen before the comparison

The reviewer's question is whether the treatment yields materially more
informative estimates. That is operationalised as three quantities, **reported
separately and never combined into a single score**.

**P1 — per-window spread (primary).** For each mode and each method, take the
per-window frequency series, remove that method's own amplitude trend, and
report the robust standard deviation (MAD × 1.4826) of the residual. Smaller is
more precise. This is the only directly comparable quantity: every method
produces per-window frequencies.

**P2 — agreement of the fitted amplitude trend (primary).** Fit the same power
law `f(a) = f_0[1 − α(a/a_ref)^γ]` to each method's per-window estimates and
compare (f_0, α, γ), each with its own uncertainty, stating how that
uncertainty was obtained.

The earlier text required every method to report a contraction ratio κ_j for α
and γ. That requirement is withdrawn: κ is defined against a posterior, and
M0, M1 and M3 do not produce one, so the metric was not applicable to three of
the four methods it was being asked of. **κ_j is now optional**, reported by a
method that has a posterior to define it against, and its absence is not a
failure to comply.

**P3 — interval coverage (secondary, M2 against M3 only).** Empirical coverage
of the nominal 90% interval on held-out windows. M0 returns no uncertainty and
M1's stabilisation scatter is not an uncertainty, so neither participates.

## 3 FDD returns no damping

Peak picking on the first singular value locates frequency only. A half-power
estimate can be forced onto it, but at the 42 and 33 mHz spacings of the
closely-spaced cluster it is contaminated by the neighbouring peak and is a
fitting artefact rather than a damping ratio. Therefore M0 takes no part in any
damping comparison, its damping cells read `n/a` rather than being left blank or
filled with zero, and any sentence involving M0 is restricted to frequency
estimates.

## 4 The closely-spaced band, and mode identity

TX2a (0.851 Hz), TX2b (0.893 Hz) and TY2 (0.926 Hz) share the 0.80–0.97 Hz band
at spacings of 42 and 33 mHz. The four methods handle that band differently and
the difference must be stated in any comparison table, or the comparison is not
valid:

| Method | Treatment of the band | Risk |
|---|---|---|
| M0 FDD | three-peak pick on the first singular value; if peaks merge, recorded as unresolved and **not** interpolated | peaks merge at high amplitude |
| M1 SSI-COV | model order must admit three modes; identity assigned by stabilisation-diagram clustering | identity exchange |
| M2 Bayesian FFT | three modes fitted jointly within the band | — |
| M3 AD-BOMA | trace-PSD projection and multi-Lorentzian fit | — |

**Mode-identity rule, fixed before running.** Anchor labels on the low-amplitude
windows (1–10) and bind them in ascending frequency order. If two adjacent modes
cross in frequency in a given window, that window is marked
`identity_ambiguous` for that cluster and excluded from P1 and P2, with the
number excluded reported. This paper's own data show a modal assurance criterion
of 0.071 for TX2a at high amplitude, so identity exchange is a known fact rather
than a hypothetical, and it must not be allowed to contaminate a baseline
comparison.

## 5 Admissible wording, chosen from this table only

Once the measurement exists, the sentence describing it must be taken from this
ladder. The ladder was fixed before the measurement so that a conclusion cannot
be selected after the result is seen.

| Tier | Condition | Admissible wording |
|---|---|---|
| **T3** | M3 better than every baseline on P1, with non-overlapping bootstrap 95% intervals | "yields materially tighter per-window estimates than X, Y and Z" |
| **T2** | M3 better than most baselines on P1, intervals overlapping | "is comparable to or tighter than", with the intervals shown |
| **T1** | no distinguishable difference on P1 | "performs on par with established methods, while additionally providing the amplitude-dependence parameters that none of them estimate" |
| **T0** | M3 worse than some baseline on P1 | report it as measured, and state why the method is still useful (P2 supplies a quantity the baselines do not) |

**T1 is the expected outcome and is not a failure.** What this paper offers is
P2 — the amplitude law and its sequential update — not per-window precision.
Describing a T1 result in T3 language is precisely what the reviewer's comment
guards against.

**Prohibited after the run:** adding a metric, changing the detrending method,
changing the robust statistic, or changing the window set. Any such change must
be recorded as a deviation, with its motivation, and cannot be made silently.

## 6 Fairness clauses

These arose from a measured fact rather than from principle. On a single-segment
periodogram (n_avg = 1), the maximum-likelihood damping estimate of a
single-degree-of-freedom resonance carries roughly a **25% negative
finite-sample bias in the damping estimate itself**, decaying approximately as
1/n_avg. (Section 7 reports a second figure that is also about 25%. The two are
different quantities and should not be read as one: this one is a bias in an
estimate, that one is the optimism of an interval.) That is a property of the spectral
configuration, not of any one pipeline: **any** method fitting on the same
configuration inherits it. A comparison that does not control for it is
comparing choices of segment length, not methods.

**F-1 — every method reports its own synthetic bias before entering the
comparison.** Any method must first report its own ξ̂/ξ_true on the same
synthetic set: ξ ∈ {0.006, 0.010, 0.014, 0.018, 0.024, 0.030}, 1800 s at 32 Hz.
These are the six levels at which the curve in this package was actually
measured; the three named in the earlier text were chosen before the
measurement existed and are superseded by it. The
true value is recovered from the parameters the generator actually realised, not
from the nominal input. The ratio is stored alongside the method's results, and
**a method without it does not enter the comparison.**

**F-2 — a single-segment configuration must be declared.** An EM Bayesian FFT
implementation working on a single FFT segment is at n_avg = 1. If it is
configured that way, the inherited finite-sample bias must be declared
explicitly; its own measured ratio must be reported under F-1; and it may not be
placed beside an n_avg = 6 result as evidence of method quality until both are
brought to the same n_avg.

**F-3 — the configuration table travels with the comparison table.** Any
comparison table carries a column for n_avg and a column for each method's own
synthetic ratio, so that a reader can see at once whether a difference is a
method difference or a configuration difference.

**F-4 — these clauses bind us equally.** AD-BOMA and the two-step estimator are
subject to F-1 on the same terms. We may not report a baseline's bias without
reporting our own.

F-1 to F-4 apply uniformly to FDD, SSI-COV, EM Bayesian FFT, the two-step
estimator and AD-BOMA.

## 7 Status of this paper's own entry

The materials in this package include our own measured curves, so that F-1 and
F-4 are satisfied on our side before anything is asked of anyone else:
frequency bias against amplitude, recovery of the amplitude coefficient and the
curvature exponent, the joint frequency-and-damping recovery, and the damping
bias curve across the six true damping levels of F-1. One of them runs against us and is
reported with the rest: the posterior coefficient of variation the fit reports
is about **25% smaller than the scatter the same fit shows across seeds** —
that is, the interval is optimistic by roughly that much. It is a different
quantity from the 25% finite-sample bias of Section 6, which is a bias in the
damping estimate rather than an optimism in an interval; the numerical
coincidence is not a relationship.

No ranking is presented in the paper or in this package, because **no baseline
has reported an F-1 ratio on this synthetic set.** That is a statement about
what has been measured, not about what the answer would be.

One qualification, so the sentence is not read as stronger than it is. Bayesian
FFT implementations do publish damping recovery on synthetic data, and we do
not claim otherwise. What is absent is a ratio on *this* set: the published
measurements we are aware of are at damping levels and a spectral configuration
that do not match the six levels of F-1, and a ratio measured elsewhere cannot
be placed in an F-3 table beside one measured here. The clause asks for a
number on a stated set, not for a first-ever measurement.
