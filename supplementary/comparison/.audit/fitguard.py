"""Domain guards for per-window modal fits.

Three concepts were collapsed into one flag. recipe_v3 sets

    out['ok'] = bool(best.success or np.isfinite(best.fun))

and every consumer read that as "this estimate is usable". It is not: it says
the objective evaluated to a finite number. Thirty-seven invalid estimates
carried ok = True, including R1 at 54.8 and 87.4 Hz against a mode at 0.447.

The three are separated here and never recombined:

  numerical_completion  what the optimiser reported. RECORDED ONLY. It is not
                        an input to any validity decision.
  fit_valid_frequency   the frequency estimate lies in its own band, on BOTH
                        the measured peak and the fitted value, and the fit
                        tracks the peak -- see below.
  damping_resolved      the damping lies in the physical range, above the
                        resolution floor for this record length, and the
                        resonance's own half-power points fall inside the
                        band.

Why the frequency gate takes two references. Window 54's measured peak is
outside the band while its fit sits inside; window 62 is the reverse. A guard
keyed on the fit alone misses the first, and one keyed on the peak alone misses
the second. Both are checked, and |f_hat - peak| is reported as a diagnostic
in its own right.

The edge margin is ONE linewidth, computed from the window's own fitted
damping: the half-power points of the resonance must fall inside the band, or
the width the fit reports is the width of the band rather than of the mode.

The margin governs DAMPING only. An earlier version let it veto the frequency
too, and 19 of 141 R1 windows failed -- not because the frequency was wrong but
because R1's damping at the amplitude peak is 0.023, so three linewidths came
to 61 mHz against a 60 mHz half-width. Frequency validity is decided by the two
things that actually bear on it: the estimate is in the band, and it tracks the
independently measured peak.

Verdicts, from the flags:

  frequency invalid                      -> the whole fit is invalid
  frequency valid, damping unresolved    -> keep the frequency, report no
                                            damping
  frequency invalid AND damping -> 0     -> degenerate solution, whole fit
                                            invalid, flagged as such
"""

XI_FLOOR = 0.004          # resolution floor at this record length
XI_LO, XI_HI = 0.001, 0.05        # physical range, Table 7
XI_CAL_HI = 0.030                 # highest calibrated level; above it the
                                  # bias ratio would be extrapolated
XI_STARTS = (0.004, 0.010, 0.020)  # recipe_v3.fit multi-start values
SEED_EPS = 1e-9
EDGE_LINEWIDTHS = 3.0
XI_DEGENERATE = 1e-4      # below this, log xi has run away
EDGE_LINEWIDTHS_DAMPING = 1.0     # half-power points must be inside the band


def linewidth(f, xi):
    """Half-power bandwidth 2 xi f, the natural scale for 'near the edge'."""
    return 2.0 * xi * f


def guard(f_hat, xi_hat, peak_hz, centre, half_width, ok_flag=None,
          xi_for_width=None):
    """Judge one window-mode fit. Returns a dict of flags and reasons."""
    res = {'numerical_completion': None if ok_flag is None else bool(ok_flag),
           'fit_valid_frequency': False, 'damping_resolved': False,
           'reasons': []}

    def finite(x):
        return x is not None and x == x and abs(x) != float('inf')

    # width used for the edge margin: the window's own estimate where it is
    # sane, otherwise the floor, so a runaway xi cannot buy itself margin
    xw = xi_for_width if finite(xi_for_width) else xi_hat
    if not finite(xw) or xw < XI_FLOOR:
        xw = XI_FLOOR
    margin = EDGE_LINEWIDTHS_DAMPING * linewidth(centre, xw)
    lo, hi = centre - half_width, centre + half_width

    fit_in = finite(f_hat) and lo <= f_hat <= hi
    peak_in = finite(peak_hz) and lo <= peak_hz <= hi

    if not finite(f_hat):
        res['reasons'].append('frequency not finite')
    elif not fit_in:
        res['reasons'].append('fitted frequency outside the band')
    if finite(peak_hz) and not peak_in:
        res['reasons'].append('measured peak outside the band')

    tracks = True
    if finite(f_hat) and finite(peak_hz):
        res['fit_minus_peak_mhz'] = (f_hat - peak_hz) * 1e3
        # the fit must sit on the peak that is there, to within the margin the
        # band itself allows
        tracks = abs(f_hat - peak_hz) <= max(margin, 0.5 * linewidth(centre, xw))
        if not tracks:
            res['reasons'].append('fitted frequency does not track the measured '
                                  'peak')

    res['fit_valid_frequency'] = bool(fit_in and peak_in and tracks)

    # the half-power points must lie inside the band, or the reported width is
    # the band's rather than the mode's. This governs damping only.
    half_power_inside = (finite(f_hat)
                         and lo + margin <= f_hat <= hi - margin)
    if finite(f_hat) and not half_power_inside:
        res['reasons'].append('half-power points reach outside the band')

    if not finite(xi_hat):
        res['reasons'].append('damping not finite')
    elif not (XI_LO <= xi_hat <= XI_HI):
        res['reasons'].append('damping outside the physical range')
    elif xi_hat < XI_FLOOR:
        res['reasons'].append('damping below the resolution floor')
    elif any(abs(xi_hat - s) < SEED_EPS for s in XI_STARTS):
        res['reasons'].append('damping equals an optimiser start value')
    elif xi_hat > XI_CAL_HI:
        res['reasons'].append('damping above the highest calibrated level; the '
                              'bias ratio would be extrapolated')
    elif not half_power_inside:
        pass
    else:
        res['damping_resolved'] = True

    degenerate = (not res['fit_valid_frequency'] and finite(xi_hat)
                  and xi_hat < XI_DEGENERATE)
    res['degenerate'] = bool(degenerate)
    if degenerate:
        res['reasons'].append('degenerate solution: frequency out of band with '
                              'damping driven to zero')

    if not res['fit_valid_frequency']:
        res['verdict'] = 'invalid'
    elif not res['damping_resolved']:
        res['verdict'] = 'frequency only'
    else:
        res['verdict'] = 'valid'
    return res
