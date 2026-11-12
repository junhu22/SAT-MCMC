"""AC-9's field damping table, computed from the released per-window file.

Three rulings this implements.

  window set   140 windows. Window 144 averages one segment in the
               16,384-point pipeline against six for the rest, and the debias
               curve is calibrated at n_avg = 6, so it falls outside the
               calibration's coverage. B_revision_list section 5: a damping
               number from the 141-window set is an error.
  uncertainty  bootstrap over windows, not the analytic OLS error. The two
               agree on the point estimates to 0.10 SE and disagree on the
               errors by up to a factor of two; the bootstrap is the
               conservative one and it is what the significance statement uses.
  entry point  the field damping values are computed here from
               results/per_window_modes_release.csv, not read from a constant.
               The hard-coded FIELD dict in calibration_curve.py is an
               evaluation point for the bias curve, not a measurement.

The bias ratio and its band are interpolated in log xi between the MEASURED
calibration levels of .audit/calibration_curve.json, which is how that file
computes them; nothing is extrapolated beyond the measured range.

The bootstrap resamples windows independently. It does not model the serial
correlation between adjacent windows, so it may still be optimistic; the
manuscript and the letter say so where they use it.

Writes .audit/field_damping_table.json and prints the LaTeX table body.
"""
import csv
import hashlib
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
REL = os.path.join(ROOT, *(sys.argv[1].split('/') if len(sys.argv) > 1
                           else ('results', 'per_window_modes_release.csv')))
CAL = os.path.join(ROOT, *(sys.argv[2].split('/') if len(sys.argv) > 2
                           else ('.audit', 'calibration_curve.json')))
DEST_NAME = sys.argv[3] if len(sys.argv) > 3 else 'field_damping_table.json'
REL_NAME = sys.argv[1] if len(sys.argv) > 1 else 'results/per_window_modes_release.csv'
CAL_NAME = sys.argv[2] if len(sys.argv) > 2 else '.audit/calibration_curve.json'
DEST = os.path.join(ROOT, '.audit', DEST_NAME)

MODES = ['TX1', 'TY1', 'R1']
W144 = 144
A_REF_UG = 100.0
NB = 20000
SEED = 20260825
# NB is 20,000 rather than the 2,000 used elsewhere because the quantity being
# reported is the RATIO alpha_xi / SE, and at 2,000 resamples that ratio moved
# by 0.54 across six seeds for TX1 -- enough to change the first decimal of a
# number the significance statement quotes. At 20,000 the seed-to-seed spread
# is under 0.2, so one decimal is defensible.


def xi_model(u, xi):
    """xi0 + alpha_xi * ln u, by OLS. Same form as .audit/window144_impact.py."""
    X = np.column_stack([np.ones_like(u), np.log(u)])
    b, *_ = np.linalg.lstsq(X, xi, rcond=None)
    return float(b[0]), float(b[1])


def boot(fn, u, y, rng, n=NB):
    out = []
    for _ in range(n):
        j = rng.integers(0, len(u), len(u))
        try:
            out.append(fn(u[j], y[j]))
        except Exception:
            pass
    return np.asarray(out)


def main():
    rows = [r for r in csv.DictReader(open(REL, encoding='utf-8'))
            if r['is_outlier'] != 'True' and int(r['window_id']) != W144]
    cal = json.load(open(CAL, encoding='utf-8'))
    xs, rs, los, his = [], [], [], []
    for k in sorted(cal['levels'], key=float):
        v = cal['levels'][k]
        xs.append(float(k)); rs.append(v['ratio'])
        los.append(v['ci_lo']); his.append(v['ci_hi'])
    lx = np.log(np.array(xs))

    rng = np.random.default_rng(SEED)
    amp = np.array([float(r['amplitude_ug']) for r in rows])
    u = amp / A_REF_UG

    out = {'schema': 'field_damping_table/1',
           'n_windows': len(rows), 'window_excluded': W144,
           'a_ref_ug': A_REF_UG, 'n_bootstrap': NB, 'seed': SEED,
           'source': REL_NAME,
           'calibration': CAL_NAME,
           'se_estimator': 'bootstrap over windows; no serial-correlation model',
           'calibration_levels': dict(zip(['%.4f' % x for x in xs], rs)),
           'modes': {}}

    print('%-5s %8s %9s %-16s %10s %-20s %11s %9s %7s'
          % ('mode', 'xi_hat', 'ratio', '95% band', 'xi_debias',
             '95% interval', 'alpha_xi', 'boot SE', 'z'))
    print('-' * 104)
    for m in MODES:
        xi = np.array([float(r['xi_%s_raw' % m])
                       if str(r['xi_%s_raw' % m]).strip() else np.nan
                       for r in rows])
        ok = np.isfinite(xi) & (xi > 0)
        xh = float(np.median(xi[ok]))
        if not (min(xs) <= xh <= max(xs)):
            raise SystemExit('%s: xi_hat %.5f outside the measured range' % (m, xh))
        ratio = float(np.interp(np.log(xh), lx, rs))
        blo = float(np.interp(np.log(xh), lx, los))
        bhi = float(np.interp(np.log(xh), lx, his))

        x0, ax = xi_model(u[ok], xi[ok])
        bs_ax = boot(lambda uu, yy: xi_model(uu, yy)[1], u[ok], xi[ok], rng)
        bs_x0 = boot(lambda uu, yy: xi_model(uu, yy)[0], u[ok], xi[ok], rng)
        se_ax = float(np.std(bs_ax, ddof=1))
        se_x0 = float(np.std(bs_x0, ddof=1))
        z = ax / se_ax if se_ax > 0 else float('nan')

        out['modes'][m] = {
            'n_ok': int(ok.sum()),
            'xi_hat': xh, 'bias_ratio': ratio, 'band_lo': blo, 'band_hi': bhi,
            'xi_debiased': xh / ratio,
            'xi_debiased_lo': xh / bhi, 'xi_debiased_hi': xh / blo,
            'xi0': x0, 'se_xi0_boot': se_x0,
            'alpha_xi': ax, 'se_alpha_xi_boot': se_ax, 'alpha_xi_z_boot': z,
        }
        print('%-5s %8.5f %9.4f [%.3f, %.3f] %10.6f [%.6f, %.6f] %11.6f %9.6f %7.2f'
              % (m, xh, ratio, blo, bhi, xh / ratio, xh / bhi, xh / blo,
                 ax, se_ax, z))

    out['script_sha256'] = hashlib.sha256(
        open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16]
    tmp = DEST + '.tmp'
    with open(tmp, 'wb') as fh:
        fh.write(json.dumps(out, indent=1).encode('utf-8'))
        fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, DEST)
    print('\nwrote %s  (%d windows)' % (DEST, len(rows)))

    print('\n--- LaTeX table body ---')
    for m in MODES:
        d = out['modes'][m]
        # xi_hat to five decimals: TX1 and TY1 differ by 1.7e-5 and printed
        # identically at four, which reads as a copy error rather than as the
        # coincidence it is.
        print('%s & %.5f & %.3f [%.3f, %.3f] & %.4f [%.4f, %.4f] & '
              '$%.5f \\pm %.5f$ & %.1f \\\\'
              % (m, d['xi_hat'], d['bias_ratio'], d['band_lo'], d['band_hi'],
                 d['xi_debiased'], d['xi_debiased_lo'], d['xi_debiased_hi'],
                 d['alpha_xi'], d['se_alpha_xi_boot'], d['alpha_xi_z_boot']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
