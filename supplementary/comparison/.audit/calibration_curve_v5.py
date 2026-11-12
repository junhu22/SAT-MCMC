"""Damping bias curve for configuration A, in the shape field_damping_table
consumes, with a bootstrap interval at each level.

Same generator and same recipe as damping_calibration_ab, restricted to the
configuration that was selected (fixed centre, 60 mHz half-width), and
extended to carry a 95 per cent interval per level so the field bias band can
be interpolated rather than asserted.

Six levels spanning 0.006 to 0.030 replace the four of the earlier curve,
which stopped at 0.012 and so did not cover R1's re-fitted damping.

Writes .audit/calibration_curve_v5.json.
Usage:  python .audit/calibration_curve_v5.py
"""
import hashlib
import json
import os
import sys
import time
import warnings

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'baselines'))
sys.path.insert(0, os.path.join(ROOT, '.audit'))
import recipe_v3 as RV                                    # noqa: E402
from bfft_gate import simulate                            # noqa: E402

F_REF = 0.447
HW = 0.060
LEVELS = (0.006, 0.010, 0.014, 0.018, 0.024, 0.030)
NREP = 96
NBOOT = 4000
SEED0 = 91000
BOOT_SEED = 20260826
SNR_DB = 60.0
DRIFT_MAX = 0.030


def main():
    n = int(1800.0 * RV.FS)
    rng = np.random.default_rng(BOOT_SEED)
    out = {'schema': 'calibration_curve_v5/1',
           'written_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'configuration': 'A: fixed centre at f_ref, half-width 60 mHz',
           'f_ref_hz': F_REF, 'half_width_hz': HW, 'n_rep': NREP,
           'n_boot': NBOOT, 'seed_base': SEED0, 'boot_seed': BOOT_SEED,
           'snr_db': SNR_DB, 'drift_max_hz': DRIFT_MAX,
           'supersedes': '.audit/calibration_curve.json, four levels to 0.012',
           'levels': {}}

    print('%-9s %8s %11s %11s %11s %9s'
          % ('xi_true', 'n_ok', 'median xi', 'ratio', 'ci_lo', 'ci_hi'))
    print('-' * 64)
    t0 = time.time()
    for li, xt in enumerate(LEVELS):
        drift = DRIFT_MAX * li / (len(LEVELS) - 1.0)
        f_true = F_REF - drift
        est = []
        for s in range(NREP):
            r = np.random.default_rng(SEED0 + 1000 * li + s)
            f, p, na = RV.spectrum(simulate(f_true, xt, n, RV.FS, SNR_DB, r))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                fit = RV.fit(f, p, F_REF, band_hw=HW, n_avg=na)
            if fit.get('ok') and np.isfinite(fit.get('xi', np.nan)):
                est.append(fit['xi'])
        est = np.array(est)
        med = float(np.median(est))
        bs = [np.median(est[rng.integers(0, len(est), len(est))])
              for _ in range(NBOOT)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        out['levels']['%.4f' % xt] = {
            'xi_true': xt, 'drift_mhz': drift * 1e3, 'n_ok': int(len(est)),
            'median_xi': med, 'ratio': med / xt,
            'ci_lo': float(lo / xt), 'ci_hi': float(hi / xt)}
        print('%-9.4f %8d %11.5f %11.4f %11.4f %9.4f'
              % (xt, len(est), med, med / xt, lo / xt, hi / xt))

    rs = [out['levels'][k]['ratio'] for k in sorted(out['levels'], key=float)]
    xs = [float(k) for k in sorted(out['levels'], key=float)]
    b = np.polyfit(xs, rs, 1)
    out['ratio_slope_per_unit_xi'] = float(b[0])
    out['ratio_range'] = [min(rs), max(rs)]
    print('\n  ratio spans %.3f to %.3f; slope on xi %+.3f'
          % (min(rs), max(rs), b[0]))
    print('  wall %.0f s' % (time.time() - t0))

    out['script_sha256'] = hashlib.sha256(
        open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16]
    dest = os.path.join(ROOT, '.audit', 'calibration_curve_v5.json')
    json.dump(out, open(dest, 'w', encoding='utf-8'), indent=1)
    print('wrote %s' % dest)
    return 0


if __name__ == '__main__':
    sys.exit(main())
