"""Extend the damping bias curve below its floor: xi = 0.004 and 0.005.

The reported curve starts at 0.006 while the resolution floor is 0.004, so
between those two there is a band of damping values the field data reach and
the calibration does not cover. This closes that gap, or shows that it cannot
be closed.

Everything is held to the existing curve: the same generator, the same recipe,
the same 60 mHz half-width about the same 0.447 Hz reference, 96 realisations,
the same bootstrap for the interval.

One thing needs saying about drift. The six-level curve assigns a band-centre
offset by level *index*, 0 mHz at the lowest level rising to 30 mHz at the
highest, so re-running the whole set with two extra levels prepended would give
every existing level a different offset and change six archived ratios. The two
new levels are therefore run at 0 mHz, which is what the level adjacent to them
carries and what the pattern implies for the lowest damping ratios. The
schedule is recorded per level in the output so the two sets can be told apart.

Writes .audit/calibration_curve_low.json.
Usage:  python .audit/calibration_curve_low.py
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
NEW_LEVELS = (0.004, 0.005)
DRIFT_MHZ = 0.0
NREP = 96
NBOOT = 4000
SEED0 = 91000            # same base; the level offset keeps the streams apart
BOOT_SEED = 20260826
SNR_DB = 60.0


def main():
    n = int(1800.0 * RV.FS)
    rng = np.random.default_rng(BOOT_SEED)
    out = {'schema': 'calibration_curve_low/1',
           'written_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'extends': '.audit/calibration_curve_v5.json',
           'configuration': 'A: fixed centre at f_ref, half-width 60 mHz',
           'f_ref_hz': F_REF, 'half_width_hz': HW, 'n_rep': NREP,
           'n_boot': NBOOT, 'boot_seed': BOOT_SEED, 'snr_db': SNR_DB,
           'drift_mhz': DRIFT_MHZ,
           'drift_note': ('run at the offset the lowest existing level carries; '
                          'the six-level schedule assigns offset by index, so '
                          'prepending levels there would have moved six '
                          'archived ratios'),
           'levels': {}}

    print('%-9s %8s %11s %11s %11s %9s'
          % ('xi_true', 'n_ok', 'median xi', 'ratio', 'ci_lo', 'ci_hi'))
    print('-' * 64)
    t0 = time.time()
    for li, xt in enumerate(NEW_LEVELS):
        f_true = F_REF - DRIFT_MHZ / 1e3
        est = []
        for s in range(NREP):
            # level index continues below zero so no stream is shared with the
            # six-level run
            r = np.random.default_rng(SEED0 - 1000 * (li + 1) + s)
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
            'xi_true': xt, 'drift_mhz': DRIFT_MHZ, 'n_ok': int(len(est)),
            'median_xi': med, 'ratio': med / xt,
            'ci_lo': float(lo / xt), 'ci_hi': float(hi / xt)}
        print('%-9.4f %8d %11.5f %11.4f %11.4f %9.4f'
              % (xt, len(est), med, med / xt, lo / xt, hi / xt))

    out['script_sha256'] = hashlib.sha256(
        open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16]
    dest = os.path.join(ROOT, '.audit', 'calibration_curve_low.json')
    json.dump(out, open(dest, 'w', encoding='utf-8'), indent=1)
    print('  wall %.0f s' % (time.time() - t0))
    print('wrote %s' % dest)
    return 0


if __name__ == '__main__':
    sys.exit(main())
