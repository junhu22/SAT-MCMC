"""Item 2: bias band for alpha_xi.

The calibration curve was measured at FIXED xi and therefore calibrates xi_hat
only. alpha_xi is fitted from the SHAPE of xi_hat(a), so it needs its own
calibration -- and the generator has to match the model being fitted, which
means BOTH f and xi vary with amplitude, exactly as acceptance rule 2 requires.

Pre-registration (this docstring, committed before the run):

  truth, fixed before running:
    f(a)  = f0 [1 - alpha_f (a/a_ref)^gamma],  f0 = 0.247, alpha_f = 0.020,
            gamma = 0.40                       (v12's own synthetic table, TX1)
    xi(a) = xi0 + alpha_xi ln(a/a_ref),        xi0 = 0.009, alpha_xi = 0.0012
            chosen near the v3 field estimates so the band applies AT the
            operating point, not somewhere else on the curve
    a_ref = 100 ug, amplitudes = the real 141-window distribution
    SNR   = 30 + 20 log10(a/a_min)

  seeds 71000+, 96 repetitions, disjoint from every previous block
  output: median and bootstrap 95% band on alpha_xi_hat/alpha_xi, and the same
          for xi0, alpha_f and gamma, since this generator exercises all four

  not a gate: no threshold, no pass/fail. If alpha_f and gamma come back biased
  here, that is a SECOND route to the same question stage B asks, and both are
  reported whichever way they point.
"""
import csv
import hashlib
import json
import os
import sys
import time
import warnings

import numpy as np
from scipy.optimize import curve_fit

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'baselines'))
sys.path.insert(0, os.path.join(ROOT, '.audit'))
from bfft_gate import simulate
import recipe_v3 as RV

NREP, BASE = 96, 71000
A_REF = 100.0
F0, ALPHA_F, GAMMA = 0.247, 0.020, 0.40
XI0, ALPHA_XI = 0.009, 0.0012


def field_amps():
    rows = [r for r in csv.DictReader(
        open(os.path.join(ROOT, 'results/per_window_modes_v3.csv'), encoding='utf-8'))
        if r['is_outlier'] != 'True']
    return np.array([float(r['amplitude_ug']) for r in rows])


def f_of(a):
    return F0 * (1.0 - ALPHA_F * (a / A_REF) ** GAMMA)


def xi_of(a):
    return XI0 + ALPHA_XI * np.log(a / A_REF)


def powerlaw(a, f0, al, ga):
    return f0 * (1.0 - al * (a / A_REF) ** ga)


def loglin(a, x0, ax):
    return x0 + ax * np.log(a / A_REF)


def band(v, truth, rb):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    bs = np.array([np.median(rb.choice(v, len(v), replace=True)) for _ in range(4000)]) / truth
    return {'n': int(len(v)), 'median': float(np.median(v)),
            'ratio': float(np.median(v) / truth),
            'ci_lo': float(np.percentile(bs, 2.5)),
            'ci_hi': float(np.percentile(bs, 97.5))}


def main():
    rb = np.random.default_rng(1618)
    amps = field_amps()
    amin = amps.min()
    n = int(1800.0 * RV.FS)
    ft, xt = f_of(amps), xi_of(amps)
    out = {'schema': 'alpha_xi_calibration/1', 'written_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'not_a_gate': True, 'n_windows': int(len(amps)), 'n_rep': NREP,
           'seeds': '%d+' % BASE,
           'truth': {'f0': F0, 'alpha_f': ALPHA_F, 'gamma': GAMMA,
                     'xi0': XI0, 'alpha_xi': ALPHA_XI, 'a_ref_ug': A_REF}}
    print('alpha_xi calibration: generator varies BOTH f and xi with amplitude')
    print('%d windows, %.1f-%.1f ug, %d repetitions, seeds %d+' %
          (len(amps), amps.min(), amps.max(), NREP, BASE))
    print('truth: alpha_f %.3f gamma %.2f | xi0 %.4f alpha_xi %.5f' %
          (ALPHA_F, GAMMA, XI0, ALPHA_XI))
    print('true xi spans %.5f to %.5f over the amplitude range\n' % (xt.min(), xt.max()))

    f0s, als, gas, x0s, axs = [], [], [], [], []
    t0 = time.time()
    for rep in range(NREP):
        fh, xh = [], []
        for i, (a, f_i, x_i) in enumerate(zip(amps, ft, xt)):
            rng = np.random.default_rng(BASE + rep * 1000 + i)
            snr = 30.0 + 20.0 * np.log10(max(a / amin, 1.0))
            f, p, na = RV.spectrum(simulate(f_i, max(x_i, 5e-4), n, RV.FS, snr, rng))
            r = RV.fit(f, p, F0, n_avg=na)
            ok = r['ok'] and np.isfinite(r['f']) and np.isfinite(r['xi'])
            fh.append(r['f'] if ok else np.nan)
            xh.append(r['xi'] if ok else np.nan)
        fh, xh = np.array(fh), np.array(xh)
        m = np.isfinite(fh) & np.isfinite(xh) & (xh > 0)
        if m.sum() < 30:
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pf, _ = curve_fit(powerlaw, amps[m], fh[m], p0=[F0, ALPHA_F, GAMMA],
                                  bounds=([0.2, 0.001, 0.05], [0.3, 0.5, 2.0]), maxfev=60000)
                px, _ = curve_fit(loglin, amps[m], xh[m], p0=[XI0, ALPHA_XI], maxfev=60000)
            f0s.append(pf[0]); als.append(pf[1]); gas.append(pf[2])
            x0s.append(px[0]); axs.append(px[1])
        except Exception:
            pass
    wall = time.time() - t0

    print('%d of %d repetitions fitted (%.0f s)\n' % (len(axs), NREP, wall))
    print('%-10s %11s %11s %9s %20s' % ('param', 'truth', 'median', 'ratio', '95% CI on ratio'))
    print('-' * 68)
    out['bands'] = {}
    for name, v, t in (('f0', f0s, F0), ('alpha_f', als, ALPHA_F), ('gamma', gas, GAMMA),
                       ('xi0', x0s, XI0), ('alpha_xi', axs, ALPHA_XI)):
        if len(v) < 10:
            out['bands'][name] = {'status': 'too few fits', 'n': int(len(v))}
            print('%-10s  too few fits (%d)' % (name, len(v)))
            continue
        b = band(v, t, rb)
        out['bands'][name] = b
        print('%-10s %11.5f %11.5f %9.3f     [%.3f, %.3f]'
              % (name, t, b['median'], b['ratio'], b['ci_lo'], b['ci_hi']))

    ax = out['bands'].get('alpha_xi', {})
    if 'ratio' in ax:
        print('\nalpha_xi bias band: %.3f [%.3f, %.3f]' % (ax['ratio'], ax['ci_lo'], ax['ci_hi']))
        print('sign of alpha_xi recovered: %s'
              % ('yes' if ax['median'] * ALPHA_XI > 0 else 'NO -- sign reversed'))
        out['alpha_xi_sign_recovered'] = bool(ax['median'] * ALPHA_XI > 0)
    al = out['bands'].get('alpha_f', {})
    if 'ratio' in al:
        print('\nSECOND ROUTE to the alpha question: this generator also carries a known')
        print('alpha_f, so alpha_f ratio %.3f [%.3f, %.3f] speaks to whether the'
              % (al['ratio'], al['ci_lo'], al['ci_hi']))
        print('estimator inflates alpha -- independently of stage B.')
        out['alpha_f_inflated'] = bool(al['ci_lo'] > 1.10)

    out['script_sha256'] = hashlib.sha256(
        open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16]
    dest = os.path.join(ROOT, '.audit', 'alpha_xi_calibration.json')
    tmp = dest + '.tmp'
    with open(tmp, 'wb') as fh:
        fh.write(json.dumps(out, indent=1).encode('utf-8'))
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, dest)
    print('\nwrote %s' % dest)
    return 0


if __name__ == '__main__':
    warnings.simplefilter('ignore')
    sys.exit(main())
