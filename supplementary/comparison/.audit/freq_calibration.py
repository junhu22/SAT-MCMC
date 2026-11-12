"""Item 1: is the frequency bias amplitude-dependent? Pre-reg freq_calibration_prereg.md."""
import csv, hashlib, json, os, sys, time, warnings
import numpy as np
from scipy.optimize import curve_fit
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'baselines'))
sys.path.insert(0, os.path.join(ROOT, '.audit'))
from bfft_gate import simulate
import recipe_v3 as RV
NS_A, BASE_A = 96, 61000
NS_B, BASE_B = 96, 62000
F0, XI = 0.25, 0.007
A_REF = 100.0
ALPHA, GAMMA, F0B = 0.020, 0.40, 0.247


def field_amps():
    rows = [r for r in csv.DictReader(open(os.path.join(ROOT, 'results/per_window_modes_v3.csv'),
                                           encoding='utf-8')) if r['is_outlier'] != 'True']
    return np.array([float(r['amplitude_ug']) for r in rows])


def snr_of(a, amin):
    return 30.0 + 20.0 * np.log10(max(a / amin, 1.0))


def stage_a(amps):
    lv = np.percentile(amps, [5, 20, 40, 60, 80, 97])
    amin = amps.min()
    n = int(1800.0 * RV.FS)
    rows = []
    for a in lv:
        est = []
        for s in range(NS_A):
            rng = np.random.default_rng(BASE_A + s)
            f, p, na = RV.spectrum(simulate(F0, XI, n, RV.FS, snr_of(a, amin), rng))
            r = RV.fit(f, p, F0, n_avg=na)
            if r['ok'] and np.isfinite(r['f']):
                est.append(r['f'])
        e = np.array(est)
        rows.append({'amplitude_ug': float(a), 'snr_db': float(snr_of(a, amin)),
                     'n_ok': int(len(e)), 'median_f': float(np.median(e)),
                     'ratio_f': float(np.median(e) / F0),
                     'bias_mhz': float((np.median(e) - F0) * 1e3),
                     'sd_mhz': float(np.std(e, ddof=1) * 1e3)})
    return rows


def powerlaw(a, f0, al, ga):
    return f0 * (1.0 - al * (a / A_REF) ** ga)


def stage_b(amps):
    n = int(1800.0 * RV.FS)
    f_true = powerlaw(amps, F0B, ALPHA, GAMMA)
    amin = amps.min()
    als, gas, f0s = [], [], []
    for rep in range(NS_B):
        fh = []
        for i, (a, ft) in enumerate(zip(amps, f_true)):
            rng = np.random.default_rng(BASE_B + rep * 1000 + i)
            f, p, na = RV.spectrum(simulate(ft, XI, n, RV.FS, snr_of(a, amin), rng))
            r = RV.fit(f, p, F0B, n_avg=na)
            fh.append(r['f'] if (r['ok'] and np.isfinite(r['f'])) else np.nan)
        fh = np.array(fh); ok = np.isfinite(fh)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                p_, _ = curve_fit(powerlaw, amps[ok], fh[ok],
                                  p0=[F0B, ALPHA, GAMMA],
                                  bounds=([0.2, 0.001, 0.05], [0.3, 0.5, 2.0]),
                                  maxfev=40000)
            f0s.append(p_[0]); als.append(p_[1]); gas.append(p_[2])
        except Exception:
            pass
    return np.array(f0s), np.array(als), np.array(gas)


def band(v, truth, rb):
    bs = np.array([np.median(rb.choice(v, len(v), replace=True)) for _ in range(4000)]) / truth
    return {'n': int(len(v)), 'median': float(np.median(v)),
            'ratio': float(np.median(v) / truth),
            'ci_lo': float(np.percentile(bs, 2.5)), 'ci_hi': float(np.percentile(bs, 97.5))}


def main():
    rb = np.random.default_rng(2718)
    amps = field_amps()
    out = {'schema': 'freq_calibration/1', 'written_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'prereg': '.audit/freq_calibration_prereg.md',
           'n_windows': int(len(amps)), 'amp_range_ug': [float(amps.min()), float(amps.max())]}
    print('stage A: does the frequency bias vary with amplitude? truth has NO a-dependence')
    print('%d field windows, amplitude %.1f-%.1f ug, %d seeds per level\n'
          % (len(amps), amps.min(), amps.max(), NS_A))
    ra = stage_a(amps); out['stage_a'] = ra
    print('%-12s %8s %6s %12s %11s %9s' %
          ('amp (ug)', 'SNR dB', 'n_ok', 'median f', 'bias (mHz)', 'sd (mHz)'))
    print('-' * 64)
    for r in ra:
        print('%-12.1f %8.1f %6d %12.6f %11.4f %9.4f'
              % (r['amplitude_ug'], r['snr_db'], r['n_ok'], r['median_f'],
                 r['bias_mhz'], r['sd_mhz']))
    la = np.log(np.array([r['amplitude_ug'] for r in ra]))
    bi = np.array([r['bias_mhz'] for r in ra])
    sl = float(np.polyfit(la, bi, 1)[0])
    sd = np.array([r['sd_mhz'] for r in ra]) / np.sqrt(NS_A)
    se_sl = float(np.sqrt(np.sum(sd ** 2)) / (la.max() - la.min()))
    out['stage_a_slope_mhz_per_ln_a'] = sl
    out['stage_a_slope_se'] = se_sl
    out['stage_a_flat'] = bool(abs(sl) < 2 * se_sl)
    print('\nbias slope %+.4f mHz per ln(a), SE %.4f -> %s'
          % (sl, se_sl, 'FLAT (no amplitude dependence in the bias)'
             if out['stage_a_flat'] else 'SLOPED -- the estimator manufactures a-dependence'))

    print('\nstage B: alpha/gamma bias bands. truth alpha=%.3f gamma=%.2f f0=%.3f (v12 table)'
          % (ALPHA, GAMMA, F0B))
    f0s, als, gas = stage_b(amps)
    if len(als) < 10:
        print('  too few successful fits (%d) -- stage B inconclusive' % len(als))
        out['stage_b'] = {'n': int(len(als)), 'status': 'inconclusive'}
    else:
        out['stage_b'] = {'f0': band(f0s, F0B, rb), 'alpha': band(als, ALPHA, rb),
                          'gamma': band(gas, GAMMA, rb)}
        print('\n%-8s %10s %10s %10s %18s' % ('param', 'truth', 'median', 'ratio', '95% CI on ratio'))
        print('-' * 62)
        for k, t in (('f0', F0B), ('alpha', ALPHA), ('gamma', GAMMA)):
            b = out['stage_b'][k]
            print('%-8s %10.4f %10.4f %10.3f     [%.3f, %.3f]'
                  % (k, t, b['median'], b['ratio'], b['ci_lo'], b['ci_hi']))
    out['script_sha256'] = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16]
    dest = os.path.join(ROOT, '.audit', 'freq_calibration.json')
    tmp = dest + '.tmp'
    with open(tmp, 'wb') as fh:
        fh.write(json.dumps(out, indent=1).encode('utf-8')); fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, dest)
    print('wrote %s' % dest)
    return 0


if __name__ == '__main__':
    warnings.simplefilter('ignore')
    sys.exit(main())
