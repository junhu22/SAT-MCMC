"""Synthetic identifiability gate for damping at this record length.

The question is binary and the answer decides whether the damping line
continues at all: at 1800 s, 32 Hz, and the field's own noise level, is a
damping ratio of order 1% recoverable?

Estimator. Single-mode Bayesian FFT most-probable-value: fit

    E[|Y(f_k)|^2] = S |H(f_k; f_0, xi)|^2 + S_e

to the FULL-RESOLUTION periodogram over the resonance band by maximising the
exponential (Whittle) likelihood, which for a single averaged periodogram
ordinate is  -sum_k [ y_k/mu_k + log mu_k ]. No Welch averaging, so no window
broadening enters -- which is the whole point, since window broadening is what
corrupts the current pipeline. This is the single-mode core of the Bayesian FFT
method; it is NOT the full Li and Au EM for multiple modes and setups, and the
docstring says so rather than implying more than was built.

Uncertainty. sigma_post(xi) from the observed information at the MPV (second
difference of the negative log-likelihood). kappa_xi = sigma_post/sigma_prior
with a uniform prior on [1e-4, 0.05], sigma_prior = 0.05/sqrt(12).

Gate, both conditions required:
    kappa_xi <= 0.316          the posterior is data-dominated
    median |xi_hat/xi - 1| <= 0.20   and it is in the right place

Failing this at 1% is not a null result. It converts 'the damping numbers are
unreliable' into 'damping at this level is not identifiable from a record of
this length at this SNR', which is a stronger and cleaner statement.
"""
import hashlib
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
FS, DUR = 32.0, 1800.0
F0 = 0.25
XIS = [0.005, 0.0075, 0.01, 0.015, 0.02]
N_SEED = 20
BAND_HW = 0.03
XI_PRIOR_LO, XI_PRIOR_HI = 1e-4, 0.05
SIGMA_PRIOR = (XI_PRIOR_HI - XI_PRIOR_LO) / np.sqrt(12.0)


def simulate(f0, xi, n, fs, snr_db, rng):
    """White-noise-driven SDOF, discretised, plus measurement noise."""
    w = 2 * np.pi * f0
    dt = 1.0 / fs
    wd = w * np.sqrt(max(1 - xi * xi, 1e-12))
    a1 = 2 * np.exp(-xi * w * dt) * np.cos(wd * dt)
    a2 = -np.exp(-2 * xi * w * dt)
    e = rng.standard_normal(n)
    q = np.zeros(n)
    for t in range(2, n):
        q[t] = a1 * q[t - 1] + a2 * q[t - 2] + e[t]
    q /= q.std()
    noise = rng.standard_normal(n) * 10 ** (-snr_db / 20.0)
    return q + noise


def periodogram(x, fs):
    n = len(x)
    Y = np.fft.rfft(x * np.hanning(n))
    f = np.fft.rfftfreq(n, 1.0 / fs)
    scale = 2.0 / (fs * np.sum(np.hanning(n) ** 2))
    return f, (np.abs(Y) ** 2) * scale


def nll(theta, f, y):
    S, f0, xi, Se = np.exp(theta[0]), theta[1], np.exp(theta[2]), np.exp(theta[3])
    r = f / f0
    H2 = 1.0 / ((1 - r ** 2) ** 2 + (2 * xi * r) ** 2 + 1e-300)
    mu = S * H2 + Se
    return float(np.sum(y / mu + np.log(mu)))


def fit_mpv(f, y, f0_guess):
    from scipy.optimize import minimize
    m = np.abs(f - f0_guess) <= BAND_HW
    fb, yb = f[m], y[m]
    k = int(np.argmax(yb))
    th0 = np.array([np.log(max(yb[k] * 4e-4, 1e-12)), fb[k], np.log(0.01),
                    np.log(max(np.median(yb), 1e-12))])
    res = minimize(nll, th0, args=(fb, yb), method='Nelder-Mead',
                   options={'maxiter': 20000, 'xatol': 1e-10, 'fatol': 1e-10})
    th = res.x
    xi_hat = float(np.exp(th[2]))
    # observed information in log xi, converted to sigma on xi
    h = 0.02
    lp, l0, lm = (nll(th + np.array([0, 0, h, 0]), fb, yb), nll(th, fb, yb),
                  nll(th - np.array([0, 0, h, 0]), fb, yb))
    curv = (lp - 2 * l0 + lm) / h ** 2
    sig_logxi = float(1.0 / np.sqrt(curv)) if curv > 0 else np.inf
    return xi_hat, float(th[1]), xi_hat * sig_logxi, res.success


def field_noise_floor_db():
    """SNR proxy from the field record: in-band peak over the 8-14 Hz floor."""
    try:
        import h5py
        with h5py.File(os.path.join(ROOT, 'processed/mangkhut_processed.h5'), 'r') as h5:
            x = np.asarray(h5['window_0053']['acc_filtered'])[:, 0]
        f, p = periodogram(x, FS)
        band = (f > 0.20) & (f < 0.30)
        hi = (f > 8.0) & (f < 14.0)
        with open(os.path.join(ROOT, '.audit', 'quarantine_access.log'), 'ab') as fh:
            fh.write(('%s\tbfft_gate.py\tprocessed/mangkhut_processed.h5\t'
                      'noise floor for the synthetic SNR\n'
                      % time.strftime('%Y-%m-%dT%H:%M:%S')).encode('utf-8'))
        return float(10 * np.log10(p[band].max() / np.median(p[hi]))), True
    except Exception:
        return 30.0, False


def main():
    snr_db, measured = field_noise_floor_db()
    print('field in-band peak over 8-14 Hz floor: %.1f dB %s'
          % (snr_db, '(measured)' if measured else '(fallback)'))
    n = int(DUR * FS)
    out = {'schema': 'bfft_damping_gate/1', 'written_at': time.strftime('%Y-%m-%dT%H:%M:%S'),
           'fs': FS, 'duration_s': DUR, 'f0': F0, 'n_seeds': N_SEED,
           'snr_db_from_field': snr_db, 'snr_measured': measured,
           'sigma_prior_xi': SIGMA_PRIOR, 'gate': {'kappa_max': 0.316, 'bias_max': 0.20},
           'levels': {}}
    print('\n%-8s %9s %11s %11s %10s %9s %s'
          % ('xi_true', 'n_ok', 'med xi_hat', 'med bias', 'med kappa', 'frac_ok', 'GATE'))
    print('-' * 78)
    xi_min = None
    for xi in XIS:
        hats, kaps = [], []
        for s in range(N_SEED):
            rng = np.random.default_rng(1000 * s + int(xi * 1e5))
            x = simulate(F0, xi, n, FS, snr_db, rng)
            f, y = periodogram(x, FS)
            xh, fh, sig, ok = fit_mpv(f, y, F0)
            if ok and np.isfinite(xh) and np.isfinite(sig):
                hats.append(xh)
                kaps.append(sig / SIGMA_PRIOR)
        hats, kaps = np.array(hats), np.array(kaps)
        if len(hats) == 0:
            print('%-8.4f %9d  all fits failed' % (xi, 0))
            continue
        bias = np.median(np.abs(hats / xi - 1.0))
        kap = float(np.median(kaps))
        passed = bool(kap <= 0.316 and bias <= 0.20)
        if passed and xi_min is None:
            xi_min = xi
        out['levels'][str(xi)] = {'n_ok': int(len(hats)), 'median_xi_hat': float(np.median(hats)),
                                  'median_abs_bias': float(bias), 'median_kappa': kap,
                                  'frac_within_20pct': float(np.mean(np.abs(hats / xi - 1) <= 0.20)),
                                  'pass': passed}
        print('%-8.4f %9d %11.5f %10.1f%% %10.4f %8.2f  %s'
              % (xi, len(hats), np.median(hats), 100 * bias, kap,
                 np.mean(np.abs(hats / xi - 1) <= 0.20), 'PASS' if passed else 'FAIL'))
    out['xi_min_identifiable'] = xi_min
    print('\nlowest identifiable damping at this record length and SNR: %s'
          % (('%.4f' % xi_min) if xi_min else 'NONE of the tested levels'))
    out['verdict'] = ('identifiable from xi >= %.4f' % xi_min) if xi_min else \
        'NOT IDENTIFIABLE at any tested level, including 2%'
    out['config_hash'] = hashlib.sha256(json.dumps(
        {'fs': FS, 'dur': DUR, 'f0': F0, 'xis': XIS, 'seeds': N_SEED,
         'snr': round(snr_db, 3)}, sort_keys=True).encode()).hexdigest()[:16]
    out['script_sha256'] = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()[:16]
    dest = os.path.join(ROOT, '.audit', 'bfft_damping_gate.json')
    tmp = dest + '.tmp'
    with open(tmp, 'wb') as fh:
        fh.write(json.dumps(out, indent=1).encode('utf-8'))
        fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, dest)
    print('wrote %s  config_hash=%s' % (dest, out['config_hash']))
    return 0 if xi_min else 1


if __name__ == '__main__':
    sys.exit(main())
