"""The finalised per-window recipe (batch dc64c2b, item 1).

Four elements, each traceable to the finding that forced it:

  nperseg 16384      "sufficient", not "full", resolution. dF = 1.95 mHz, n_avg = 6.
                     Full resolution (n_avg = 1) is HARMFUL: dc64c2b measured a
                     ~25% negative finite-sample bias in the MLE there, decaying
                     as 1/n_avg. Two biases oppose each other and both are
                     monotone in nperseg, so the choice is an interior optimum,
                     not a limit.
  Gamma(n_avg)       the averaged periodogram ordinate is Gamma(n) distributed;
                     least squares on it is neither the MLE nor unbiased.
  normalisation      divide by the band mean so the objective is O(1). This was
                     the root cause of seed-sticking (bc1f707): PSD ordinates
                     around 1e-17 sit below curve_fit's default tolerances.
  no peak picking    fit f_i freely from the nominal centre. Locating a peak
                     first quantises f onto the bin grid and biases the band.

Exported so the v3 rerun and every diagnostic use ONE fitting path; a recipe
that differs between the calibration and the application is not a recipe.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.signal import get_window, welch

FS = 32.0
NPERSEG = 16384
BAND_HW = 0.020


def spectrum(x, fs=FS, nperseg=NPERSEG):
    """Welch PSD and the number of averaged segments actually realised."""
    if nperseg >= len(x):
        w = np.hanning(len(x))
        Y = np.fft.rfft(x * w)
        f = np.fft.rfftfreq(len(x), 1.0 / fs)
        return f, (np.abs(Y) ** 2) * (2.0 / (fs * np.sum(w ** 2))), 1
    f, p = welch(x, fs=fs, window=get_window('hann', nperseg),
                 nperseg=nperseg, noverlap=nperseg // 2, scaling='density')
    return f, p, (len(x) - nperseg // 2) // (nperseg // 2)


def _nll(q, fb, Xs, n, xi_fix=None):
    A, fi, C = np.exp(q[0]), q[1], np.exp(q[2])
    xi = xi_fix if xi_fix is not None else np.exp(q[3])
    r = fb / fi
    m = A / ((1.0 - r ** 2) ** 2 + (2.0 * xi * r) ** 2) + C
    return float(np.sum(n * np.log(m) + n * Xs / m))


def fit(f, p, f_center, band_hw=BAND_HW, n_avg=None, xi_fix=None,
        starts=(0.004, 0.010, 0.020), hessian=False):
    """Gamma(n) MLE for a single SDOF plus a flat floor. No peak picking."""
    n = 1 if n_avg is None else n_avg
    m = (f >= f_center - band_hw) & (f <= f_center + band_hw)
    if m.sum() < 8:
        return {'ok': False, 'reason': 'band too narrow'}
    fb, X = f[m], p[m]
    sc = float(np.mean(X))
    if not np.isfinite(sc) or sc <= 0:
        return {'ok': False, 'reason': 'degenerate band'}
    Xs = X / sc
    best = None
    for s in starts:
        q0 = [np.log(max(Xs.max() * 4 * s ** 2, 1e-12)), f_center,
              np.log(max(np.median(Xs), 1e-12))]
        if xi_fix is None:
            q0.append(np.log(s))
        r = minimize(_nll, q0, args=(fb, Xs, n, xi_fix), method='Nelder-Mead',
                     options={'maxiter': 20000, 'xatol': 1e-11, 'fatol': 1e-13})
        if best is None or r.fun < best.fun:
            best = r
    out = {'ok': bool(best.success or np.isfinite(best.fun)),
           'f': float(best.x[1]), 'nll': float(best.fun),
           'n_avg': int(n), 'n_bins': int(m.sum()), 'band_hw': float(band_hw),
           'xi': float(xi_fix) if xi_fix is not None else float(np.exp(best.x[3]))}
    if hessian and xi_fix is None:
        out['sigma_xi'] = _sigma_xi(best.x, fb, Xs, n)
    return out


def _sigma_xi(q, fb, Xs, n):
    """Laplace posterior SD for xi under a flat prior, via a numerical Hessian.

    The fit is in log xi, so the SD is transformed back: sd(xi) = xi * sd(log xi).
    """
    k, h = len(q), 1e-4
    H = np.zeros((k, k))
    f0 = _nll(q, fb, Xs, n)
    for i in range(k):
        for j in range(i, k):
            qa, qb, qc, qd = q.copy(), q.copy(), q.copy(), q.copy()
            qa[i] += h; qa[j] += h
            qb[i] += h; qb[j] -= h
            qc[i] -= h; qc[j] += h
            qd[i] -= h; qd[j] -= h
            H[i, j] = H[j, i] = ((_nll(qa, fb, Xs, n) - _nll(qb, fb, Xs, n)
                                 - _nll(qc, fb, Xs, n) + _nll(qd, fb, Xs, n))
                                / (4 * h * h)) if i != j else \
                (_nll(qa, fb, Xs, n) - 2 * f0 + _nll(qd, fb, Xs, n)) / (4 * h * h)
    try:
        cov = np.linalg.inv(H)
        v = float(cov[3, 3])
        if not np.isfinite(v) or v <= 0:
            return float('nan')
        return float(np.exp(q[3]) * np.sqrt(v))
    except np.linalg.LinAlgError:
        return float('nan')


def fit_adaptive(f, p, f_center, k_hpbw=5.0, n_avg=None, seed_band=BAND_HW,
                 hessian=False):
    """One iteration of band adaptation: fit fixed, then refit on +-k*HPBW.

    A FIXED band is not xi-neutral. Half-power bandwidth is 2*xi*f, so at
    f = 0.25 Hz a +-20 mHz band spans +-10 HPBW at xi = 0.004 but only +-2.7
    HPBW at xi = 0.015 -- the highest damping gets the narrowest relative band
    and loses the most skirt. If that truncation drives the residual slope in
    the 1.02 / 0.95 / 0.89 ratios, adapting the band should flatten it.
    """
    r0 = fit(f, p, f_center, band_hw=seed_band, n_avg=n_avg)
    if not r0['ok']:
        return r0
    hw = float(np.clip(k_hpbw * 2.0 * r0['xi'] * r0['f'], 0.004, 0.060))
    r = fit(f, p, r0['f'], band_hw=hw, n_avg=n_avg, hessian=hessian,
            starts=(r0['xi'] * 0.6, r0['xi'], r0['xi'] * 1.6))
    r['xi_seed'] = r0['xi']
    r['k_hpbw'] = float(k_hpbw)
    return r
