"""Frequency Domain Decomposition, baseline M0 of the head-to-head study.

FDD (Brincker, Zhang and Andersen 2000): form the output PSD matrix G(f) from
the measured channels, take its SVD at every frequency line, and read modal
frequencies off the peaks of the first singular value. The first singular
vector at a peak estimates the mode shape.

What this implementation deliberately does NOT do, per the frozen
preregistration in .audit/h2h_preregistration.md section 3:

  FDD produces NO damping estimate. The half-power method could be bolted on,
  but at the 42 and 33 mHz separations of the TX2 cluster it would be reading
  the neighbouring peak's skirt, and the number it returned would be a fitting
  artefact rather than a damping ratio. M0 therefore reports frequencies only,
  the damping column is n/a, and every comparison sentence involving M0 is
  restricted to frequency estimates.

Closely-spaced handling, also from the preregistration: peaks are picked from
the first singular value; where two modes merge into one peak the window is
recorded as unresolved for that pair and is NOT interpolated or split.
Identity is assigned by ascending frequency against a low-amplitude anchor,
and a window where adjacent modes cross is flagged identity_ambiguous and
excluded from P1/P2 rather than silently relabelled.
"""
import numpy as np
from scipy.signal import csd, get_window


def psd_matrix(x, fs, nperseg, noverlap=None):
    """x: (n_samples, n_ch) -> f, G with G shape (n_freq, n_ch, n_ch)."""
    n_ch = x.shape[1]
    noverlap = nperseg // 2 if noverlap is None else noverlap
    win = get_window('hann', nperseg)
    f = None
    G = None
    for i in range(n_ch):
        for j in range(i, n_ch):
            f, p = csd(x[:, i], x[:, j], fs=fs, window=win,
                       nperseg=nperseg, noverlap=noverlap, scaling='density')
            if G is None:
                G = np.zeros((len(f), n_ch, n_ch), complex)
            G[:, i, j] = p
            if j != i:
                G[:, j, i] = np.conj(p)
    return f, G


def fdd_singular(f, G):
    """SVD at every line. Returns s (n_freq, n_ch) and U (n_freq, n_ch, n_ch)."""
    n_f, n_ch, _ = G.shape
    s = np.empty((n_f, n_ch))
    U = np.empty((n_f, n_ch, n_ch), complex)
    for k in range(n_f):
        u, sv, _ = np.linalg.svd(G[k])
        s[k] = sv
        U[k] = u
    return s, U


def pick_peaks(f, s1, bands, min_prominence_db=3.0):
    """One peak per band, from the first singular value.

    A band whose maximum is not at least `min_prominence_db` above the band
    edges is reported unresolved rather than given a frequency.
    """
    out = []
    logs = 10.0 * np.log10(np.maximum(s1, 1e-300))
    for (lo, hi) in bands:
        m = (f >= lo) & (f <= hi)
        if m.sum() < 5:
            out.append({'f': None, 'resolved': False, 'reason': 'band too narrow'})
            continue
        idx = np.where(m)[0]
        k = idx[np.argmax(s1[idx])]
        edge = max(logs[idx[0]], logs[idx[-1]])
        prom = logs[k] - edge
        if k in (idx[0], idx[-1]) or prom < min_prominence_db:
            out.append({'f': float(f[k]), 'resolved': False,
                        'reason': 'no interior peak with %.1f dB prominence' % min_prominence_db,
                        'prominence_db': float(prom)})
            continue
        # parabolic refinement on the log spectrum
        y0, y1, y2 = logs[k - 1], logs[k], logs[k + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        df = f[1] - f[0]
        out.append({'f': float(f[k] + delta * df), 'resolved': True,
                    'prominence_db': float(prom), 'bin': int(k)})
    return out


def mac(a, b):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()
    num = abs(np.vdot(a, b)) ** 2
    return float(num / (np.vdot(a, a).real * np.vdot(b, b).real))


def run_fdd(x, fs, bands, nperseg=4096):
    f, G = psd_matrix(x, fs, nperseg)
    s, U = fdd_singular(f, G)
    peaks = pick_peaks(f, s[:, 0], bands)
    shapes = []
    for p in peaks:
        if not p['resolved']:
            shapes.append(None)
            continue
        v = U[p['bin'], :, 0]
        v = v / np.linalg.norm(v)
        ref = np.argmax(np.abs(v))
        v = v * np.exp(-1j * np.angle(v[ref]))
        shapes.append(np.real(v))
    return {'f_axis': f, 's1': s[:, 0], 'peaks': peaks, 'shapes': shapes,
            'damping': [None] * len(peaks)}
