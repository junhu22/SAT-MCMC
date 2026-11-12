"""Covariance-driven Stochastic Subspace Identification, baseline M1.

Standard SSI-COV: build the block Toeplitz matrix of output covariances, take
its SVD, split into observability and controllability factors, recover the
state matrix A and output matrix C, and read modal parameters off the
eigenvalues of A. Sweep the model order to build a stabilisation diagram and
keep only poles that persist across orders.

Unlike FDD, SSI-COV DOES estimate damping, so M1 participates in the damping
comparison. The preregistration requires the order range and the stability
criteria to be fixed before any run, and they are the module constants below.

Closely-spaced handling matches the FDD rule and the preregistration: a band
is assigned a mode only if a stable cluster is found inside it. If two nominal
modes fall in one cluster, or a band yields no stable cluster, the result is
UNRESOLVED for that band. Nothing is interpolated and no cluster is split to
make the count come out right.
"""
import numpy as np

# frozen before any run, per .audit/h2h_preregistration.md
ORDER_MIN, ORDER_MAX, ORDER_STEP = 10, 60, 2
N_BLOCK_ROWS = 20                 # block rows i in the Toeplitz matrix
TOL_FREQ = 0.01                   # 1% frequency change between adjacent orders
TOL_DAMP = 0.05                   # 5% damping change
TOL_MAC = 0.98                    # mode-shape MAC between adjacent orders
DAMP_MIN, DAMP_MAX = 1e-4, 0.20   # physically admissible damping window


def block_toeplitz(y, i):
    """y: (n_samples, n_ch) -> T1 of shape (i*n_ch, i*n_ch) from covariances."""
    n, nc = y.shape
    y = y - y.mean(0)
    maxlag = 2 * i
    R = np.empty((maxlag + 1, nc, nc))
    for k in range(maxlag + 1):
        R[k] = (y[k:].T @ y[:n - k]) / (n - k)
    T = np.empty((i * nc, i * nc))
    for a in range(i):
        for b in range(i):
            T[a * nc:(a + 1) * nc, b * nc:(b + 1) * nc] = R[i + a - b]
    return T, R


def modes_at_order(T, nc, i, order, fs):
    U, s, Vt = np.linalg.svd(T)
    U1 = U[:, :order]
    S1 = np.diag(np.sqrt(s[:order]))
    O = U1 @ S1                              # observability
    C = O[:nc, :]
    A = np.linalg.pinv(O[:-nc, :]) @ O[nc:, :]
    lam, V = np.linalg.eig(A)
    out = []
    dt = 1.0 / fs
    for k in range(len(lam)):
        if np.imag(lam[k]) <= 0:
            continue                          # keep one of each conjugate pair
        mu = np.log(lam[k]) / dt
        wn = abs(mu)
        if wn <= 0:
            continue
        f = wn / (2 * np.pi)
        z = -np.real(mu) / wn
        if not (DAMP_MIN < z < DAMP_MAX) or not np.isfinite(f):
            continue
        phi = C @ V[:, k]
        out.append({'f': float(f), 'z': float(z), 'phi': phi})
    return out


def mac_c(a, b):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()
    return float(abs(np.vdot(a, b)) ** 2 / (np.vdot(a, a).real * np.vdot(b, b).real))


def stabilisation(y, fs, orders=None, i=N_BLOCK_ROWS):
    nc = y.shape[1]
    T, _ = block_toeplitz(y, i)
    orders = orders or list(range(ORDER_MIN, ORDER_MAX + 1, ORDER_STEP))
    per_order = {o: modes_at_order(T, nc, i, o, fs) for o in orders}

    stable = []
    for a, b in zip(orders[:-1], orders[1:]):
        for m in per_order[b]:
            for p in per_order[a]:
                if (abs(m['f'] - p['f']) / max(p['f'], 1e-12) < TOL_FREQ
                        and abs(m['z'] - p['z']) / max(p['z'], 1e-12) < TOL_DAMP
                        and mac_c(m['phi'], p['phi']) > TOL_MAC):
                    stable.append({'order': b, 'f': m['f'], 'z': m['z'], 'phi': m['phi']})
                    break
    return stable, per_order


def cluster(stable, tol_hz):
    """Group stable poles into clusters by frequency proximity."""
    if not stable:
        return []
    s = sorted(stable, key=lambda d: d['f'])
    groups, cur = [], [s[0]]
    for p in s[1:]:
        if p['f'] - cur[-1]['f'] <= tol_hz:
            cur.append(p)
        else:
            groups.append(cur)
            cur = [p]
    groups.append(cur)
    return groups


def run_ssi_cov(y, fs, bands, min_stable=4):
    """Assign one mode per band from stable clusters; UNRESOLVED otherwise."""
    stable, per_order = stabilisation(y, fs)
    results = []
    for (lo, hi) in bands:
        inband = [p for p in stable if lo <= p['f'] <= hi]
        if len(inband) < min_stable:
            results.append({'resolved': False,
                            'reason': 'only %d stable poles in band, need %d'
                                      % (len(inband), min_stable),
                            'n_stable': len(inband)})
            continue
        groups = cluster(inband, tol_hz=0.4 * (hi - lo))
        groups = [g for g in groups if len(g) >= min_stable]
        if len(groups) != 1:
            results.append({'resolved': False,
                            'reason': '%d stable clusters in band, expected exactly 1'
                                      % len(groups),
                            'n_clusters': len(groups), 'n_stable': len(inband)})
            continue
        g = groups[0]
        f = float(np.median([p['f'] for p in g]))
        z = float(np.median([p['z'] for p in g]))
        ref = max(g, key=lambda p: p['order'])
        phi = np.real(ref['phi'] * np.exp(-1j * np.angle(ref['phi'][np.argmax(np.abs(ref['phi']))])))
        results.append({'resolved': True, 'f': f, 'z': z, 'n_stable': len(g),
                        'f_scatter': float(np.std([p['f'] for p in g], ddof=1)) if len(g) > 1 else 0.0,
                        'z_scatter': float(np.std([p['z'] for p in g], ddof=1)) if len(g) > 1 else 0.0,
                        'phi': phi})
    return {'bands': results, 'n_stable_total': len(stable),
            'orders': sorted(per_order)}
