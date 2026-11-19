"""semliki_dynrating.py -- a dynamic stage-discharge rating: a random-walk stage shift on an installed rating,
estimated from gaugings with a Kalman filter and smoother, plus the Jones loop correction for unsteady flow.

Part of the reproducibility package of Kakooza, C., Nafula, M., Kasibante, G., Waiswa, F. and Zaake, B. (2026),
"A dynamic stage-discharge rating for a migrating tropical river", preprint, https://doi.org/10.5281/zenodo.24367492.

THE MODEL
  Channel change, a time-varying stage shift delta(t) on a base rating Q_base(h):
      Q(h, t) = Q_base(h - delta(t))
      delta(t) = delta(t_prev) + w,   w ~ N(0, q * dt_years)                     random walk
      gauging i:  h_i - h_base(Q_i) = delta(t_i) + e_i,   e_i ~ N(0, (sigma_Q * dh/dlnQ)^2)
  A Kalman filter gives delta and its variance at every gauging; between gaugings the mean is held and the
  variance grows by q per year. q and the gauging error sigma_Q are estimated by maximum likelihood.
  Unsteady flow, the Jones loop correction:
      Q = Q_steady(h) * sqrt(1 + (dh/dt) / (S0 * c)),   c = (dQ/dh) / B        kinematic-wave celerity
Licence: MIT (see LICENSE).
"""
import math
import numpy as np

YEAR = 365.25
Z90 = 1.6449


def kalman(t, d, sd, q, x0, p0):
    """1-D random walk + noise.  t (years), d (observed shifts, m), sd (their error, m), q (m2/yr),
    x0, p0: prior mean and variance.  -> loglik, [(x_pred, P_pred, innovation, S, x_filt, P_filt)]"""
    x, P, tl, ll, rec = x0, p0, t[0], 0.0, []
    for i in range(len(d)):
        Pp = P + q * (t[i] - tl); tl = t[i]
        S = Pp + sd[i] ** 2; v = d[i] - x; K = Pp / S
        ll += -0.5 * (math.log(2 * math.pi * S) + v * v / S)
        xp = x
        x = x + K * v; P = (1 - K) * Pp
        rec.append((xp, Pp, v, S, x, P))
    return ll, rec


def smooth(t, rec, q):
    """Rauch-Tung-Striebel smoother -> smoothed means and variances at the gaugings"""
    n = len(rec)
    xs = [r[4] for r in rec]; Ps = [r[5] for r in rec]
    for i in range(n - 2, -1, -1):
        Pp = rec[i][5] + q * (t[i + 1] - t[i])
        J = rec[i][5] / Pp if Pp > 0 else 0.0
        xs[i] = rec[i][4] + J * (xs[i + 1] - rec[i][4])
        Ps[i] = rec[i][5] + J * J * (Ps[i + 1] - Pp)
    return xs, Ps


def power_rating(a, h0, c):
    """Q = a (h - h0)^c and its inverse h = h0 + (Q/a)^(1/c); dh/dlnQ = (h - h0)/c"""
    q_of_h = lambda h: a * np.clip(np.asarray(h, float) - h0, 0, None) ** c
    h_of_q = lambda Q: h0 + (np.asarray(Q, float) / a) ** (1.0 / c)
    dh_dlnq = lambda h: (np.asarray(h, float) - h0) / c
    return q_of_h, h_of_q, dh_dlnq


def fit_drift(t, h, Q, h_of_q, dh_dlnq, prior_sd=0.1, verbose=True):
    """maximum likelihood of the drift rate q (m2/yr) and the relative gauging error, grid then refinement.
    t in years, h observed stages, Q measured flows; the base rating supplies h_of_q and dh_dlnq."""
    d = [hi - float(h_of_q(qi)) for hi, qi in zip(h, Q)]
    k = [float(dh_dlnq(hi)) for hi in h]

    def ll(lq, ls):
        s = math.exp(ls)
        return kalman(t, d, [s * x for x in k], math.exp(lq), d[0], prior_sd ** 2)[0]
    best = None
    for lq in np.linspace(math.log(1e-5), math.log(0.1), 41):
        for ls in np.linspace(math.log(0.01), math.log(0.2), 31):
            v = ll(lq, ls)
            if best is None or v > best[0]:
                best = (v, lq, ls)
    v, lq, ls = best
    step = 0.25
    for _ in range(40):
        moved = False
        for dq, ds in ((step, 0), (-step, 0), (0, step), (0, -step)):
            w = ll(lq + dq, ls + ds)
            if w > v + 1e-9:
                v, lq, ls, moved = w, lq + dq, ls + ds, True
        if not moved:
            step /= 2
            if step < 1e-4:
                break
    grid = np.linspace(math.log(0.02), math.log(0.3), 200)
    static = max(ll(math.log(1e-12), x) for x in grid)
    s_static = math.exp(grid[int(np.argmax([ll(math.log(1e-12), x) for x in grid]))])
    q, s = math.exp(lq), math.exp(ls)
    _, rec = kalman(t, d, [s * x for x in k], q, d[0], prior_sd ** 2)
    z = [r[2] / math.sqrt(r[3]) for r in rec][5:]
    xs, Ps = smooth(t, rec, q)
    out = dict(q_m2_per_yr=q, drift_sd_m_per_sqrt_yr=math.sqrt(q), sigma_gauging=s, loglik=v, loglik_static=static,
               sigma_static=s_static, aic=-2 * v + 4, aic_static=-2 * static + 2,
               z_sd=float(np.std(z)), cover90=float(np.mean(np.abs(z) < Z90)),
               obs_shift=d, smoothed_shift=xs, smoothed_sd=[math.sqrt(p) for p in Ps], filter=rec)
    if verbose:
        print("drift rate %.3f m per sqrt(year); gauging error %.1f %%; loglik %.1f (fixed rating %.1f); AIC %.0f vs %.0f; "
              "one-step z sd %.2f, %.0f %% inside the 90 %% band"
              % (out["drift_sd_m_per_sqrt_yr"], 100 * s, v, static, out["aic"], out["aic_static"], out["z_sd"], 100 * out["cover90"]))
    return out


class Rating:
    """a base rating table (h ascending, Q) with a stage shift and its standard deviation"""
    def __init__(self, base_h, base_q, shift, shift_sd):
        self.h, self.q, self.shift, self.sd = np.asarray(base_h, float), np.asarray(base_q, float), shift, shift_sd

    def _q(self, h, sh):
        return np.interp(np.asarray(h, float) - sh, self.h, self.q, left=0.0, right=np.nan)

    def flow(self, h, dhdt=None, slope=None, width=None):
        """stage -> (Q, Q5, Q95, Jones factor)"""
        q = self._q(h, self.shift)
        lo = self._q(h, self.shift + Z90 * self.sd)
        hi = self._q(h, self.shift - Z90 * self.sd)
        f = jones(self, h, dhdt, slope, width) if dhdt is not None else 1.0
        return q * f, lo * f, hi * f, f

    def stage(self, Q):
        """flow -> (h, h5, h95)"""
        hb = float(np.interp(Q, self.q, self.h))
        return hb + self.shift, hb + self.shift - Z90 * self.sd, hb + self.shift + Z90 * self.sd

    def dqdh(self, h, e=0.01):
        return (float(self._q(h + e, self.shift)) - float(self._q(h - e, self.shift))) / (2 * e)


def jones(R, h, dhdt_m_per_day, slope, width, cap=1.0):
    """Jones loop factor Q/Q_steady for a rate of stage change (m/day), capped at +-cap m/day"""
    if dhdt_m_per_day is None:
        return 1.0
    r = max(-cap, min(cap, float(dhdt_m_per_day))) / 86400.0
    ck = R.dqdh(h) / width
    if not ck or ck <= 0:
        return 1.0
    return math.sqrt(max(0.05, 1.0 + r / (slope * ck)))


def band_growth(sd_last, q, days_since, months=(0, 3, 6, 12)):
    """90 % half-width of the stage band (m) if no gauging is made for the given extra months"""
    return {m: Z90 * math.sqrt(sd_last ** 2 + q * (days_since / YEAR + m / 12.0)) for m in months}


def band_after_gauging(sd_now, sd_gauging):
    """90 % half-width after one new gauging of error sd_gauging (m of stage)"""
    return Z90 * math.sqrt(1.0 / (1.0 / sd_now ** 2 + 1.0 / sd_gauging ** 2))
