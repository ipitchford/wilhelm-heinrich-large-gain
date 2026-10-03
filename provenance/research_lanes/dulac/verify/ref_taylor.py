#!/usr/bin/env python3
"""Referee check of lanes/dulac/NOTES.md (Theorem D, Prop. CD inputs, O1-O2) -- independent of the lane's scripts.

Method: 400-bit Taylor-series integration (order 36, Picard recurrences with python-flint arb_series used as a
multiprecision float: every step is re-centred on exact midpoints) of the state (U,p,y,x) together with the Riccati
pair (eta,zeta).  The cycle is seeded by refcyc.find (own double-precision Newton shooting in log variables) and
polished in high precision.  All derived quantities (V,Q,V3,K,W,...) are computed from their DEFINITIONS in
400-bit arithmetic, so deep troughs (U ~ -8000) cause no cancellation problem.  Not interval-rigorous.

usage: ref_taylor.py MU_FRACTION RHO_FRACTION [U0guess]
"""
import sys, os, math, time, json, bisect
from fractions import Fraction
import numpy as np
from flint import arb, arb_series, arb_poly, ctx, fmpq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import refcyc

PREC, N = 400, 36
ctx.prec = PREC
ctx.cap = N + 1
LOGTOL = -32.0 * math.log(10.0)
NG, PAN = 8, 3
GL = sorted([arb.legendre_p_root(NG, i, weight=True) for i in range(NG)], key=lambda nw: float(nw[0].mid()))   # ascending nodes
INF = float('inf')


def M(v):
    return arb(v.mid())


def lg(v):
    a = abs(v)
    if not (a > 0):
        return -INF
    return float(a.log().mid())


def fl(v):
    """float of an arb (via log, so no overflow/underflow exception)"""
    a = abs(v)
    if not (a > 0):
        return 0.0
    L = float(a.log().mid())
    if L < -700:
        return 0.0
    if L > 700:
        return math.copysign(INF, 1 if v > 0 else -1)
    return float(v.mid())


def sgnlog10(v):
    """(sign, log10|v|) for reporting extreme magnitudes"""
    L = lg(v)
    s = 1 if v > 0 else (-1 if v < 0 else 0)
    return s, (L / math.log(10.0) if L > -INF else -INF)


def taylor(mu, k, st):
    U0, p0, y0, x0, e0, z0 = st
    U = arb_series([U0]); p = arb_series([p0]); y = arb_series([y0])
    x = arb_series([x0]); e = arb_series([e0]); z = arb_series([z0])
    mu1 = mu - 1
    for _ in range(N + 1):
        U, p, y, x, e, z = (p.integral() + U0,
                            (k - mu * p - y).integral() + p0,
                            (x - y).integral() + y0,
                            (x * p).integral() + x0,
                            (z + mu1 * e - e * e).integral() + e0,
                            (x - z - e * z).integral() + z0)
    return U, p, y, x, e, z


FLOORS = [math.log(1e-2), math.log(1e-3), -INF, -INF, -INF, -INF]


def stepsize(S, st, hmax):
    best = INF
    for s, c0, flo in zip(S, st, FLOORS):
        co = s.coeffs()
        sc = max(lg(c0), flo)
        for n in (N - 1, N):
            if n < len(co):
                l = lg(co[n])
                if l > -INF:
                    best = min(best, (LOGTOL + sc - l) / n)
    return min(0.9 * math.exp(best), hmax)


def bisect_root(fun, a, b, iters=170):
    """root of a monotone arb function on [a,b] (arb endpoints); sign change assumed"""
    fa = fun(a)
    sa = fa > 0
    for _ in range(iters):
        m = M((a + b) / 2)
        fm = fun(m)
        if (fm > 0) == sa:
            a = m
        else:
            b = m
    return M((a + b) / 2)


def derived(mu, k, s):
    U, p, y, x, e, z = s
    A = 1 + mu
    pd = k - mu * p - y
    Hd = x - y
    H = y - k
    V = p * p - U * pd
    Q = -pd * H + p * Hd
    V3 = U * Hd - p * H
    K = z * U - e * p - H
    ed = z + (mu - 1) * e - e * e
    zd = x - (1 + e) * z
    W = zd * p - ed * pd
    O = z * p - e * pd - Hd
    kpsi = x * U - x + k
    return dict(U=U, p=p, y=y, x=x, e=e, z=z, pd=pd, Hd=Hd, H=H, V=V, Q=Q, V3=V3, K=K, ed=ed, zd=zd, W=W, O=O, kpsi=kpsi, A=A)


INTEGRANDS = ['one', 'eta', 'pK/V', 'pW/Q', 'eta0', 'eta1', 'eta*p', 'zeta*pd', 'zeta*p', 'eta*pd', 'x', 'p', 'pd^2', 'p^2', 'x*p^2']


def integrands(mu, k, d):
    one = d['V'] * 0 + 1
    return [one, d['e'], d['p'] * d['K'] / d['V'], d['p'] * d['W'] / d['Q'], d['V3'] / d['V'],
            (d['x'] * d['V'] + d['kpsi'] * d['pd']) / d['Q'] - 1, d['e'] * d['p'], d['z'] * d['pd'], d['z'] * d['p'],
            d['e'] * d['pd'], d['x'], d['p'], d['pd'] * d['pd'], d['p'] * d['p'], d['x'] * d['p'] * d['p']]


class Rec:
    def __init__(self, mu, k):
        self.mu, self.k = mu, k
        self.lam = k / mu
        self.steps = []           # (t0, h, polys)
        self.I_gl = [arb(0)] * len(INTEGRANDS)
        self.I_ser = [arb(0)] * len(INTEGRANDS)
        self.mins = {}
        self.maxs = {}
        self.samples = []         # (t, state tuple) at GL nodes
        self.Umin = None; self.Umax = None
        self.sign_counts = dict(resolved=0, unresolved=0, violations=0)

    def mn(self, key, v):
        if key not in self.mins or v < self.mins[key]:
            self.mins[key] = v

    def mx(self, key, v):
        if key not in self.maxs or v > self.maxs[key]:
            self.maxs[key] = v

    def process(self, t0, S, polys, h):
        mu, k = self.mu, self.k
        self.steps.append((t0, h, polys))
        # ---- series-level identities at the step start (derivatives taken from the flow's Taylor series) ----
        d = derived(mu, k, S)
        s0 = [s[0] for s in S]
        d0 = derived(mu, k, s0)
        U, p, y, x, e, z = s0
        V, Q, V3, K, W, pd, kpsi, A = d0['V'], d0['Q'], d0['V3'], d0['K'], d0['W'], d0['pd'], d0['kpsi'], d0['A']
        dV, dQ, dV3, dK, dW = d['V'][1], d['Q'][1], d['V3'][1], d['K'][1], d['W'][1]
        kap_s = d['K'] / d['V']
        kap = K / V; eta0 = V3 / V
        def rel(r, *terms):
            den = sum(abs(t) for t in terms)
            return abs(r) / den if den > 0 else abs(r)
        self.mx('id: V\'=(eta-mu)V+pK', rel(dV - ((e - mu) * V + p * K), dV, (e - mu) * V, p * K))
        self.mx('id: Q\'=(eta-mu)Q+pW', rel(dQ - ((e - mu) * Q + p * W), dQ, (e - mu) * Q, p * W))
        self.mx('id: V\'=-muV+V3', rel(dV - (-mu * V + V3), dV, mu * V, V3))
        self.mx('id: V3\'=-V3+Q+kpsi p', rel(dV3 - (-V3 + Q + kpsi * p), dV3, V3, Q, kpsi * p))
        self.mx('id: Q\'=xV-AQ+kpsi pd', rel(dQ - (x * V - A * Q + kpsi * pd), dQ, x * V, A * Q, kpsi * pd))
        self.mx('id: K\'+(1+eta)K=kpsi', rel(dK + (1 + e) * K - kpsi, dK, (1 + e) * K, kpsi))
        self.mx('id: W\'+(1+eta)W=xp^2', rel(dW + (1 + e) * W - x * p * p, dW, (1 + e) * W, x * p * p))
        self.mx('id: kappa\'', rel(kap_s[1] - (kpsi / V - (1 + e + eta0 - mu) * kap), kap_s[1], kpsi / V, (1 + e + eta0 - mu) * kap))
        # ---- series integrals over [0,h] ----
        for i, g in enumerate(integrands(mu, k, d)):
            G = arb_poly(g.integral().coeffs())
            self.I_ser[i] = M(self.I_ser[i] + G(h))
        # ---- Gauss-Legendre nodes: pointwise checks and quadrature ----
        for j in range(PAN):
            for xi, w in GL:
                for sg in (1, -1) if False else (1,):
                    tau = h * (j + (1 + xi) / 2) / PAN
                    wt = w * h / (2 * PAN)
                    s = [P(tau) for P in polys]
                    dd = derived(mu, k, s)
                    self.point(t0 + tau, s, dd)
                    for i, g in enumerate(integrands(mu, k, dd)):
                        self.I_gl[i] = self.I_gl[i] + wt * g
        self.I_gl = [M(v) for v in self.I_gl]

    def point(self, t, s, d):
        mu, k = self.mu, self.k
        U, p, y, x, e, z = s
        V, Q, V3, K, W, pd, O, kpsi = d['V'], d['Q'], d['V3'], d['K'], d['W'], d['pd'], d['O'], d['kpsi']
        self.mn('V/(U^2+p^2)', V / (U * U + p * p))
        self.mn('Q/(p^2+pd^2)', Q / (p * p + pd * pd))
        self.mn('K/(|zU|+|ep|+|H|)', K / (abs(z * U) + abs(e * p) + abs(d['H'])))
        self.mn('W/(|zd p|+|ed pd|)', W / (abs(d['zd'] * p) + abs(d['ed'] * pd)))
        self.mn('eta-max(0,mu-1)', e - (mu - 1 if mu > 1 else 0))
        self.mn('zeta', z)
        self.mn('(y-zeta)/y', (y - z) / y)
        self.mn('(lam-p-eta)/lam', (self.lam - p - e) / self.lam)
        self.mn('eta0', V3 / V)
        self.mx('|O| rel', abs(O) / (abs(z * p) + abs(e * pd) + abs(d['Hd'])))
        eta0 = V3 / V
        eta1 = (x * V + kpsi * pd) / Q - 1
        zeta0 = Q / V
        kap = K / V
        self.mx('|eta-eta0+pK/V| rel', abs(e - eta0 + p * kap) / (abs(e) + abs(eta0) + abs(p * kap)))
        self.mx('|eta-eta1+pW/Q| rel', abs(e - eta1 + p * W / Q) / (abs(e) + abs(eta1) + abs(p * W / Q)))
        self.mx('|zeta-zeta0+pd K/V| rel', abs(z - zeta0 + pd * kap) / (abs(z) + abs(zeta0) + abs(pd * kap)))
        self.mx('|V3-pK-etaV| rel', abs(V3 - p * K - e * V) / (abs(V3) + abs(p * K) + abs(e * V)))
        # W = W0(z) + kappa*Q  with W0 := x p - (1+eta0-mu)(x-y) + zeta0 (y-k)     (NOTES O1)
        W0 = x * p - (1 + eta0 - mu) * d['Hd'] + zeta0 * d['H']
        self.mx('|W-W0-kappa Q| rel', abs(W - W0 - kap * Q) / (abs(W) + abs(W0) + abs(kap * Q)))
        # sign pattern  (eta-eta0)*sign(p) <= 0 , (eta-eta1)*sign(p) <= 0   (strict off p=0)
        # The direct differences eta-eta0, eta-eta1 are unresolvable at 400 bits where |pK/V| or |pW/Q| is below
        # 1e-120*eta (deep troughs: W ~ e^{-8000}).  So: count certain violations / certain confirmations / unresolved,
        # and separately test the factorised form (K/V>0, W/Q>0 together with the identity residuals above).
        for nm, df in (('eta0', e - eta0), ('eta1', e - eta1)):
            if p > 0:
                good, bad, key, val = df < 0, df > 0, f'sign: (eta-{nm})/(1+eta) on p>0', df / (1 + e)
            elif p < 0:
                good, bad, key, val = df > 0, df < 0, f'sign: ({nm}-eta)/(1+eta) on p<0', -df / (1 + e)
            else:
                continue
            if bad:
                self.sign_counts['violations'] += 1
            elif good:
                self.sign_counts['resolved'] += 1
                self.mx(key, val)
            else:
                self.sign_counts['unresolved'] += 1
        self.mx('kappa', kap); self.mn('kappa', kap)
        self.mx('varpi', W / Q)
        self.mx('eta-mu', e - mu); self.mn('eta-mu', e - mu)
        if self.Umin is None or U < self.Umin:
            self.Umin = U
        if self.Umax is None or U > self.Umax:
            self.Umax = U
        self.samples.append((t, [M(v) for v in s]))


def integrate_period(mu, k, st, hmax, rec=None):
    """from a crest (p=0 exactly) to the next crest. returns (state at next crest, T, nsteps)"""
    t = arb(0)
    n = 0
    while True:
        S = taylor(mu, k, st)
        h = arb(stepsize(S, st, hmax))
        polys = [arb_poly(s.coeffs()) for s in S]
        event = False
        if st[1] > 0:
            pe = polys[1](h)
            if pe <= 0 or not (pe > 0):
                h = bisect_root(polys[1], arb(0), h, iters=PREC + 10)
                event = True
        if rec is not None:
            rec.process(t, S, polys, h)
        st = [M(P(h)) for P in polys]
        t = M(t + h)
        n += 1
        if event:
            st[1] = arb(0)
            return st, t, n
        if n > 400000:
            raise RuntimeError('too many steps')


def eq_data(mu, k):
    muf, kf = float(mu.mid()), float(k.mid())
    r = np.roots([1.0, 1.0 + muf, muf, kf])
    g = -min(v.real for v in r if abs(v.imag) < 1e-12)
    g = arb(g)
    for _ in range(12):                       # Newton in high precision for gamma(gamma-1)(gamma-mu) = k
        f = g * (g - 1) * (g - mu) - k
        fp = (g - 1) * (g - mu) + g * (g - mu) + g * (g - 1)
        g = M(g - f / fp)
    etaE = g - 1
    zetaE = etaE * (etaE + 1 - mu)
    # complex pair: lambda^2 + (A - gamma) lambda + k/gamma = 0  ->  sigma = (gamma - A)/2
    sigma = (g - (1 + mu)) / 2
    lam_u = [v for v in r if v.imag > 1e-12][0]
    return dict(gammaE=g, etaE=M(etaE), zetaE=M(zetaE), sigma=sigma, sigma_np=lam_u.real, omega_np=lam_u.imag)


def run(mu_fr, rho_fr, guess=1.0, out=None, verbose=True):
    t_start = time.time()
    mu = M(arb(fmpq(mu_fr.numerator, mu_fr.denominator)))
    k_fr = rho_fr * mu_fr * (1 + mu_fr)
    k = M(arb(fmpq(k_fr.numerator, k_fr.denominator)))
    muf, kf = float(mu_fr), float(k_fr)
    c = refcyc.find(muf, kf, U0=guess)
    ed = eq_data(mu, k)
    res = dict(mu=muf, rho=float(rho_fr), k=kf, T_double=c['T'], Umin_double=c['Umin'], m_t_double=c['m_t'], lam_M_double=c['lam_M'],
               U0_double=c['U0'], sigma=float(ed['sigma'].mid()), sigma_np=ed['sigma_np'],
               gammaE_minus_A_minus_2sigma_np=float((ed['gammaE'] - (1 + mu)).mid()) - 2 * ed['sigma_np'])
    hmax = c['T'] / 80.0
    U0 = arb(c['U0']); a0 = arb(c['a0'])
    e0, z0 = ed['etaE'], ed['zetaE']
    Jm = c['J'] - np.eye(2)
    closure = []
    # --- polish: quasi-Newton with the double-precision Jacobian; the Riccati pair relaxes meanwhile ---
    for it in range(3):
        y0 = M(k * a0.exp())
        st = [U0, arb(0), y0, M(k * U0.exp()), e0, z0]
        st1, T, n = integrate_period(mu, k, st, hmax)
        rU = st1[0] - U0; ra = (st1[2] / y0).log()
        closure.append(max(abs(float(rU.mid())), abs(float(ra.mid()))))
        dlt = np.linalg.solve(Jm, -np.array([float(rU.mid()), float(ra.mid())]))
        U0 = M(U0 + arb(float(dlt[0]))); a0 = M(a0 + arb(float(dlt[1])))
        e0, z0 = st1[4], st1[5]
        if verbose:
            print(f'   polish {it}: residual {closure[-1]:.2e}  T={float(T.mid()):.10f} steps={n}  [{time.time()-t_start:.0f}s]', flush=True)
        if closure[-1] < 1e-24:
            break
    # --- let the Riccati pair relax to its periodic solution on the polished cycle ---
    y0 = M(k * a0.exp())
    ric = []
    for it in range(14):
        st = [U0, arb(0), y0, M(k * U0.exp()), e0, z0]
        st1, T, n = integrate_period(mu, k, st, hmax)
        dr = max(abs(float((st1[4] / e0).log().mid())), abs(float((st1[5] / z0).log().mid())))
        ric.append(dr)
        e0, z0 = st1[4], st1[5]
        if dr < 1e-30:
            break
    res['riccati_relaxation'] = ric
    if verbose:
        print(f'   riccati relaxation defects: {["%.1e" % v for v in ric]}  [{time.time()-t_start:.0f}s]', flush=True)
    # --- recorded period ---
    st = [U0, arb(0), y0, M(k * U0.exp()), e0, z0]
    rec = Rec(mu, k)
    st1, T, n = integrate_period(mu, k, st, hmax, rec)
    res['closure_polish'] = closure
    res['closure_final'] = [abs(float((st1[0] - st[0]).mid())), abs(float((st1[2] / st[2]).log().mid())),
                            abs(float((st1[4] / st[4]).log().mid())), abs(float((st1[5] / st[5]).log().mid()))]
    res['x_consistency'] = float(abs((st1[3] / (k * st1[0].exp())).log()).mid())
    res['T'] = float(T.mid()); res['steps'] = n
    res['Umin'] = float(rec.Umin.mid()); res['Umax'] = float(rec.Umax.mid())
    Tm = T
    means_gl = {nm: M(v / Tm) for nm, v in zip(INTEGRANDS, rec.I_gl)}
    means_se = {nm: M(v / Tm) for nm, v in zip(INTEGRANDS, rec.I_ser)}
    res['quad_discrepancy(GL vs series), rel'] = {nm: fl(abs(means_gl[nm] - means_se[nm]) / (abs(means_se[nm]) + abs(means_gl[nm]) + arb(10) ** (-300)))
                                                  for nm in INTEGRANDS}
    mg = means_gl
    lam_eta = means_se['eta'] - mu                   # polynomial integrand: series integral is exact to integrator tolerance
    res['lambda_c: <eta>-mu'] = fl(lam_eta)
    res['lambda_c: -<pK/V> (GL)'] = fl(-mg['pK/V']); res['lambda_c: -<pK/V> (series)'] = fl(-means_se['pK/V'])
    res['lambda_c: -<pW/Q> (GL)'] = fl(-mg['pW/Q']); res['lambda_c: -<pW/Q> (series)'] = fl(-means_se['pW/Q'])
    res['<eta0>-mu (GL)'] = fl(mg['eta0'] - mu); res['<eta1>-mu (GL)'] = fl(mg['eta1'] - mu)
    res['<eta p>'] = fl(means_se['eta*p']); res['<zeta pd>'] = fl(means_se['zeta*pd']); res['<zeta p>'] = fl(means_se['zeta*p'])
    res['<eta pd>'] = fl(means_se['eta*pd'])
    det = means_se['eta*p'] * means_se['zeta*pd'] - means_se['zeta*p'] ** 2
    res['det'] = fl(det)
    res['<x>/k-1'] = fl(means_se['x'] / k - 1); res['<p>'] = fl(means_se['p'])
    res['<pd^2>/(mu<p^2>)-1'] = fl(means_se['pd^2'] / (mu * means_se['p^2']) - 1)
    res['<x p^2>/(A<pd^2>)-1'] = fl(means_se['x*p^2'] / ((1 + mu) * means_se['pd^2']) - 1)
    res['mins'] = {kk: fl(v) for kk, v in rec.mins.items()}
    res['mins_log10'] = {kk: sgnlog10(v) for kk, v in rec.mins.items() if kk in ('V/(U^2+p^2)', 'Q/(p^2+pd^2)')}
    res['maxs'] = {kk: fl(v) for kk, v in rec.maxs.items()}
    # exact (arb) sign tests -- floats underflow in deep troughs (quantities ~ e^{-8000}), so never judge signs from fl()
    res['signs'] = {'V>0': bool(rec.mins['V/(U^2+p^2)'] > 0), 'Q>0': bool(rec.mins['Q/(p^2+pd^2)'] > 0),
                    'K>0': bool(rec.mins['K/(|zU|+|ep|+|H|)'] > 0), 'W>0': bool(rec.mins['W/(|zd p|+|ed pd|)'] > 0),
                    'eta>max(0,mu-1)': bool(rec.mins['eta-max(0,mu-1)'] > 0), 'zeta>0': bool(rec.mins['zeta'] > 0),
                    'y>zeta': bool(rec.mins['(y-zeta)/y'] > 0), 'eta<lam-p': bool(rec.mins['(lam-p-eta)/lam'] > 0),
                    'kappa>0': bool(rec.mins['kappa'] > 0),
                    'sign pattern': bool(rec.sign_counts['violations'] == 0 and all(v < 0 for kk, v in rec.maxs.items() if kk.startswith('sign:')))}
    res['sign_counts'] = dict(rec.sign_counts)
    res['log10_abs'] = {kk: sgnlog10(v) for kk, v in list(rec.mins.items()) + [(q, w) for q, w in rec.maxs.items() if q.startswith('sign:')]}
    res['nsamples'] = len(rec.samples)
    if verbose:
        print(f'   recorded: T={res["T"]:.10f} steps={n} samples={len(rec.samples)} [{time.time()-t_start:.0f}s]', flush=True)

    # ------------------------------------------------------------------ global-time evaluation
    t0s = [float(s[0].mid()) for s in rec.steps]

    def ev(t):
        i = bisect.bisect_right(t0s, float(t.mid())) - 1
        i = max(0, min(i, len(rec.steps) - 1))
        t0, h, polys = rec.steps[i]
        tau = t - t0
        return [P(tau) for P in polys]

    def crossing_time(idx_fun, ta, tb):
        """zero of a function with a sign change between ta and tb"""
        return bisect_root(idx_fun, ta, tb)

    # locate special times by scanning the samples
    ts = [s[0] for s in rec.samples]
    def scan(fun, sign_from):
        """first sample interval where fun changes from sign_from to -sign_from"""
        prev = None
        for t, s in rec.samples:
            v = fun(s)
            sg = 1 if v > 0 else -1
            if prev is not None and prev[1] == sign_from and sg == -sign_from:
                return prev[0], t
            prev = (t, sg)
        return None
    pfun = lambda s: s[1]
    pdfun = lambda s: k - mu * s[1] - s[2]
    iv = scan(pfun, -1); t_tr = crossing_time(lambda t: ev(t)[1], iv[0], iv[1])
    iv = scan(pdfun, -1); t1 = crossing_time(lambda t: pdfun(ev(t)), iv[0], iv[1])       # p minimum
    iv = scan(pdfun, +1)
    # the first + -> - change of pd after t1 is the p maximum
    prev = None; iv2 = None
    for t, s in rec.samples:
        if not (t > t1):
            continue
        sg = 1 if pdfun(s) > 0 else -1
        if prev is not None and prev[1] == 1 and sg == -1:
            iv2 = (prev[0], t); break
        prev = (t, sg)
    t2 = crossing_time(lambda t: pdfun(ev(t)), iv2[0], iv2[1])                               # p maximum
    res['times'] = dict(t_pmin=float(t1.mid()), t_trough=float(t_tr.mid()), t_pmax=float(t2.mid()), T=res['T'])
    order_ok = bool(arb(0) < t1 and t1 < t_tr and t_tr < t2 and t2 < T)
    res['phase order crest<pmin<trough<pmax<crest'] = order_ok
    # number of sign changes of p and pd over the samples (unimodality inputs from v2.0 Thm 4.1)
    def nchanges(fun):
        c_ = 0; prev = None
        for t, s in rec.samples:
            sg = 1 if fun(s) > 0 else -1
            if prev is not None and sg != prev:
                c_ += 1
            prev = sg
        return c_
    res['sign changes of p (interior)'] = nchanges(pfun); res['sign changes of pd'] = nchanges(pdfun)
    Utr = ev(t_tr)[0]; Ucr = st[0]
    pmin = ev(t1)[1]; pmax = ev(t2)[1]
    res['p_min'] = fl(pmin); res['p_max'] = fl(pmax); res['lam'] = fl(k / mu)
    res['eta(crest)-mu'] = fl(st[4] - mu); res['eta(trough)-mu'] = fl(ev(t_tr)[4] - mu)

    # ------------------------------------------------------------------ equal-U branch orderings
    nl = 160
    us = [0.5 * (1 - math.cos(math.pi * j / (nl + 1))) for j in range(1, nl + 1)] + [1e-3, 1e-6, 1e-9, 1e-12, 1 - 1e-3, 1 - 1e-6, 1 - 1e-9]
    worst = {}
    def upd(key, v):
        if key not in worst or v < worst[key]:
            worst[key] = v
    for u in us:
        cU = M(Ucr - (Ucr - Utr) * arb(u))
        tl = bisect_root(lambda t: ev(t)[0] - cU, arb(0), t_tr)      # lower branch (p<0), U decreasing
        tu = bisect_root(lambda t: ev(t)[0] - cU, t_tr, T)           # upper branch (p>0), U increasing
        sl, su = ev(tl), ev(tu)
        Ul, pl, yl, xl, el, zl = sl; Uu, pu, yu, xu, eu, zu = su
        upd('equalU: (eta_low-eta_up)/(sum)', (el - eu) / (el + eu))
        upd('equalU: (y_low-y_up)/(sum)', (yl - yu) / (yl + yu))
        dp = pu - pl
        upd('equalU: p_up-p_low', dp)
        upd('equalU: [(y_low-y_up) - eta_up dp]/(y_low+y_up)', ((yl - yu) - eu * dp) / (yl + yu))
        upd('equalU: [eta_low dp - (y_low-y_up)]/(y_low+y_up)', (el * dp - (yl - yu)) / (yl + yu))
    # ------------------------------------------------------------------ equal-p branch orderings
    for u in us:
        cp = M(pmax - (pmax - pmin) * arb(u))
        tL = bisect_root(lambda t: ev(t)[1] - cp, t1, t2)            # trough side, p increasing
        if cp >= 0:
            tR = bisect_root(lambda t: ev(t)[1] - cp, t2, T)
        else:
            tR = bisect_root(lambda t: ev(t)[1] - cp, arb(0), t1)
        sL, sR = ev(tL), ev(tR)
        UL, pL, yL, xL, eL, zL = sL; UR, pR, yR, xR, eR, zR = sR
        dU = UR - UL
        upd('equalp: U_right-U_left', dU)
        upd('equalp: (y_r-y_l)/(sum)', (yR - yL) / (yR + yL))
        upd('equalp: (zeta_r-zeta_l)/(sum)', (zR - zL) / (zR + zL))
        upd('equalp: [(y_r-y_l) - zeta_l dU]/(y_r+y_l)', ((yR - yL) - zL * dU) / (yR + yL))
        upd('equalp: [zeta_r dU - (y_r-y_l)]/(y_r+y_l)', (zR * dU - (yR - yL)) / (yR + yL))
    res['orderings(min, must be >0)'] = {kk: fl(v) for kk, v in worst.items()}
    res['orderings_all_positive'] = bool(all(v > 0 for v in worst.values()))
    if verbose:
        print(f'   orderings done [{time.time()-t_start:.0f}s]', flush=True)

    # ------------------------------------------------------------------ pair tests (Thm 3.1 on one cycle, convexity, O2)
    ns = len(rec.samples)
    npair = 300
    idx = sorted(set(int(round(i * (ns - 1) / (npair - 1))) for i in range(npair)))
    pts = []
    for i in idx:
        t, s = rec.samples[i]
        d = derived(mu, k, s)
        pts.append((s[0], s[1], s[2], s[4], s[5], d['pd'], d['V3'] / d['V'], d['Q'] / d['V'], d['K'] / d['V']))
    mxD = None; mxC = None; mnMono = None; mxId = None; allD = True; allC = True; allM = True
    for ia, (UA, pA, yA, eA, zA, pdA, e0A, z0A, kA) in enumerate(pts):
        for ib, (UB, pB, yB, eB, zB, pdB, e0B, z0B, kB) in enumerate(pts):
            if ia == ib:
                continue
            dU = UB - UA; dp = pB - pA; dy = yB - yA
            Dl = zA * dU - eA * dp - dy
            scD = abs(zA * dU) + abs(eA * dp) + abs(dy)
            cv = pA * dp - pdA * dU
            scC = abs(pA * dp) + abs(pdA * dU)
            if not (Dl < 0):
                allD = False
            if not (cv < 0):
                allC = False
            rD = Dl / scD; rC = cv / scC
            if mxD is None or rD > mxD:
                mxD = rD
            if mxC is None or rC > mxC:
                mxC = rC
            # O2 identity: Delta = [zeta0 dU - eta0 dp - dy] + kappa * (F ^ dq)
            idr = abs(Dl - ((z0A * dU - e0A * dp - dy) + kA * cv)) / (scD + abs(kA * cv))
            if mxId is None or idr > mxId:
                mxId = idr
            if ib > ia:
                mono = (zA - zB) * (UA - UB) - (eA - eB) * (pA - pB)
                scM = abs((zA - zB) * (UA - UB)) + abs((eA - eB) * (pA - pB))
                if not (mono > 0):
                    allM = False
                rM = mono / scM
                if mnMono is None or rM < mnMono:
                    mnMono = rM
    res['pairs'] = dict(n=len(pts), all_Delta_negative=allD, max_Delta_rel=fl(mxD), all_convex_negative=allC, max_FwedgeDq_rel=fl(mxC),
                        all_gradient_monotone=allM, min_monotone_rel=fl(mnMono), max_O2_identity_residual=fl(mxId))
    res['seconds'] = time.time() - t_start
    if out:
        with open(out, 'w') as fh:
            json.dump(res, fh, indent=1, default=str)
    return res


def verdicts(r):
    """list of (name, bool) referee checks for one cycle"""
    m, x = r['mins'], r['maxs']
    lam = r['lambda_c: <eta>-mu']
    tol = 1e-9 * max(1.0, abs(lam)) + 2e-7 * abs(lam)
    v = []
    v.append(('closure<1e-15', max(r['closure_final']) < 1e-15))
    v.append(('O=0 (<1e-20)', x['|O| rel'] < 1e-20))
    sg = r['signs']
    v.append(('V>0', sg['V>0'])); v.append(('Q>0', sg['Q>0']))
    v.append(('K>0', sg['K>0'])); v.append(('W>0', sg['W>0']))
    v.append(('identities<1e-20', max(val for kk, val in x.items() if kk.startswith('id:') or kk.startswith('|')) < 1e-20))
    v.append(('lam_c(V) agrees', abs(r['lambda_c: -<pK/V> (GL)'] - lam) < tol))
    v.append(('lam_c(Q) agrees', abs(r['lambda_c: -<pW/Q> (GL)'] - lam) < tol))
    v.append(('<eta0>=mu', abs(r['<eta0>-mu (GL)']) < tol)); v.append(('<eta1>=mu', abs(r['<eta1>-mu (GL)']) < tol))
    v.append(('sign pattern', sg['sign pattern']))
    v.append(('<eta p><0', r['<eta p>'] < 0)); v.append(('<zeta pd><0', r['<zeta pd>'] < 0)); v.append(('det>0', r['det'] > 0))
    v.append(('<eta pd>=<zeta p>', abs(r['<eta pd>'] - r['<zeta p>']) < 1e-12 * (abs(r['<zeta p>']) + 1e-30) + 1e-25))
    v.append(('orderings', r['orderings_all_positive']))
    v.append(('pairs: Thm3.1', r['pairs']['all_Delta_negative'])); v.append(('pairs: convex', r['pairs']['all_convex_negative']))
    v.append(('pairs: monotone grad', r['pairs']['all_gradient_monotone']))
    v.append(('y>zeta', sg['y>zeta']))
    if r['mu'] <= 0.5:
        v.append(('eta<lam-p (mu<=1/2)', sg['eta<lam-p']))
    v.append(('lam_c<0', lam < 0))
    if r['m_t_double'] > 1e-5:
        # the double-precision section derivative is a central difference with absolute noise ~1e-7 in m_t
        v.append(('lam_c vs monodromy (double)', abs(r['lam_M_double'] - lam) * r['T'] < 5e-7 / r['m_t_double'] + 1e-7))
    return v


if __name__ == '__main__':
    mu_fr = Fraction(sys.argv[1]); rho_fr = Fraction(sys.argv[2])
    guess = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'case_mu{float(mu_fr):g}_rho{float(rho_fr):g}.json')
    print(f'=== mu={float(mu_fr):g} rho={float(rho_fr):g} ===', flush=True)
    r = run(mu_fr, rho_fr, guess=guess, out=out)
    print('  max quad discrepancy (GL vs series):', max(r['quad_discrepancy(GL vs series), rel'].values()), r['quad_discrepancy(GL vs series), rel'])
    for kk in ('T', 'T_double', 'steps', 'Umin', 'Umax', 'closure_polish', 'closure_final', 'x_consistency', 'lambda_c: <eta>-mu', 'lambda_c: -<pK/V> (GL)',
               'lambda_c: -<pK/V> (series)', 'lambda_c: -<pW/Q> (GL)', 'lambda_c: -<pW/Q> (series)', 'lam_M_double', '<eta0>-mu (GL)', '<eta1>-mu (GL)',
               '<eta p>', '<zeta pd>', '<zeta p>', '<eta pd>', 'det', '<x>/k-1', '<p>', '<pd^2>/(mu<p^2>)-1', '<x p^2>/(A<pd^2>)-1',
               'mins', 'mins_log10', 'maxs', 'times', 'phase order crest<pmin<trough<pmax<crest', 'sign changes of p (interior)', 'sign changes of pd',
               'orderings(min, must be >0)', 'pairs', 'signs', 'sign_counts', 'log10_abs', 'eta(crest)-mu', 'eta(trough)-mu', 'riccati_relaxation', 'sigma', 'gammaE_minus_A_minus_2sigma_np', 'seconds'):
        print(f'  {kk}: {r[kk]}')
    vs = verdicts(r)
    bad = [nm for nm, ok in vs if not ok]
    print('  CHECKS:', ' '.join(f'{nm}={"ok" if ok else "FAIL"}' for nm, ok in vs))
    print('  RESULT:', 'ALL OK' if not bad else 'FLAGGED: ' + ', '.join(bad), flush=True)
