#!/usr/bin/env python3
"""Job 2(c)/(Dulac): pointwise sign tests of  DLambda + eta - mu  on the attractor disc for explicit candidates,
and an LP search over a 14-function basis for a Lambda with  DLambda + eta - mu <= t < 0  on the sampled disc."""
import sys, json
import numpy as np
from scipy.optimize import linprog
from job2_lib import sample
from whlib import kH, find_cycle

def features(S):
    mu, k, ed = S['mu'], S['k'], S['ed']
    U, p, y, x, pp, eta, zeta, K, W = (S[n] for n in ('U', 'p', 'y', 'x', 'pp', 'eta', 'zeta', 'K', 'W'))
    etap, zetap = S['etap'], S['zetap']
    lam = k/mu; A = 1+mu; gE = ed['gammaE']; etaE = ed['etaE']; zetaE = ed['zetaE']; sig = 0.5*(etaE-mu)
    G = U*np.exp(U)-np.exp(U)+1
    LamE = (y-k)-zetaE*U+etaE*p
    Kp = zetap*U-etap*p
    Q = W-2*sig*Kp+2*zetaE*K
    SQ = x*p*p-2*sig*x*U*p+2*zetaE*k*G+2*sig*etap*K
    D = {
        'logW': x*p*p/W-(1+eta), 'logK': k*G/K-(1+eta), 'logzeta': x/zeta-(1+eta), 'logeta': zeta/eta+mu-1-eta,
        'logy': x/y-1, 'log(lam-p)': y/(lam-p)-mu, 'U': p, 'p': pp, 'y/k': (x-y)/k,
        'logLamE': k*(np.exp(U)-1-U)/LamE-gE, 'eta': etap, 'zeta': zetap, 'log(1+eta)': etap/(1+eta),
    }
    if Q.min() > 0:
        D['logQ'] = SQ/Q-(1+eta)
    return D, dict(minQ=float(Q.min()))

def run(mu, rho, rfrac):
    k = rho*kH(mu)
    S = sample(mu, k, np.linspace(0, 2*np.pi, 8, endpoint=False), n=2500)
    cyc = find_cycle(mu, k)
    D, info = features(S)
    eta = S['eta']; base = eta-mu
    m = S['r'] > rfrac*cyc['U0']
    res = dict(mu=mu, rho=rho, rfrac=rfrac, n=int(m.sum()), **info)
    # explicit candidates:  coefficient dict -> expression
    cands = {
        'L0: 0': {},
        'L1: -logW+2logzeta': {'logW': -1, 'logzeta': 2},
        'L2: -logK+2logzeta': {'logK': -1, 'logzeta': 2},
        'L3: -logLamE+logzeta': {'logLamE': -1, 'logzeta': 1},
        'L4: -logW+logzeta+logeta': {'logW': -1, 'logzeta': 1, 'logeta': 1},
        'L5: -logK+logzeta+logeta': {'logK': -1, 'logzeta': 1, 'logeta': 1},
        'L6: logzeta-logy-log(lam-p)': {'logzeta': 1, 'logy': -1, 'log(lam-p)': -1},
        'L7: logzeta': {'logzeta': 1}, 'L8: logeta': {'logeta': 1},
    }
    if 'logQ' in D:
        cands['L9: -logQ+2logzeta'] = {'logQ': -1, 'logzeta': 2}
        cands['L10: -2logQ+3logzeta'] = {'logQ': -2, 'logzeta': 3}
    res['candidates'] = {}
    for name, co in cands.items():
        e = base.copy()
        for kk, c in co.items(): e = e+c*D[kk]
        e = e[m]
        res['candidates'][name] = dict(min=float(e.min()), max=float(e.max()), frac_pos=float(np.mean(e > 0)))
    # LP
    names = list(D.keys())
    idx = np.nonzero(m)[0]
    rng = np.random.default_rng(0)
    if len(idx) > 6000: idx = rng.choice(idx, 6000, replace=False)
    Amat = np.column_stack([D[nm][idx] for nm in names]+[-np.ones(len(idx))])
    b = -base[idx]
    c = np.zeros(len(names)+1); c[-1] = 1
    bounds = [(-30, 30)]*len(names)+[(-50, 50)]
    lp = linprog(c, A_ub=Amat, b_ub=b, bounds=bounds, method='highs')
    res['LP_status'] = int(lp.status); res['LP_t'] = float(lp.x[-1]) if lp.status == 0 else None
    res['LP_coeffs'] = {nm: float(v) for nm, v in zip(names, lp.x[:-1])} if lp.status == 0 else None
    # validate on all points
    if lp.status == 0:
        e = base.copy()
        for nm, v in zip(names, lp.x[:-1]): e = e+v*D[nm]
        res['LP_validation_max_all_points'] = float(e[m].max())
    return res

if __name__ == '__main__':
    out = []
    for mu, rho in [(1.0, 1.1), (1.0, 2.0), (0.3, 1.5), (3.0, 1.3)]:
        for rfrac in (0.002, 0.3):
            r = run(mu, rho, rfrac); out.append(r)
            print(f"mu={mu} rho={rho} r>{rfrac}*crest n={r['n']} minQ={r['minQ']:.2e}")
            for nm, v in r['candidates'].items():
                print(f"    {nm:32s} min={v['min']:+9.3f} max={v['max']:+9.3f} frac>0={v['frac_pos']:.3f}")
            print(f"    LP: status={r['LP_status']} t*={r['LP_t']} validation max={r.get('LP_validation_max_all_points')}")
            if r['LP_coeffs']:
                print('       coeffs:', {a: round(b, 3) for a, b in r['LP_coeffs'].items() if abs(b) > 1e-6}, flush=True)
    json.dump(out, open('job2_c.json', 'w'), indent=1)
