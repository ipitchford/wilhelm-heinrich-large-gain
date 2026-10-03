#!/usr/bin/env python3
"""Job 2(a): Zhang-type test.  On the crest half-axis {v=0,u>0} of the reduced planar flow let R1 be the return map
and iota(u) = int over the turn of (eta-mu) dt.  In the energy coordinate e = G~(u) = int_0^u (h(s,0)-k) ds the
return map has derivative exactly exp(iota).  Uniqueness follows if iota is strictly decreasing (concave energy map).
We sample iota on many W^u(E) spirals (event-refined) and test monotonicity; same on the trough half-axis and for
the two half-turn maps (crest->trough, trough->crest)."""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, equilibrium_data, kH

def ext(mu, k):
    def f(t, s):
        U, p, y, eta, zeta, I = s
        x = k*np.exp(U)
        return [p, k-mu*p-y, x-y, zeta+(mu-1)*eta-eta*eta, x-(1+eta)*zeta, eta-mu]
    return f
def ev_c(t, s): return s[1]
ev_c.direction = -1
def ev_t(t, s): return s[1]
ev_t.direction = 1

def turns(mu, k, nth=12, r0=1e-6):
    ed = equilibrium_data(mu, k)
    tend = (np.log(1/r0)+16)/ed['lam_u'].real
    C = []; Tr = []; HC = []; HT = []
    for th in np.linspace(0, 2*np.pi, nth, endpoint=False):
        v = np.real(np.exp(1j*th)*ed['vec_u']); z0 = np.array([0, 0, k])+r0*v/np.linalg.norm(v)
        sol = solve_ivp(ext(mu, k), [0, tend], list(z0)+[ed['etaE'], ed['zetaE'], 0.0], method='DOP853', rtol=1e-12, atol=1e-15, events=[ev_c, ev_t])
        tc, yc = sol.t_events[0], sol.y_events[0]; tt, yt = sol.t_events[1], sol.y_events[1]
        for i in range(len(tc)-1):
            C.append((yc[i][0], yc[i+1][0], yc[i+1][5]-yc[i][5], tc[i+1]-tc[i], yc[i][2]))
        for i in range(len(tt)-1):
            Tr.append((yt[i][0], yt[i+1][0], yt[i+1][5]-yt[i][5], tt[i+1]-tt[i], yt[i][2]))
        # half turns: crest -> next trough, trough -> next crest
        ev = sorted([(t, 0, y) for t, y in zip(tc, yc)]+[(t, 1, y) for t, y in zip(tt, yt)], key=lambda a: a[0])
        for a, b in zip(ev[:-1], ev[1:]):
            if a[1] == 0 and b[1] == 1: HC.append((a[2][0], b[2][0], b[2][5]-a[2][5]))
            if a[1] == 1 and b[1] == 0: HT.append((a[2][0], b[2][0], b[2][5]-a[2][5]))
    return np.array(C), np.array(Tr), np.array(HC), np.array(HT), ed

def mono_report(u, val, umin, tol_sep=1e-5):
    """val as function of |u| sorted increasing; report whether val is decreasing, ignoring pairs closer than tol_sep"""
    o = np.argsort(np.abs(u)); a = np.abs(u)[o]; v = val[o]
    m = a > umin; a = a[m]; v = v[m]
    keep = [0]
    for i in range(1, len(a)):
        if a[i]-a[keep[-1]] > tol_sep: keep.append(i)
    a = a[keep]; v = v[keep]
    d = np.diff(v)
    return dict(n=int(len(a)), decreasing=bool(np.all(d < 0)), max_increment=float(d.max()), first=float(v[0]), last=float(v[-1]),
                n_increasing_pairs=int((d >= 0).sum()))

def run(mu, rho):
    k = rho*kH(mu)
    C, Tr, HC, HT, ed = turns(mu, k)
    cyc = find_cycle(mu, k)
    res = dict(mu=mu, rho=rho, iota_at_E=float(2*np.pi*(ed['etaE']-mu)/ed['lam_u'].imag), log_m_t=float(np.log(cyc['m_t'])))
    res['crest'] = mono_report(C[:, 0], C[:, 2], 1e-3)
    res['trough'] = mono_report(Tr[:, 0], Tr[:, 2], 1e-3)
    res['half_crest_to_trough'] = mono_report(HC[:, 0], HC[:, 2], 1e-3)
    res['half_trough_to_crest'] = mono_report(HT[:, 0], HT[:, 2], 1e-3)
    # R1(u)/u decreasing?  (sublinearity in the raw coordinate u) and in coordinate y-k on the crest
    res['R1_over_u_decreasing'] = mono_report(C[:, 0], C[:, 1]/C[:, 0], 1e-3)
    return res, C

if __name__ == '__main__':
    out = []
    grid = [(1.0, 1.02), (1.0, 1.1), (1.0, 1.5), (1.0, 2.0), (1.0, 3.0), (0.3, 1.1), (0.3, 1.5), (0.3, 2.5), (0.1, 1.2), (0.1, 2.0),
            (3.0, 1.3), (3.0, 2.0), (10.0, 1.5), (0.03, 1.5)]
    for mu, rho in grid:
        try:
            r, C = run(mu, rho)
        except Exception as e:
            print(mu, rho, 'FAILED', repr(e)); continue
        out.append(r)
        print(f"mu={mu} rho={rho}: iota(E)={r['iota_at_E']:+.4f} log m_t={r['log_m_t']:+.4f} | crest: n={r['crest']['n']} decr={r['crest']['decreasing']} maxinc={r['crest']['max_increment']:.1e} "
              f"| trough: decr={r['trough']['decreasing']} maxinc={r['trough']['max_increment']:.1e} | half c->t decr={r['half_crest_to_trough']['decreasing']} ({r['half_crest_to_trough']['first']:+.3f}->{r['half_crest_to_trough']['last']:+.3f}) "
              f"| half t->c decr={r['half_trough_to_crest']['decreasing']} ({r['half_trough_to_crest']['first']:+.3f}->{r['half_trough_to_crest']['last']:+.3f}) | R1/u decr={r['R1_over_u_decreasing']['decreasing']}", flush=True)
    json.dump(out, open('job2_a.json', 'w'), indent=1)
