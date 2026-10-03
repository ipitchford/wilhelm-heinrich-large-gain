#!/usr/bin/env python3
"""One-line-per-cycle summary of the case_*.json files written by ref_taylor.py."""
import json, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
here = os.path.dirname(os.path.abspath(__file__))
rows = []
for fn in glob.glob(os.path.join(here, 'case_*.json')):
    r = json.load(open(fn))
    rows.append(r)
rows.sort(key=lambda r: (r['mu'], r['rho']))
print(f'{"mu":>6} {"rho":>5} {"T":>10} {"Umin":>9} {"lambda_c=<eta>-mu":>19} {"|dV|rel":>8} {"|dQ|rel":>8} {"<eta0>-mu":>9} {"<eta1>-mu":>9} '
      f'{"minV/r2":>9} {"minQ/s2":>9} {"<eta p>":>10} {"<zeta pd>":>10} {"det":>9} {"idmax":>7} {"ord":>3} {"3.1":>3} {"cvx":>3} {"mono":>4} {"sgn":>3}')
for r in rows:
    lam = r['lambda_c: <eta>-mu']
    dV = abs(r['lambda_c: -<pK/V> (GL)'] - lam) / abs(lam); dQ = abs(r['lambda_c: -<pW/Q> (GL)'] - lam) / abs(lam)
    x = r['maxs']; m = r['mins']
    idmax = max(val for kk, val in x.items() if kk.startswith('id:') or kk.startswith('|'))
    if 'signs' in r:
        sgn = r['signs']['sign pattern'] and r['signs']['V>0'] and r['signs']['Q>0'] and r['signs']['K>0'] and r['signs']['W>0']
        sg = 'ok' if sgn else 'NO'
    else:
        sgn = all(val < 0 for kk, val in x.items() if kk.startswith('sign:')) and m['V/(U^2+p^2)'] > 0 and m['Q/(p^2+pd^2)'] > 0
        sg = 'ok' if sgn else 'fl?'
    p = r['pairs']
    print(f'{r["mu"]:>6g} {r["rho"]:>5g} {r["T"]:>10.4f} {r["Umin"]:>9.2f} {lam:>+19.11e} {dV:>8.0e} {dQ:>8.0e} {r["<eta0>-mu (GL)"]:>+9.0e} {r["<eta1>-mu (GL)"]:>+9.0e} '
          f'{m["V/(U^2+p^2)"]:>9.2e} {m["Q/(p^2+pd^2)"]:>9.2e} {r["<eta p>"]:>+10.3e} {r["<zeta pd>"]:>+10.3e} {r["det"]:>9.2e} {idmax:>7.0e} '
          f'{"ok" if r["orderings_all_positive"] else "NO":>3} {"ok" if p["all_Delta_negative"] else "NO":>3} {"ok" if p["all_convex_negative"] else "NO":>3} {"ok" if p["all_gradient_monotone"] else "NO":>4} {sg:>3}')
print('columns: |dV|rel, |dQ|rel = relative difference of -<pK/V>, -<pW/Q> from <eta>-mu; idmax = largest relative residual of all pointwise identities;')
print('         ord = branch orderings at equal U / equal p (167+167 exact level matches); 3.1 = Delta<0 on 300x299 ordered pairs; cvx = F^(q2-q1)<0 on the same pairs;')
print('         mono = gradient monotonicity on pairs; sgn = V,Q,K,W>0 and (eta-eta0)p<0, (eta-eta1)p<0 at every sample ("fl?" = judged from floats, see VERIFY.md).')
