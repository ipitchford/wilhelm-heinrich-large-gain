#!/usr/bin/env python3
"""Extra referee cases beyond the task grid:
  * rho = 20 at mu = 1, 0.7  -- the two cycles on which the lane first saw Q<0 (NOTES sec. 2 'numerical trap');
  * very small / very large mu;
  * re-run of mu=0.02, rho=12 with exact (arb) sign tests (the first grid run judged signs from underflowing floats)."""
import sys, os, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_taylor as rt

here = os.path.dirname(os.path.abspath(__file__))
allok = True
cases = [('1', '20', 6.0), ('7/10', '20', 6.0), ('1/200', '101/100', 1.0), ('1/200', '3', 3.0), ('50', '2', 2.0), ('1/50', '12', 8.0)]
if len(sys.argv) > 1:
    cases = [c for c in cases if c[0] + ',' + c[1] in sys.argv[1:]]
for mu, rho, guess in cases:
    mu_fr, rho_fr = Fraction(mu), Fraction(rho)
    print(f'=== mu={float(mu_fr):g} rho={float(rho_fr):g} ===', flush=True)
    try:
        r = rt.run(mu_fr, rho_fr, guess=guess, out=os.path.join(here, f'case_mu{float(mu_fr):g}_rho{float(rho_fr):g}.json'))
    except Exception as ex:
        import traceback; traceback.print_exc()
        print('  RESULT: EXCEPTION', repr(ex), flush=True); allok = False; continue
    vs = rt.verdicts(r)
    bad = [nm for nm, ok in vs if not ok]
    allok &= not bad
    m, x = r['mins'], r['maxs']
    idmax = max(val for kk, val in x.items() if kk.startswith('id:') or kk.startswith('|'))
    print(f'  T={r["T"]:.8f} (double {r["T_double"]:.8f}) steps={r["steps"]} Umin={r["Umin"]:.3f} closure={max(r["closure_final"]):.1e}')
    print(f'  lambda_c: <eta>-mu={r["lambda_c: <eta>-mu"]:+.12e}  -<pK/V>={r["lambda_c: -<pK/V> (GL)"]:+.12e}  -<pW/Q>={r["lambda_c: -<pW/Q> (GL)"]:+.12e}  monodromy(double)={r["lam_M_double"]:+.6e} (m_t double {r["m_t_double"]:.2e})')
    print(f'  <eta0>-mu={r["<eta0>-mu (GL)"]:+.1e} <eta1>-mu={r["<eta1>-mu (GL)"]:+.1e}  quad discrepancy GL/series (pK/V, pW/Q, eta0, eta1): '
          + ' '.join(f'{r["quad_discrepancy(GL vs series), rel"][q]:.0e}' for q in ('pK/V', 'pW/Q', 'eta0', 'eta1')))
    print(f'  signs (exact arb tests): {r["signs"]}  direct sign-pattern counts: {r["sign_counts"]}')
    print(f'  log10|.| of extremes (sign, log10): ' + ' | '.join(f'{kk}: ({v[0]:+d}, {v[1]:.2f})' for kk, v in r['log10_abs'].items()))
    print(f'  identity residual max={idmax:.1e}  kappa in [{m["kappa"]:.3e},{x["kappa"]:.3e}] varpi_max={x["varpi"]:.3e}  eta-mu in [{m["eta-mu"]:+.3f},{x["eta-mu"]:+.3f}]'
          f'  eta(crest)-mu={r["eta(crest)-mu"]:+.4f} eta(trough)-mu={r["eta(trough)-mu"]:+.4f}')
    print(f'  <eta p>={r["<eta p>"]:+.6e} <zeta pd>={r["<zeta pd>"]:+.6e} <zeta p>={r["<zeta p>"]:+.6e} <eta pd>={r["<eta pd>"]:+.6e} det={r["det"]:+.4e}')
    print(f'  orderings min: ' + ' | '.join(f'{kk.split(":")[0]}:{kk.split(":")[1].strip()}={v:.2e}' for kk, v in r['orderings(min, must be >0)'].items()))
    print(f'  pairs: {r["pairs"]}')
    print('  CHECKS:', ' '.join(f'{nm}={"ok" if ok else "FAIL"}' for nm, ok in vs))
    print(f'  RESULT: {"ALL OK" if not bad else "FLAGGED: " + ", ".join(bad)}  [{r["seconds"]:.0f}s]', flush=True)
print('EXTRA:', 'ALL OK' if allok else 'SOME FLAGGED')
