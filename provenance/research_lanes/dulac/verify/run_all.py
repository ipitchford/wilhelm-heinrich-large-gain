#!/usr/bin/env python3
"""Referee grid: mu in {0.02, 0.3, 1} x rho in {1.01, 1.5, 4, 12} (task grid) plus mu=5 (Theorem D is claimed for all mu>0)."""
import sys, os, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_taylor as rt

here = os.path.dirname(os.path.abspath(__file__))
allok = True
rows = []
for mu in ('1/50', '3/10', '1', '5'):
    guess = 1.0
    for rho in ('101/100', '3/2', '4', '12'):
        mu_fr, rho_fr = Fraction(mu), Fraction(rho)
        print(f'=== mu={float(mu_fr):g} rho={float(rho_fr):g} ===', flush=True)
        try:
            r = rt.run(mu_fr, rho_fr, guess=guess, out=os.path.join(here, f'case_mu{float(mu_fr):g}_rho{float(rho_fr):g}.json'))
        except Exception as ex:
            import traceback; traceback.print_exc()
            print('  RESULT: EXCEPTION', repr(ex), flush=True); allok = False; continue
        guess = max(r['U0_double'], 0.5)
        vs = rt.verdicts(r)
        bad = [nm for nm, ok in vs if not ok]
        allok &= not bad
        m, x = r['mins'], r['maxs']
        idmax = max(val for kk, val in x.items() if kk.startswith('id:') or kk.startswith('|'))
        sgn = max(val for kk, val in x.items() if kk.startswith('sign:'))
        print(f'  T={r["T"]:.8f} (double {r["T_double"]:.8f}) steps={r["steps"]} Umin={r["Umin"]:.3f} closure={max(r["closure_final"]):.1e}')
        print(f'  lambda_c: <eta>-mu={r["lambda_c: <eta>-mu"]:+.12e}  -<pK/V>={r["lambda_c: -<pK/V> (GL)"]:+.12e}  -<pW/Q>={r["lambda_c: -<pW/Q> (GL)"]:+.12e}  monodromy(double)={r["lam_M_double"]:+.6e}')
        print(f'  <eta0>-mu={r["<eta0>-mu (GL)"]:+.1e} <eta1>-mu={r["<eta1>-mu (GL)"]:+.1e}  quad discrepancy GL/series (pK/V, pW/Q, eta0, eta1): '
              + ' '.join(f'{r["quad_discrepancy(GL vs series), rel"][q]:.0e}' for q in ('pK/V', 'pW/Q', 'eta0', 'eta1')))
        print(f'  min V/(U^2+p^2)={m["V/(U^2+p^2)"]:.3e} (log10 {r["mins_log10"]["V/(U^2+p^2)"][1]:.2f})  min Q/(p^2+pd^2)={m["Q/(p^2+pd^2)"]:.3e} (sign {r["mins_log10"]["Q/(p^2+pd^2)"][0]}, log10 {r["mins_log10"]["Q/(p^2+pd^2)"][1]:.2f})'
              f'  minK_rel={m["K/(|zU|+|ep|+|H|)"]:.2e} minW_rel={m["W/(|zd p|+|ed pd|)"]:.2e} min eta0={m["eta0"]:+.3e}')
        print(f'  identity residual max={idmax:.1e}  sign-pattern max={sgn:+.2e}  kappa in [{m["kappa"]:.3e},{x["kappa"]:.3e}] varpi_max={x["varpi"]:.3e}  eta-mu in [{m["eta-mu"]:+.3f},{x["eta-mu"]:+.3f}]')
        print(f'  <eta p>={r["<eta p>"]:+.6e} <zeta pd>={r["<zeta pd>"]:+.6e} <zeta p>={r["<zeta p>"]:+.6e} <eta pd>={r["<eta pd>"]:+.6e} det={r["det"]:+.4e}')
        print(f'  moment sanity: <x>/k-1={r["<x>/k-1"]:.0e} <p>={r["<p>"]:.0e} <pd^2>/(mu<p^2>)-1={r["<pd^2>/(mu<p^2>)-1"]:.0e} <xp^2>/(A<pd^2>)-1={r["<x p^2>/(A<pd^2>)-1"]:.0e}')
        print(f'  orderings min: ' + ' | '.join(f'{kk.split(":")[0]}:{kk.split(":")[1].strip()}={v:.2e}' for kk, v in r['orderings(min, must be >0)'].items()))
        print(f'  pairs: {r["pairs"]}')
        print(f'  extra: min (y-zeta)/y={m["(y-zeta)/y"]:.3e}  min (lam-p-eta)/lam={m["(lam-p-eta)/lam"]:+.3e}  sign changes p={r["sign changes of p (interior)"]} pd={r["sign changes of pd"]} order_ok={r["phase order crest<pmin<trough<pmax<crest"]}'
              f'  2sigma={2*r["sigma"]:.6e} (gammaE-A-2sigma_np={r["gammaE_minus_A_minus_2sigma_np"]:.1e})')
        print('  CHECKS:', ' '.join(f'{nm}={"ok" if ok else "FAIL"}' for nm, ok in vs))
        print(f'  RESULT: {"ALL OK" if not bad else "FLAGGED: " + ", ".join(bad)}  [{r["seconds"]:.0f}s]', flush=True)
print('GRID:', 'ALL OK' if allok else 'SOME FLAGGED')
