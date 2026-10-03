#!/usr/bin/env python3
"""Referee: compact table from the out_v1b_*.json files (battery v1b)."""
import json, glob, math
rows = []
for fn in glob.glob('/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify/out_v1b_*.json'):
    rows.append(json.load(open(fn)))
rows.sort(key=lambda o: (o['mu'], o['rho']))
print("| mu | rho | U_min | T | lambda_c Riccati | lambda_c QR | pairs | #Delta>=0 | signif. | worst (sep>1e-6) | max|O| | min K | min W | cones | P8.1 var | P8.1 par | P8.1 FD |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for o in rows:
    p = [r for r in o['prop81'] if r['var'] == r['var'] and abs(r['pred']) > 1e3 * r['noise_floor']]
    var = max((r['rel'] for r in p), default=float('nan'))
    par = max((abs(r['parallel_defect']) for r in p), default=float('nan'))
    fd = max((r['rel'] for r in o.get('prop81_fd', [])), default=float('nan'))
    print(f"| {o['mu']:g} | {o['rho']:g} | {o['Umin']:.1f} | {o['T']:.3f} | {o['lamc_riccati']:.9f} | {o['lamc_qr']:.9f} | {o['support_pairs']} |"
          f" {o['support_nonneg']} | {o['support_significant']} | {o['support_worst_rel_far']:.1e} | {o['O_rel_max']:.0e} | {o['K_rel_min']:.0e} |"
          f" {o['W_rel_min']:.0e} | {o['n_in_Cplus']}+{o['n_in_Cminus']} | {var:.0e} ({len(p)}) | {par:.0e} | {fd:.0e} ({len(o.get('prop81_fd', []))}) |")
