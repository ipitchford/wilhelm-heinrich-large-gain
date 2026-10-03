#!/usr/bin/env python3
"""Write CLAIMS.json for the v3.0 bundle. Statement numbers are resolved from paper/manuscript.aux."""
import json, re, sys
from pathlib import Path
B = Path('/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/v3/wilhelm_heinrich_v3_0')
aux = (B/'paper'/'manuscript.aux').read_text()
lab = {m.group(1): m.group(2) for m in re.finditer(r'\\newlabel\{([^}]+)\}\{\{([^}]*)\}', aux)}

def num(label):
    if label not in lab:
        raise SystemExit(f'label {label} not found in manuscript.aux')
    return lab[label]

P = 'proved'
C = 'proved_conditional'
claims = [
 # --- headline
 dict(id='L1', label='thm:largegain', kind='Theorem', status=P,
      statement='For 0 < mu <= 1 and k >= 1000 every cycle satisfies T >= 0.1833 lambda/mu, has slope p >= 0.819 lambda at the up-crossing of U = 0, has a crest interval of explicit short duration, and satisfies <eta> <= 0.83 mu < mu; it is hyperbolic and attracting.',
      evidence=['evidence/v3/v3_largegain_constants.json', 'evidence/v3/v3_largegain_lemmas.json']),
 dict(id='L2', label='cor:lgunique', kind='Corollary', status=P,
      statement='For 0 < mu <= 1 and k >= 1000 there is exactly one cycle and it attracts every orbit of the positive domain off the one-dimensional stable manifold of the equilibrium.',
      depends=['thm:largegain', 'cor:criterion']),
 dict(id='L3', label='cor:lgreciprocal', kind='Corollary', status=P,
      statement='For mu >= 1 and k >= 1000 mu^3 the same conclusions hold. Only cycles and Floquet multipliers are transported by the reciprocal-rate correspondence; the positive domains do not correspond.',
      depends=['thm:largegain', 'cor:criterion']),
 dict(id='L4', label='thm:mainlarge', kind='Theorem', status=P,
      statement='In original rates: if r0 >= 1000 max(a0,b0)^2 / min(a0,b0) the system has exactly one positive cycle; it is hyperbolic, attracting, and attracts every positive solution off the stable manifold of the equilibrium.',
      depends=['thm:largegain', 'cor:lgunique', 'cor:lgreciprocal']),
 # --- attractor
 dict(id='A1', label='thm:attractor', kind='Theorem', status=P,
      statement='(A1)-(A5), all mu,k > 0: the log-bounded entire orbits form a compact connected global attractor; positive adjoint ratios exist uniquely along each; strict supporting planes hold between any two of its points; the attractor is the graph of the restriction of a convex function over its planar projection, C^1 in the interior with gradient (zeta, -eta); the reduced planar field has divergence eta - mu. (A6)-(A8), k > k_H: E is a source in the attractor; the projection is the closed Jordan domain bounded by the outermost projected cycle; planar multipliers equal m_t; innermost and outermost cycles have lambda_c <= 0.',
      imports=['Birkhoff-Hopf contraction theorem', 'planar Poincare-Bendixson theory', 'stable/unstable manifold theorems', 'Jordan-Schoenflies theorem', 'Whitney extension theorem (C^1)'],
      evidence=['evidence/v3/v3_attractor.json']),
 dict(id='A2', label='thm:twoorbit', kind='Theorem', status=P,
      statement='Strict supporting-plane comparison for any two distinct log-bounded entire solutions.'),
 dict(id='A3', label='cor:criterion', kind='Corollary', status=P,
      statement='Parameter-wise criterion: if every cycle at (mu,k) has lambda_c < 0, the cycle is unique, hyperbolic, attracting, and is the omega-limit set of every positive orbit off the stable manifold of E; if there are two cycles, one has lambda_c >= 0.'),
 dict(id='A4', label='prop:crestaxis', kind='Proposition', status=P,
      statement='Return map on the crest half-axis of the attractor: derivative formula through the one-turn integral of eta - mu; conditional uniqueness criteria. The two monotonicity hypotheses of its last part are NOT proved and fail in floating point at extreme rate ratios.'),
 # --- cycle structure
 dict(id='G1', label='thm:gauges', kind='Theorem', status=P,
      statement='On every cycle Sigma = p^2 - U p\' > 0 and Xi = p\'^2 - p p\'\' > 0: the (U,p) projection is strictly star-shaped about the origin and strictly convex; lambda_c = -<pK/Sigma> = -<pW/Xi>; an explicit state function eta_0 has mean mu and eta - eta_0 = -pK/Sigma.',
      evidence=['evidence/v3/v3_gauges.json', 'evidence/v3/v3_exact_identities.json']),
 dict(id='G2', label='prop:kappa', kind='Proposition', status=P,
      statement='kappa = K/Sigma solves a scalar Riccati equation with explicit coefficients; lambda_c T = -integral of (kappa_up - kappa_down) dU.'),
 dict(id='G3', label='prop:polar', kind='Proposition', status=P, statement='Polar form of the projected flow and of the weak exponent.'),
 dict(id='G4', label='prop:nestedradii', kind='Proposition', status=P, statement='For two nested cycles the radial ratio is monotone on each half-plane.'),
 dict(id='G5', label='prop:monotoneintegrals', kind='Proposition', status=P, statement='<(alpha p + beta p\')(alpha eta + beta zeta)> < 0 for (alpha,beta) != 0.'),
 dict(id='G6', label='prop:holder', kind='Proposition', status=P, statement='<x^(1/3)>^3 < gamma(gamma-1)(gamma-mu) < k on every cycle.'),
 dict(id='G7', label='prop:lvchain', kind='Proposition', status=P, statement='Cyclic Lotka-Volterra form; second-order equation for eta; an explicit pair solving the first Riccati equation with mean mu.'),
 dict(id='G8', label='prop:strongidentity', kind='Proposition', status=P, statement='Exact pairing identity for lambda_c with the strong-stable solution; comparison (supersolution) bound.'),
 dict(id='G9', label='prop:gaussmap', kind='Proposition', status=P, statement='The map from cycle points to (eta, zeta) is injective.'),
 dict(id='W1', label='lem:strongstable', kind='Lemma', status=P, statement='Strong-stable Floquet vector e_s = (a,-b,-c) with a,b,c > 0 and the equations of its ratios; lambda_c < 0 iff the mean of b/a is below A.'),
 dict(id='W2', label='thm:twelve', kind='Theorem', status=P, statement='Each component of the weak adjoint has exactly two simple zeros per period, in a fixed order with the six extrema.',
      evidence=['evidence/v3/v3_weak_adjoint.json']),
 dict(id='W3', label='thm:dualsupport', kind='Theorem', status=P, statement='The plane spanned by the tangent and the strong-stable vector strictly supports the cycle, the equilibrium and nested inner cycles.'),
 dict(id='W4', label='prop:crestmonotone', kind='Proposition', status=P, statement='The crest return map exists, is C^1 and has an entrywise positive derivative on the region {p = 0, k e^U >= y > k}.'),
 dict(id='W5', label='prop:rankone', kind='Proposition', status=P, statement='Rank-one window identity: every periodic test function gives the same single unknown sign.'),
 # --- framework carried over from v2.0
 dict(id='F1', label='prop:multipliers', kind='Proposition', status=P, statement='Positive transverse multipliers m_s < m_t with exact exponents; lambda_c = <eta> - mu.'),
 dict(id='F2', label='prop:factor', kind='Proposition', status=P, statement='Tangent orthogonality and factorisation of the scalar variational operator.'),
 dict(id='F3', label='thm:support', kind='Theorem', status=P, statement='Strict supporting planes with a cycle or the equilibrium as base.'),
 dict(id='F4', label='cor:cones', kind='Corollary', status=P, statement='K > 0, Lambda_E > 0 and exclusion of both closed cones.'),
 dict(id='F5', label='thm:geometry', kind='Theorem', status=P, statement='Injective regular clockwise planar projection; one crest and one trough; nested cycles.'),
 dict(id='F6', label='thm:six', kind='Theorem', status=P, statement='Six nondegenerate extrema in the fixed order X_min, Z_min, Y_min, X_max, Z_max, Y_max.'),
 dict(id='F7', label='prop:crestmap', kind='Proposition', status=P, statement='Entrywise positive crest-return derivative with determinant e^{-AT} < e^{-4 pi}; T > 2 pi / sqrt(mu).'),
 dict(id='F8', label='thm:compact', kind='Theorem', status=P, statement='Uniform compactness of cycle families on compact subsets of the oscillatory domain (state box, separation from E, bounded least periods).'),
 dict(id='F9', label='lem:collar', kind='Lemma', status=P, statement='Global uniqueness in a non-explicit right-hand neighbourhood of each Hopf point.',
      imports=['permanence and competitive Poincare-Bendixson property (Smith 2012)', 'Hopf bifurcation theorem']),
 dict(id='F10', label='thm:neutral', kind='Theorem', status=P, statement='The neutral parameter set is relatively closed; the cycle count is locally constant off it; non-uniqueness at a finite parameter requires a neutral cycle on the fixed-mu path from the collar.'),
 dict(id='F11', label='obs:cor:uniqueset', kind='Corollary', status=P, statement='The set of parameters at which every cycle has lambda_c < 0 is open, its boundary lies in the neutral set, and it contains the Hopf collars and the large-gain region.'),
 # --- wedge
 dict(id='H1', label='lem:tauidentities', kind='Lemma', status=P, statement='Exact identities on every cycle in frequency-one time.'),
 dict(id='H2', label='thm:kurtosis', kind='Theorem', status=P, statement='For a strictly convex potential, <P^4> <= (9/5)<P^2>^2 on every level; the selection function is a bijection with derivative at least 4/5.',
      evidence=['evidence/v3/v3_wedge.json']),
 dict(id='H3', label='thm:wedge', kind='Theorem', status=C,
      statement='Small rate ratio, rho - 1 = r eps^2 with r in a compact interval: exactly one cycle in a fixed energy ball, hyperbolic attracting, with an asymptotic formula for log m_t. Local in phase space; mu_0 not explicit; NOT a global uniqueness statement.',
      condition='property (F), taken from Fenichel\'s theorem; the original sources were not accessible and the theorem numbering is unchecked'),
 # --- negative results
 dict(id='N1', label='thm:failed', kind='Theorem', status=P, statement='The square-root radial contraction certificate fails near Hopf for 0 < mu < 1.'),
 dict(id='N2', label='prop:nodulac', kind='Proposition', status=P, statement='Constraints on Dulac functions for the reduced flow: none regular at E; only one sign is possible.'),
 dict(id='N3', label='prop:nomonomial', kind='Proposition', status=P, statement='No vector supersolution of monomial form in the chemical variables exists on any cycle above Hopf.'),
 dict(id='N4', label='lem:rigidity', kind='Lemma', status=P, statement='Any certificate direction must agree with the true adjoint pair to first order at Hopf.'),
 dict(id='N5', label='neg:prop:signchange', kind='Proposition', status=P, statement='Conditional: a single sign change of the one-turn integral of eta - mu along the crest half-axis implies uniqueness. The hypothesis is unproved.'),
]
for c in claims:
    c['number'] = num(c['label'])
    c['where'] = f"paper/manuscript.pdf, {c['kind']} {c['number']}"

doc = {
 'schema_version': '2.0',
 'version': '3.0',
 'date': '2026-10-03',
 'headline': 'Uniqueness of the cycle is proved for large gain (k >= 1000 max(1, mu^3)) and near each Hopf point; it is NOT proved for all oscillatory parameters.',
 'status_vocabulary': {
   'proved': 'complete proof in the paper; imported results are named',
   'proved_conditional': 'complete proof given a named hypothesis taken from the literature',
   'floating_point_evidence': 'supported by ordinary floating-point computation only',
   'refuted_floating_point': 'contradicted by floating-point computation at stated parameters',
   'open': 'neither proved nor refuted'},
 'claims': claims,
 'not_proved': [
   {'statement': 'Exactly one cycle at every parameter k > mu(1+mu) (the conjecture of Smith 2012).', 'status': 'open',
    'known': 'true for k >= 1000 max(1, mu^3) and in non-explicit Hopf collars'},
   {'statement': 'Every cycle has <eta> < mu.', 'status': 'open', 'note': 'equivalent, up to the degenerate case lambda_c = 0, to uniqueness at each parameter'},
   {'statement': 'The one-turn mean of eta - mu is strictly decreasing in crest amplitude on the attractor.', 'status': 'floating_point_evidence',
    'note': 'would imply uniqueness; held at every parameter tested'},
   {'statement': 'No second cycle and no cycle with lambda_c >= 0 exists in the open region.', 'status': 'floating_point_evidence',
    'evidence': ['evidence/v3/v3_cartography.json', 'evidence/v3/cartography/']},
   {'statement': 'Global uniqueness in the small-rate-ratio wedge.', 'status': 'open'},
   {'statement': 'Historical priority of the results marked new.', 'status': 'open', 'note': 'bounded literature search only'}],
 'refuted_candidates': 'see paper Section on certificates that fail and evidence/v3/v3_negative.json',
 'review': {
   'independent_specialist_review': False,
   'adversarial_internal_review': 'each new theorem was attacked by an independent referee agent with its own numerical tests, and the assembled sections by a second model; records in provenance/research_lanes/*/VERIFY.md and provenance/research_lanes/sol/',
   'v2_0_referee_report': 'review/v2_0_referee_report.md (assesses v2.0 only); response in review/REFEREE_RESPONSE.md'}
}
(B/'CLAIMS.json').write_text(json.dumps(doc, indent=2, ensure_ascii=False) + '\n')
print('CLAIMS.json written with', len(claims), 'claims')
