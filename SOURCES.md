# Source ledger

## Version 3.1 intake and repair update (3 October 2026)

This dated update is retrospective for the supplied v3.0 research and prospective
for the slow-graph repair. It does not backdate the publication agent's reading.

| Primary source | Exact text inspected | Consequence |
| --- | --- | --- |
| Smith (2012), https://arxiv.org/pdf/1109.5204 | Theorems 1.1 and 1.4 and the positive-interior/nonnegative attractor distinction | Supports inherited finiteness/convergence; motivates the corrected attractor wording |
| Boros–Hofbauer (2023), https://arxiv.org/pdf/2202.11034 | Section 4, Theorem 13 and adjacent model definitions | Confirms the original-model cycle-count context; the homogenized extension is not the same equation |
| Baigent (2019), DOI 10.3934/dcdsb.2018288, https://discovery.ucl.ac.uk/id/eprint/10049790/ and author PDF https://www.ucl.ac.uk/~ucess29/resources/map2convexDCDS.pdf | Section 4 standing assumptions and Theorem 6.1 | A convex carrying-simplex antecedent; its face-preserving Kolmogorov-map hypothesis fails directly in the physical coordinates here |
| Mierczynski (2020), DOI 10.1017/etds.2018.85, https://arxiv.org/pdf/1801.01032 | Sections 2–3, including forward/backward coordinate-face invariance | A regularity antecedent, not a direct application to the oscillator flow |
| Eldering (2013), DOI 10.2991/978-94-6239-003-4, author pre-proof dated 1 May 2018, https://www.jaapeldering.nl/files/NHIM-noncompact-book.pdf | Theorem 3.2, Remark 3.3, and relevant definitions; targeted checker also inspected Section 4.2 | Load-bearing persistence result for the new local slow-graph lemma. Its hypotheses are mapped explicitly in the manuscript |
| Eldering author errata, https://www.jaapeldering.nl/files/errata.pdf, 1 July 2021 | Both pages | Boundary extension not used; no correction withdrawing the invoked theorem |

These are specified full-text portions, not claims to have re-read every source
cover to cover. Targeted searches included the supplied Baigent/Mierczynski
titles, their DOIs, “Li Bingxi 1981 uniqueness third order differential equation”
and DOI 10.1016/0362-546X(81)90065-1. The Li publisher route did not return usable
full text. Whether Li's hypotheses cover this model remains unresolved.
The original Fenichel/Jones texts were not obtained. They remain historical
references, but the v3.1 proof of property (F) now relies on the inspected
Eldering theorem, not on unchecked theorem numbering.

The earlier ledger below describes the input package's source work. Its
Fenichel/Jones status is superseded only by the precise local repair above;
other access and priority limitations remain.

Search and inspection date: **3 October 2026**. This is a bounded check of primary sources, not a systematic priority certificate. The working record of the search, with queries, citation-graph dumps and fetched open texts, is in `provenance/research_lanes/literature/`.

Access levels: **FT** full text read; **AB** abstract or publisher record only; **SN** search snippet or publisher description only; **Q** statement taken from a later paper that quotes it; **META** bibliographic metadata verified (Crossref or arXiv) but text not accessed in this work.

## 1. The model and its established dynamics

| Source | Access | Use in the paper |
|---|---|---|
| Wilhelm & Heinrich (1995), *J. Math. Chem.* 17, 1–14. doi:10.1007/BF01165134 | AB | Model provenance. No theorem imported. |
| Wilhelm & Heinrich (1996), *J. Math. Chem.* 19, 111–130. doi:10.1007/BF01165179 | AB | Local Hopf analysis is acknowledged; the coefficient needed is recomputed in Appendix A. |
| Smith, H. L. (2012), *J. Math. Chem.* 50, 989–995. doi:10.1007/s10910-011-9946-9; arXiv:1109.5204 | FT (preprint) | Permanence, competitive Poincaré–Bendixson property, finiteness, convergence to a cycle off the stable manifold (Theorems 1.1, 1.3, 1.4). Used in Lemma 9.1 (Step 2) and as an alternative route in Corollary 10.2. The cycle count is stated open there. |
| Boros & Hofbauer (2023), *Nonlinear Anal. RWA* 72, 103839. doi:10.1016/j.nonrwa.2023.103839; arXiv:2202.11034 | FT (preprint) | Theorem 13: summary of the established results, including critical global stability; cycle count retained as open. |
| Ruan (2006), *Delay Differential Equations and Applications*, 477–517. doi:10.1007/1-4020-3647-7_11 | FT (author copy) | Identification of the system with the linear-chain form of the distributed-delay logistic equation; only local Hopf results exist for the strong kernel. |

No later work citing these papers was found that addresses the number of cycles (citation graphs via OpenAlex and Semantic Scholar, arXiv listings, keyword searches).

## 2. The analogous delay equation (Jones' conjecture)

| Source | Access | Use |
|---|---|---|
| Xie (1993), *J. Differential Equations* 103, 350–374. doi:10.1006/jdeq.1993.1054 | Q | Large-parameter stability of every slowly oscillating periodic solution; "all stable implies unique". Positioning only. |
| Jaquette, Lessard & Mischaikow (2017), *J. Differential Equations* 263, 7263–7286. doi:10.1016/j.jde.2017.08.018; arXiv:1705.02432 | FT (preprint) | Computer-assisted middle range. Positioning and programme. |
| van den Berg & Jaquette (2018), *J. Differential Equations* 264, 7412–7462. doi:10.1016/j.jde.2018.02.018; arXiv:1704.00029 | FT (preprint, parts) | Explicit neighbourhood of the Hopf point. Positioning and programme. |
| Jaquette (2019), *J. Differential Equations* 266, 3818–3859. doi:10.1016/j.jde.2018.09.017; arXiv:1801.09806 | FT (preprint) | Completion of the proof. Positioning. |

## 3. Structure of competitive and cyclic feedback systems

| Source | Access | Use |
|---|---|---|
| Mallet-Paret & Smith (1990), *J. Dynam. Differential Equations* 2, 367–421. doi:10.1007/BF01054041 | AB | Positioning: planar projection and sign structure for monotone cyclic feedback systems. A theorem-level priority comparison was not possible without the full text. |
| Hirsch (1982), *SIAM J. Math. Anal.* 13, 167–179. doi:10.1137/0513013 | AB | Positioning: limit sets of competitive systems project injectively. |
| Walther (1995), *Mem. Amer. Math. Soc.* 113(544). doi:10.1090/memo/0544 | SN | Positioning: two-dimensional graph attractor for a scalar delay equation, compatible with several nested cycles. |
| Baigent (2013), *Nonlinearity* 26, 1001–1029. doi:10.1088/0951-7715/26/4/1001 | AB | Positioning: convexity of carrying simplices by a related mechanism. |
| Smith, H. L. (1986), *Nonlinear Anal.* 10, 1037–1052. doi:10.1016/0362-546X(86)90087-8 | Q | Monotone concave fixed-point theorem whose hypotheses were tested and fail (Section 12). |

## 4. Planar uniqueness criteria tested in Section 12

| Source | Access | Use |
|---|---|---|
| Massera (1954), *Boll. Un. Mat. Ital.* (3) 9, 367–369 | not obtained | Named as the origin of the radial-scaling criterion; not in Crossref. |
| Zhang (1986), *Applicable Analysis* 23, 63–76. doi:10.1080/00036818608839631 | META | Zhang-type criterion; hypotheses stated from secondary accounts. |
| Sabatini & Villari (2010), *Le Matematiche* 65(2); arXiv:1101.2761 | FT (preprint) | Radial angular monotonicity criterion. |
| Sabatini (2010), *Chaos, Solitons & Fractals* 43, 25–30. doi:10.1016/j.chaos.2010.07.002; arXiv:1006.2020 | FT (preprint) | Star-shapedness criterion for inseparable equations. |
| Gasull & Giacomini (2021), *J. Differential Equations* 305, 347–367. doi:10.1016/j.jde.2021.10.011; arXiv:2101.03874 | FT (preprint) | Extended Bendixson–Dulac criterion; the curvature function. |
| Farcot & Gouzé (2010), *Int. J. Systems Science* 41, 119–130. doi:10.1080/00207720903144552; arXiv:0903.1969 | FT (preprint) | Source of the quoted statement of H. L. Smith's concave-map theorem; uniqueness for piecewise-affine loops. |
| Duff (1953), *Ann. of Math.* 57, 15–31. doi:10.2307/1969724 | META | Rotated vector fields; the hypothesis fails for the reduced field. |

## 5. Classical results quoted in standard form

These are textbook results. Their metadata were verified through Crossref; the texts were not re-read for this paper, and each is quoted in the form stated at its point of use.

| Source | Use |
|---|---|
| Birkhoff (1957), *Trans. Amer. Math. Soc.* 85, 219–227; Bushell (1973), *Arch. Rational Mech. Anal.* 52, 330–338; Seneta (1981), *Non-negative Matrices and Markov Chains* | Contraction of Hilbert's projective metric by positive matrices (Section 7, item I1). |
| Hartman (2002), *Ordinary Differential Equations*, 2nd ed., SIAM | Planar Poincaré–Bendixson theory; invariant manifolds; orbital stability (Section 7, items I2, I3). |
| Thomassen (1992), *Amer. Math. Monthly* 99, 116–130 | Jordan–Schoenflies theorem (Section 7, item I4). |
| Whitney (1934), *Trans. Amer. Math. Soc.* 36, 63–89 | C^1 extension theorem, used only at the outermost cycle (Section 7, item I5). |
| Fenichel (1979), *J. Differential Equations* 31, 53–98. doi:10.1016/0022-0396(79)90152-9; Jones (1995), in *Dynamical Systems*, Lecture Notes in Mathematics, 44–118. doi:10.1007/BFb0095239 | Persistence of normally hyperbolic slow manifolds. Theorem 11.9 takes one consequence as an explicit hypothesis. The originals returned access errors; the theorem numbering quoted in Section 11 is unchecked. |

## 6. Access gaps that matter

1. Mallet-Paret & Smith (1990), full text: needed for a theorem-by-theorem priority comparison.
2. Li, Bingxi (1981), *Nonlinear Analysis* 5, 13–19 (doi:10.1016/0362-546X(81)90065-1) and 5, 931–958: a uniqueness theorem for a third-order neuron model whose hypotheses were not obtained. It is not used. Whether it covers this system is unknown.
3. Wilhelm & Heinrich (1995, 1996), full texts.
4. Fenichel (1979) and Jones (1995), full texts.
5. Xie (1993) and the 1991 thesis, full texts.

A negative bounded search is not proof of novelty. No third-party article is redistributed in this bundle; the fetched open-access texts stay outside it.

## 7. Working inputs

The v2.0 bundle is preserved byte-for-byte in `provenance/v2_0_bundle/`, together with its own provenance (v1.2 attachments and the earlier v0.3 review). The referee report on v2.0 is in `review/v2_0_referee_report.md`.
