# Building the paper

Run the evidence suite first if the generated tables are to be refreshed (`python code/run_all.py`, which ends by writing `paper/generated/` and checking it against the receipts). Then:

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

Three passes resolve the table of contents, theorem and equation numbers, and citations. Only standard LaTeX packages are used: article, geometry, AMS mathematics, mathtools, bm, lmodern, microtype, booktabs, array, enumitem, fancyhdr, xcolor, longtable, graphicx, url and hyperref. The reference PDF was built with pdfTeX from TeX Live 2026 (see `ENVIRONMENT.json`).

Structure:

- `manuscript.tex` inputs `preamble.tex`, the section files in `sections/` in order, and `bibliography.tex`.
- `bibliography.tex` is merged from the fragments in `bib/` (one per section that cites new sources); the merge script is `provenance/research_lanes/frame/assemble.py`. Edit the fragments, not the merged file.
- `generated/` is written by `code/make_tables.py` from the JSON receipts; do not edit it by hand.
- `figures/make_fig_attractor.py` regenerates Figure 1 (an illustration, not evidence).
- `test_<section>.tex` wrappers compile one section on its own; undefined references to other sections are then expected.
- `sections/<section>.report.md` records, for each section, its sources, the corrections applied after review, and the checks run on its formulas.

Inspect the rendered PDF after any rebuild; a successful compiler exit does not establish readable typography.
