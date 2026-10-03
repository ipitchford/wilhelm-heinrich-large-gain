# Building the paper

Run the evidence suite before regenerating numerical tables. The PDF uses only standard LaTeX packages: article, geometry, AMS mathematics, mathtools, lmodern, microtype, booktabs, enumitem, fancyhdr, xcolor and hyperref.

```sh
python code/make_tables.py
cd paper
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

The second pass resolves theorem, equation and bibliography references. The reference PDF was built with pdfTeX from TeX Live 2025. No external figures, bibliographic service, downloaded fonts or shell escape are required. Render the PDF after rebuilding and inspect every page; a successful compiler exit does not establish readable typography.

Generated numerical tables and macros are committed with the source. Their source-receipt hashes appear in `evidence/manuscript_summary.json`. A rebuild or rerun mutates generated files; do it in a working copy if preserving the original manifest matters.
