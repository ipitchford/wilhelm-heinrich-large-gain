# appB_reproduce.tex — report of the evidence engineer (3 October 2026)

Label `app:reproduce`. Compile test: `paper/test_appB_reproduce.tex` (no LaTeX error, no overfull box).

## Content and sources

| Paragraph | Source |
|---|---|
| Layout of the bundle | the directories as built: `code/`, `code/v3/`, `code/v3/lib/` (provenance header in the first three lines of every library copy; checked by `check_consistency.py`), `evidence/`, `logs/`, `paper/generated/`, `provenance/research_lanes/` |
| Commands | `code/run_all.py` (options `--regenerate`, `--quick`), `code/check_consistency.py`, `code/verify_manifest.py` (unchanged from v2.0; needs the `MANIFEST.sha256` written when the release is assembled), `code/v3/v3_largegain_lemmas.py --small-mu-supplement` |
| Environment and cost | `ENVIRONMENT.json`, `requirements.txt`, `evidence/run_receipt.json` (macros `\VthreeRuntimeMinutes`, `\VthreePeakMemoryMb`); v2.0 environment from `provenance/v2_0_bundle/ENVIRONMENT.json` |
| Verification boundary | v2.0 Appendix B, extended to the three kinds of check of `num:sec:classes` |

## Notes for the editor

- The version numbers in the paragraph "Environment and cost" are typed from `ENVIRONMENT.json`; `check_consistency.py`
  verifies that `ENVIRONMENT.json` agrees with the packages recorded in the run receipt, not the typed text.
- `verify_manifest.py` fails until `MANIFEST.sha256` exists; the manifest is not written by the evidence engineer.
  A rerun of `run_all.py` after the manifest is written invalidates the digests of receipts, logs and generated files.
- `run_all.py` leaves `__pycache__` directories under `code/` (Python bytecode and the numba cache of
  `code/v3/lib/cart.py`). They were removed after the recorded run and should be excluded from the manifest.
- The appendix uses a `verbatim` block, which is safe inside `\input`.
