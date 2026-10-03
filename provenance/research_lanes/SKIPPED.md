# Research record: what was copied and what was skipped

This directory is a copy of the working record of the research phase that preceded the paper: the task brief
(`BRIEF.md`), the digest of results and referee verdicts (`attack_results.json`), the shared numerical toolkit
(`whlib.py`), the consistency checks of the framework sections (`framework_checks/`), the cross-model reviews
(`sol/`: prompts, raw responses and extracted text), and one directory per research line with its notes
(`NOTES.md`), the adversarial referee report (`VERIFY.md`, with the referee scripts in `verify/`), proofs in LaTeX,
scripts and text outputs.

The files are copies of what was executed. They contain absolute paths of the machine on which the research was
done and are not meant to be run in place. The scripts that the paper relies on are re-executed by
`code/run_all.py` from copies under `code/v3/lib/`, each of which names its original here in its first lines.
Statements in the notes are working statements; where a referee report corrects them, the report and the paper
take precedence.

Copied: 476 files, 6.46 MB. Skipped: 93 files, 128.38 MB.

Rules applied: copy notes, referee reports, LaTeX proofs, scripts and text outputs (extensions .md .py .tex .txt .out
.json .jsonl .csv .sh .jl .log .req .html); skip compiled Python caches, scratch directories, superseded runs, raw web
material, binary intermediate data, compiled PDFs of copied LaTeX files, images, empty files, and any file over 20 MB.

## Skipped files

### Binary intermediate data (numpy archive), not a script or text output (14 files, 58.01 MB)

- `cartography/hunt_results_pairs.npz` (3.49 MB)
- `cartography/portrait_mu0.01_rho3.npz` (0.63 MB)
- `cartography/portrait_mu0.1_rho1.5.npz` (0.77 MB)
- `cartography/portrait_mu1_rho2.npz` (0.77 MB)
- `dulac/data_mu0.05_rho1.02.npz` (4.69 MB)
- `dulac/data_mu0.05_rho1.5.npz` (5.91 MB)
- `dulac/data_mu0.05_rho5.0.npz` (9.62 MB)
- `dulac/data_mu0.2_rho1.05.npz` (3.90 MB)
- `dulac/data_mu0.2_rho2.0.npz` (4.87 MB)
- `dulac/data_mu0.5_rho1.2.npz` (3.89 MB)
- `dulac/data_mu0.5_rho3.0.npz` (5.20 MB)
- `dulac/data_mu1.0_rho1.02.npz` (4.11 MB)
- `dulac/data_mu1.0_rho1.5.npz` (3.90 MB)
- `dulac/data_mu1.0_rho5.0.npz` (6.28 MB)

### Single file over 20 MB (1 files, 49.21 MB)

- `supersolution/cycles.pkl` (49.21 MB)

### Raw web material (third-party papers and API responses; 15.4 MB in total) (54 files, 15.43 MB)

Papers and API responses downloaded during the literature search. They are third-party material and are not
redistributed; the bibliographic record of what was consulted is in `literature/NOTES.md`.

- `literature/raw/0903.1969.pdf` (0.32 MB)
- `literature/raw/0903.1969.txt` (70.6 kB)
- `literature/raw/0909.4107.pdf` (0.44 MB)
- `literature/raw/0909.4107.txt` (0.11 MB)
- `literature/raw/1006.2020.pdf` (0.43 MB)
- `literature/raw/1006.2020.txt` (35.1 kB)
- `literature/raw/1011.1785.pdf` (0.47 MB)
- `literature/raw/1011.1785.txt` (58.7 kB)
- `literature/raw/1101.2761.pdf` (0.29 MB)
- `literature/raw/1101.2761.txt` (32.6 kB)
- `literature/raw/1305.3402.pdf` (0.21 MB)
- `literature/raw/1305.3402.txt` (51.2 kB)
- `literature/raw/1704.00029.pdf` (0.72 MB)
- `literature/raw/1704.00029.txt` (0.17 MB)
- `literature/raw/1705.02432.pdf` (0.60 MB)
- `literature/raw/1705.02432.txt` (82.0 kB)
- `literature/raw/1801.09806.pdf` (0.59 MB)
- `literature/raw/1801.09806.txt` (0.16 MB)
- `literature/raw/1906.03528.pdf` (1.65 MB)
- `literature/raw/1906.03528.txt` (0.17 MB)
- `literature/raw/2101.03874.pdf` (0.28 MB)
- `literature/raw/2101.03874.txt` (63.1 kB)
- `literature/raw/2407.00461.pdf` (1.45 MB)
- `literature/raw/2407.00461.txt` (60.6 kB)
- `literature/raw/2503.09155.pdf` (0.43 MB)
- `literature/raw/2503.09155.txt` (77.9 kB)
- `literature/raw/abstracts1.txt` (6.6 kB)
- `literature/raw/bh2023.json` (15.2 kB)
- `literature/raw/bh2023.pdf` (0.73 MB)
- `literature/raw/bh2023.txt` (64.4 kB)
- `literature/raw/cites_W1980509781.json` (25.8 kB)
- `literature/raw/cites_W2030734556.json` (44.8 kB)
- `literature/raw/cites_W2047410993.json` (0.16 MB)
- `literature/raw/cites_W4317940786.json` (20.9 kB)
- `literature/raw/dai1979.pdf` (1.21 MB)
- `literature/raw/hastings1977.pdf` (1.2 kB)
- `literature/raw/match2025.html` (8.1 kB)
- `literature/raw/match2025.txt` (56.4 kB)
- `literature/raw/math_0307372.pdf` (0.16 MB)
- `literature/raw/math_0307372.txt` (30.9 kB)
- `literature/raw/math_0409387.pdf` (0.16 MB)
- `literature/raw/math_0409387.txt` (23.1 kB)
- `literature/raw/rasmith1987.pdf` (1.21 MB)
- `literature/raw/ruan_nato.pdf` (0.97 MB)
- `literature/raw/ruan_nato.txt` (0.12 MB)
- `literature/raw/s2_bh_arxiv_cites.json` (0.2 kB)
- `literature/raw/s2_bh_cites.json` (3.8 kB)
- `literature/raw/s2_smith_cites.json` (3.8 kB)
- `literature/raw/sanchez2006.pdf` (1.21 MB)
- `literature/raw/smith2012.json` (9.9 kB)
- `literature/raw/smith2012.pdf` (0.12 MB)
- `literature/raw/smith2012.txt` (20.2 kB)
- `literature/raw/wh95.json` (10.5 kB)
- `literature/raw/wh96.json` (13.2 kB)

### Superseded run, replaced by the files that were copied (3 files, 4.08 MB)

- `cartography/superseded/hunt_results_v1_before_event_filter.jsonl` (0.51 MB)
- `cartography/superseded/hunt_results_v1_pairs.npz` (3.53 MB)
- `cartography/superseded/hunt_stdout_v1.txt` (47.9 kB)

### Compiled PDF of a .tex file that was copied (3 files, 1.16 MB)

- `largegain/theorem.pdf` (0.43 MB)
- `monotone/ds_theorem.pdf` (0.29 MB)
- `planar/attractor_theorem.pdf` (0.44 MB)

### Image output (3 files, 0.37 MB)

- `cartography/portrait_mu0.01_rho3.png` (0.11 MB)
- `cartography/portrait_mu0.1_rho1.5.png` (0.13 MB)
- `cartography/portrait_mu1_rho2.png` (0.13 MB)

### Scratch directory (debugging scripts) (10 files, 98.8 kB)

- `cartography/scratch/branch_table_partial.csv` (59.3 kB)
- `cartography/scratch/dbg_attr.py` (0.9 kB)
- `cartography/scratch/dbg_attr2.py` (0.8 kB)
- `cartography/scratch/dbg_r.py` (0.9 kB)
- `cartography/scratch/dbg_s3.py` (1.1 kB)
- `cartography/scratch/explore1.py` (1.6 kB)
- `cartography/scratch/grid_log_partial.txt` (8.9 kB)
- `cartography/scratch/grid_table_partial.csv` (23.6 kB)
- `cartography/scratch/prof2.py` (0.6 kB)
- `cartography/scratch/prof_hunt.py` (1.0 kB)

### Binary intermediate data (pickle), not a script or text output (5 files, 37.1 kB)

- `supersolution/famA2_U_p_ylin.pkl` (7.0 kB)
- `supersolution/famA2_U_y_P.pkl` (7.0 kB)
- `supersolution/famA2_U_y_P_LamE.pkl` (7.9 kB)
- `supersolution/famA2_U_y_P_Pstar.pkl` (7.9 kB)
- `supersolution/verify/vcycles.pkl` (7.4 kB)

## Notes

- `cartography/hunt_results_pairs.npz` is skipped here as binary data but is shipped with the cartography tables under
  `evidence/v3/cartography/`.
- The 49 MB file `supersolution/cycles.pkl` and the ten `dulac/data_*.npz` files are caches of computed cycles and of
  samples of the unstable manifold; the scripts that produce them are copied (`supersolution/base.py`, `dulac/gen_data.py`).
