#!/bin/bash
cd /Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify
PY=/Users/admin/.venvs/analysis-py313/bin/python
for mu in 1 0.3 0.02; do
  for rho in 1.5 4 12 1.01; do
    $PY v6_extras.py $mu $rho 2>&1
  done
done
echo ALLDONE
