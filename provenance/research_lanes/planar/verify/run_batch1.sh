#!/bin/bash
# referee batch: Theorem A battery on the 12 requested parameter pairs (+1 extra with mu>1)
cd /Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/planar/verify
PY=/Users/admin/.venvs/analysis-py313/bin/python
for mu in 1 0.3 0.02; do
  for rho in 1.5 4 12 1.01; do
    echo "##### $mu $rho"
    $PY v1_attractor.py $mu $rho 2>&1
  done
done
echo "##### 3 1.5"
$PY v1_attractor.py 3 1.5 2>&1
echo ALLDONE
