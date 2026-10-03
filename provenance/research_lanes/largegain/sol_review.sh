#!/bin/zsh
# Cross-model adversarial review of theorem.tex (replies are hypotheses to be verified, not verdicts).
set -a; . /Users/admin/.zprofile >/dev/null 2>&1; set +a
cd /Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/largegain
/Users/admin/.venvs/analysis-py313/bin/python - <<'PY'
import json
tex = open('theorem.tex').read()
prompt = r"""You are a hostile but scrupulously honest referee for a mathematical proof. Your job is to REFUTE it: find any step that is wrong, unjustified, circular, or where a stated numerical constant does not follow. Do not summarise or praise. Work through the document lemma by lemma and re-derive each displayed inequality yourself.

Context: the ODE system is U' = p, p' = -mu p - w, w' = -w + k(e^U - 1), with 0 < mu <= 1, k large. Facts (F1)-(F7) are imported and may be assumed true exactly as stated. Everything else must follow from (F1)-(F7) and the equations.

Pay special attention to:
(1) Lemma 'short' (Sturm argument for T_+ < pi/omega; bound omega >= omega_0(k)).
(2) Lemmas 'i', 'ii', 'rise' (sign conventions; monotonicity of phi; the integral formulas for w and |p|; the use of mu <= 1).
(3) Lemma 'stage1' (the four-duration pigeonhole argument).
(4) Lemmas 'P1' and 'mass' (the energy identity, the bound m <= k(e^U-1-U)/(theta p_k), the final quadratic inequality for M_+).
(5) Lemma 'stage2' and Lemma 'boot' (the continuum bootstrap: is the logic 'T_- > G(T_-) and h quasi-convex => z > z^sharp' airtight? is every monotonicity-in-k and monotonicity-in-mu claim valid?).
(6) Lemma 'crest' (each majorant; the claim that side conditions need only be checked at lambda = K_0).
(7) Lemma 'lyap' (Dini derivative of the max-norm; matching constants; periodicity).
(8) The two corollaries.

Output format: a numbered list of findings. For each: severity (FATAL / GAP / MINOR / COSMETIC), exact location, what is wrong, and a concrete counter-computation or the missing argument. If after genuine effort you find no fatal error in a lemma, say 'no error found' for that lemma in one line. Be concrete; do not speculate without a calculation.

=== theorem.tex ===
""" + tex
payload = {"model": "gpt-5.6-sol", "input": [{"role": "user", "content": prompt}]}
json.dump(payload, open('sol_payload.json', 'w'))
PY
curl -s --max-time 1500 https://api.openai.com/v1/responses \
  -H "Authorization: Bearer $OPENAI_API_KEY" -H "Content-Type: application/json" \
  -d @sol_payload.json > sol_raw.json
jq -r '.output[] | select(.type=="message") | .content[] | select(.type=="output_text") | .text' sol_raw.json > sol_review.md 2>/dev/null
[ -s sol_review.md ] || (echo "NO TEXT; raw head:"; head -c 600 sol_raw.json) > sol_review.md
echo "sol review finished: $(wc -c < sol_review.md) bytes"
