#!/bin/bash
# usage: call_sol.sh PROMPT_FILE OUT_JSON [MODEL] [MAXTIME]
set -a; . /Users/admin/.zprofile >/dev/null 2>&1; set +a
PROMPT="$1"; OUT="$2"; MODEL="${3:-gpt-5.6-sol}"; MAXT="${4:-1800}"
/Users/admin/.venvs/analysis-py313/bin/python - "$PROMPT" "$MODEL" > "$OUT.req" <<'PY'
import json,sys
print(json.dumps({"model":sys.argv[2],"input":[{"role":"user","content":open(sys.argv[1]).read()}]}))
PY
curl -s --max-time "$MAXT" https://api.openai.com/v1/responses -H "Authorization: Bearer $OPENAI_API_KEY" -H "Content-Type: application/json" -d @"$OUT.req" > "$OUT"
jq -r '.output[]? | select(.type=="message") | .content[] | select(.type=="output_text") | .text' "$OUT" > "${OUT%.json}.md"
echo "exit=$? bytes=$(wc -c < "$OUT") text=$(wc -c < "${OUT%.json}.md")"; jq -r '.error // empty' "$OUT" | head -5
