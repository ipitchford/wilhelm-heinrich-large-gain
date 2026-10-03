import sys, json
import numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import find_cycle, kH, cycle_profile
out = []
for mu, rho in [(1.0, 1.2), (1.0, 2.0), (1.0, 8.0), (0.2, 1.5), (0.05, 1.5), (0.05, 3.0)]:
    k = rho * kH(mu)
    c = find_cycle(mu, k, warm=80)
    pr = cycle_profile(mu, k, c, n=4001)
    out.append(dict(mu=mu, rho=rho, k=k, U0=c['U0'], y0=c['y0'], T=c['T'], m_t=c['m_t'], Umin=float(pr['U'].min()),
                    M=c['M'].tolist()))
    print(out[-1])
json.dump(out, open('crest_data.json', 'w'), indent=1)
