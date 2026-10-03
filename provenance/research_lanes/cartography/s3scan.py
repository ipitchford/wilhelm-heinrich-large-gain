import sys, json
sys.path.insert(0,'.')
import numpy as np, cart as C
from attractor import sample_attractor
out=[]
print("S3a at the cycle crest: Qn = normalised d/dU[(mu-zeta)/(k(e^U-1))] along the attractor section (>0 = pass)")
for mu in (0.3,0.5,0.7,0.85,1.0,1.2,1.5,2.0,3.0):
    line=f"mu={mu:<5}"
    for rho in (2.5,3.5,4.0,4.5,5.0,6.0,8.0,10.0):
        try:
            S = sample_attractor(mu, rho, n_theta=16, n_out=0)
        except Exception as e:
            line += "   err  "; continue
        TT=S['turns']; k=S['k']; Us=S['cycle']['U0']
        conv = (np.abs(TT['yU0']/TT['zeta0']-1)<1e-6)
        U=TT['U0']; z=TT['zeta0']; zU=TT['zU0']
        Q=(z-mu)*np.exp(U)-zU*np.expm1(U); Qn=Q/(np.abs(z-mu)*np.exp(U)+np.abs(zU*np.expm1(U)))
        near = conv & (np.abs(U/Us-1)<1e-7)
        r = (1+mu)*S['cycle']['T']/abs(S['profile']['log_mt'])-1
        qc = float(np.median(Qn[near])) if near.any() else float('nan')
        qmin = float(Qn[conv].min())
        out.append(dict(mu=mu,rho=rho,Qn_at_cycle=qc,Qn_min_section=qmin,r=r))
        line += f" {qc:+.3f}({r:4.1f})"
    print(line, flush=True)
print("columns rho = 2.5,3.5,4,4.5,5,6,8,10 ; value = Qn at cycle crest, (r) = smoothness exponent log m_s/log m_t")
json.dump(out, open('s3a_cycle_scan.json','w'), indent=1)
