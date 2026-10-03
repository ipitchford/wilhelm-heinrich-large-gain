import sys, numpy as np
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work')
from whlib import *
from job1_section import sections
mu,rho=1.0,1.1
k=rho*kH(mu); ed=equilibrium_data(mu,k); cyc=find_cycle(mu,k)
tend=(np.log(1e6)+10)/ed['lam_u'].real
S1=[];S2=[]
for j,th in enumerate(np.linspace(0,2*np.pi,7,endpoint=False)):
    v=np.real(np.exp(1j*th)*ed['vec_u']); z0=np.array([0,0,k])+1e-6*v/np.linalg.norm(v)
    a,b=sections(mu,k,z0,[ed['etaE'],ed['zetaE']],tend)
    S1.append(np.c_[a,np.full(len(a),j)]); S2.append(np.c_[b,np.full(len(b),j)])
S2=np.vstack(S2); S2=S2[np.argsort(S2[:,1])]; S2=S2[S2[:,1]>1e-3]
dp=np.diff(S2[:,1]); dy=np.diff(S2[:,2]); sl=dy/dp
lo=(sl+S2[:-1,3])/S2[:-1,3]; hi=(-S2[1:,3]-sl)/S2[1:,3]
for i in np.argsort(lo)[:6]:
    print('lo',i,len(S2),'p',S2[i,1],S2[i+1,1],'dp',dp[i],'eta',S2[i,3],S2[i+1,3],'slope',sl[i],'U',S2[i,0],S2[i+1,0],'orb',S2[i,5],S2[i+1,5])
for i in np.argsort(hi)[:3]:
    print('hi',i,'p',S2[i,1],S2[i+1,1],'dp',dp[i],'eta',S2[i,3],S2[i+1,3],'slope',sl[i])
print('max p on section', S2[-5:,1])
S1=np.vstack(S1); S1=S1[np.argsort(S1[:,0])]; S1=S1[(S1[:,0]>1e-3)&(S1[:,0]<cyc['U0']-1e-5)]
dU=np.diff(S1[:,0]); dy=np.diff(S1[:,2]); sl=dy/dU
lo=(sl-S1[:-1,4])/S1[:-1,4]
for i in np.argsort(lo)[:5]:
    print('S1 lo',i,len(S1),'U',S1[i,0],'dU',dU[i],'zeta',S1[i,4],S1[i+1,4],'slope',sl[i],'rem/dU^2',(dy[i]-S1[i,4]*dU[i])/dU[i]**2, 'bracket',(S1[i+1,4]-S1[i,4]))
