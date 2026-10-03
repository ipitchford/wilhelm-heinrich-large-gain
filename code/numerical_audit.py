#!/usr/bin/env python3
"""Independent numerical consistency checks of fresh BVPs and new geometry.
The ratios, normalised linear adjoint, and state/return derivative are integrated
separately. These are floating-point diagnostics, not enclosures or exclusions.
"""
from __future__ import annotations
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.interpolate import PPoly
from scipy.optimize import brentq
from cycles import ROOT,strict_write,equations

RTOL=3e-11;ATOL=2e-13

def valid_phase_order(seq):
    """Shared scientific predicate, also called by the negative controls."""
    return (len(seq)==6 and bool(np.all(np.isfinite(seq))) and seq[0]==0
            and all(b>a for a,b in zip(seq,seq[1:])) and seq[-1]<1)

def valid_winding(value):
    return bool(np.isfinite(value) and abs(value+1)<1e-6)

def radial_certificate_pass(ratios):
    values=np.asarray(ratios,dtype=float)
    return bool(values.shape==(2,) and np.all(np.isfinite(values))
                and np.all(values>0) and np.all(values<1))

def load(path):
    q=np.load(path);pp=PPoly(q['coefficients'],q['breaks'],extrapolate=False)
    def state(s):return pp(np.mod(s,1.)).T
    return state,float(q['mu']),float(q['lam']),float(q['period_tau'])

def roots(fun,n=4001):
    grid=np.linspace(1e-8,1-1e-8,n);vals=np.asarray(fun(grid));out=[]
    for i in np.flatnonzero(vals[:-1]*vals[1:]<0):
        v=brentq(fun,grid[i],grid[i+1],xtol=5e-15,rtol=1e-14)
        if not out or abs(v-out[-1])>1e-8:out.append(float(v))
    return out

def avg(v,t):return float(simpson(v,x=t))
def rel(a,b):return abs(a-b)/max(abs(a),abs(b),1e-12)

def audit(row,method='DOP853',save_profiles=True):
    start=time.time();state,mu,lam,Ttau=load(ROOT/row['mesh']);k=mu*lam;Tf=Ttau/np.sqrt(mu);A=1+mu
    lower=max(0.,mu-1.)
    se=brentq(lambda s:s*(s+1)*(s+1-mu)-k,lower,max(mu+2,2+k),xtol=1e-14)
    ze=se*(se+1-mu)
    def rhs(q,v):
        eta,zeta=v;x=k*np.exp(state(q)[0])
        return Tf*np.array([zeta+(mu-1)*eta-eta*eta,x-(1+eta)*zeta])
    training=solve_ivp(rhs,(0,2),[se,ze],rtol=RTOL,atol=ATOL,method=method)
    if not training.success:raise RuntimeError('projective training failed: '+training.message)
    grid=np.linspace(0,1,10001)
    ratio=solve_ivp(rhs,(0,1),training.y[:,-1],t_eval=grid,dense_output=True,rtol=RTOL,atol=ATOL,method=method)
    if not ratio.success:raise RuntimeError('projective period failed: '+ratio.message)
    eta,zeta=ratio.y;U,V,z=state(grid);x=k*np.exp(U);p=np.sqrt(mu)*V;y=k*z;fp=k-mu*p-y
    ep=zeta+(mu-1)*eta-eta**2;zp=x-(1+eta)*zeta
    O=zeta*p-eta*fp-(x-y);W=zp*p-ep*fp
    lc=avg(eta,grid)-mu;ms=np.exp(-(1+avg(eta,grid))*Tf);mt=np.exp(lc*Tf)
    # Independent normalised positive linear-adjoint propagation.
    def lin_rhs(q,v):
        h=v[:3];xq=k*np.exp(state(q)[0]);bh=np.array([xq*h[2],h[0]+mu*h[1],h[1]+h[2]]);growth=np.sum(bh)
        return np.r_[Tf*(bh-growth*h),Tf*growth]
    tr=solve_ivp(lin_rhs,(0,2),[1/3,1/3,1/3,0],rtol=RTOL,atol=ATOL,method=method)
    if not tr.success:raise RuntimeError('linear adjoint training failed')
    li=solve_ivp(lin_rhs,(0,1),np.r_[tr.y[:3,-1]/np.sum(tr.y[:3,-1]),0.],rtol=RTOL,atol=ATOL,method=method)
    if not li.success:raise RuntimeError('linear adjoint measurement failed')
    lc_lin=float(li.y[3,-1]/Tf-A)
    # Actual nonlinear IVP, with two initial section variations and event correction.
    initial=state(0);K0=np.array([[1.,0],[0,0],[0,1.]])
    def aug_rhs(q,aug):
        w=aug[:3];K=aug[3:].reshape(3,2);UU,VV,zz=w;rt=np.sqrt(mu)
        J=np.array([[0,1,0],[0,-rt,-lam],[np.exp(UU)/rt,0,-1/rt]])
        return np.r_[Ttau*equations(mu,w,lam),(Ttau*J@K).ravel()]
    def event(q,aug):return aug[1]
    event.direction=-1;event.terminal=False
    iv=solve_ivp(aug_rhs,(0,1.01),np.r_[initial,K0.ravel()],events=event,dense_output=True,rtol=RTOL,atol=ATOL,method=method)
    ev=[t for t in iv.t_events[0] if t>1e-6] if iv.success else []
    M=None;return_ok=len(ev)==1
    if return_ok:
        tret=ev[0];out=iv.sol(tret);K=out[3:].reshape(3,2);F=equations(mu,out[:3],lam)
        K-=F[:,None]*(K[1,:]/F[1]);M=K[[0,2],:]
        eig=np.linalg.eigvals(M);Mt=float(np.max(eig.real))
        closure=float(np.max(np.abs(out[:3]-initial)))
        logdet=-A*Ttau*tret/np.sqrt(mu)
    else:Mt=closure=logdet=float('nan')
    # Six simple extrema in one fixed order, starting from X maximum.
    rq=roots(lambda q:state(q)[1]);rz=roots(lambda q:np.exp(state(q)[0])-state(q)[2]);rp=roots(lambda q:k-mu*np.sqrt(mu)*state(q)[1]-k*state(q)[2])
    counts=(len(rq),len(rz),len(rp))
    phase_order=False;phases={}
    if counts==(1,2,2):
        phases={'X_max':0.,'Z_max':rz[0],'Y_max':rp[0],'X_min':rq[0],'Z_min':rz[1],'Y_min':rp[1]}
        seq=list(phases.values());phase_order=valid_phase_order(seq)
    # Supporting-plane gaps for a fixed nonzero phase shift.
    shift=.231;U2,V2,z2=state(grid+shift);du=U2-U;dp=np.sqrt(mu)*V2-p;dy=k*z2-y
    gap=dy-zeta*du+eta*dp
    forcing=x*(np.expm1(du)-du)
    support_balance=rel(avg((1+eta)*gap,grid),avg(forcing,grid))
    Kgap=zeta*U-eta*p-y+k
    Egap=-(ze*U-se*p-y+k)
    winding=float((np.unwrap(np.arctan2(p,U))[-1]-np.unwrap(np.arctan2(p,U))[0])/(2*np.pi))
    xmin=k*np.exp(state(rq[0])[0]) if len(rq)==1 else float(np.min(x));xmax=k*np.exp(state(0)[0])
    splus=(mu-1+np.sqrt((mu-1)**2+4*xmax))/2
    zminus=xmin/(1+splus)
    sminus=(mu-1+np.sqrt((mu-1)**2+4*zminus))/2
    zplus=xmax/(1+sminus)
    scale=max(1.,float(np.max(x)),float(np.max(y)))
    moments={'mean_x':rel(avg(x,grid),k),'mean_y':rel(avg(y,grid),k),'p_energy':rel(avg(fp**2,grid),mu*avg(p**2,grid)),'gain_identity':rel(avg(eta*(eta+1)*(eta+1-mu),grid),k),'curvature_identity':rel(avg((1+eta)*W,grid),avg(x*p*p,grid)),'support_identity':support_balance}
    checks={'ratios_finite_positive':bool(np.all(np.isfinite(ratio.y)) and np.min(eta)>0 and np.min(zeta)>0),'ratios_periodic':float(np.max(np.abs(ratio.y[:,0]-ratio.y[:,-1])))<1e-7*scale,'adjoint_exponent_agreement':abs(lc-lc_lin)<2e-8,'orthogonality':float(np.max(np.abs(O)))<2e-7*scale,'moment_identities':max(moments.values())<2e-6,'sampled_support_gaps':min(np.min(gap),np.min(Kgap),np.min(Egap))>-2e-8*scale,'sampled_curvature':np.min(W)>-2e-8*scale,'six_extrema_order':phase_order,'winding_minus_one':valid_winding(winding),'projective_bounds':np.min(eta)>=sminus-1e-8*scale and np.max(eta)<=splus+1e-8*scale and np.min(zeta)>=zminus-1e-8*scale and np.max(zeta)<=zplus+1e-8*scale,'first_return_detected':return_ok,'return_closure':closure<2e-7*scale,'positive_return_matrix':bool(M is not None and np.min(M)>0),'return_multiplier_agreement':abs(Mt-mt)<2e-7}
    result={'id':row['id'],'method':method,'mu':mu,'k':k,'rho':row['rho'],'amplitude':row['amplitude'],'period_fast':Tf,'lambda_c_projective':lc,'lambda_c_linear_adjoint':lc_lin,'exponent_difference':abs(lc-lc_lin),'strong_multiplier':ms,'weak_multiplier':mt,'matrix_weak_multiplier':Mt,'return_matrix':M,'return_closure_error':closure,'log_determinant_exact_formula':logdet,'determinant_float':float(np.linalg.det(M)) if M is not None else None,'determinant_caution':'Floating determinant can lose relative accuracy when exponentially small; no relative enclosure is claimed.','phase_fractions':phases,'stationary_counts':counts,'winding':winding,'moment_relative_residuals':moments,'orthogonality_max_scaled':float(np.max(np.abs(O)))/scale,'minimum_sampled_support_gap':float(np.min(gap)),'minimum_sampled_equilibrium_gap':float(np.min(Kgap)),'minimum_sampled_reverse_gap':float(np.min(Egap)),'minimum_sampled_curvature':float(np.min(W)),'projective_bounds':[sminus,splus,zminus,zplus],'observed_extrema':[float(np.min(eta)),float(np.max(eta)),float(np.min(zeta)),float(np.max(zeta))],'checks':checks,'passed':bool(all(checks.values())),'seconds':time.time()-start,'interval_certified':False}
    if M is not None:
        zz=initial[2];v=np.array([2*(1-np.exp(-initial[0]/2)),2*(zz-np.sqrt(zz))]);result['square_root_test_ratios']=M@v/v
        result['trace_obstruction']=float(np.trace(M)-1-np.exp(logdet))
    if save_profiles:
        path=ROOT/'evidence/profiles'/f"{row['id']}.npz";path.parent.mkdir(exist_ok=True)
        np.savez_compressed(path,phase=grid,eta=eta,zeta=zeta,U=U,p=p,y=y,curvature=W,support_gap=gap,mu=mu,k=k,period_fast=Tf)
        result['profile']=str(path.relative_to(ROOT))
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');ap.add_argument('--limit',type=int);ap.add_argument('--method',default='DOP853',choices=['DOP853','Radau']);ap.add_argument('--out',type=Path,default=ROOT/'evidence/numerical_audit.json');args=ap.parse_args()
    inputs=json.loads((ROOT/'evidence/cycle_trials.json').read_text())['attempts'];inputs=[r for r in inputs if r['accepted']]
    if args.limit:inputs=inputs[:args.limit]
    out=[]
    if args.resume and args.out.exists():out=json.loads(args.out.read_text())['records']
    for row in inputs:
        if any(r['id']==row['id'] for r in out):continue
        try:r=audit(row,args.method)
        except Exception as exc:r={'id':row['id'],'passed':False,'error':repr(exc)}
        out.append(r);print(row['id'],'PASS' if r['passed'] else 'FAIL',r.get('exponent_difference'),r.get('phase_fractions'),r.get('error',''),flush=True)
        if not r['passed']:print({k:v for k,v in r.get('checks',{}).items() if not v},flush=True)
        data={'kind':'floating_point_identity_and_geometry_audit','records':out,'number_attempted':len(out),'number_passed':sum(r['passed'] for r in out),'passed':all(r['passed'] for r in out),'claim_boundary':'Sampled diagnostics of accepted approximations, not proof of all-cycle properties or interval enclosures.'}
        strict_write(args.out,data)
    return 0 if all(r['passed'] for r in out) else 1
if __name__=='__main__':sys.exit(main())
