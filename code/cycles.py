#!/usr/bin/env python3
"""Fresh amplitude-fixed periodic BVPs; diagnostic, not a global exclusion.
State (U,V,z), time tau=sqrt(mu)*t. All rejected solves are retained.
Exponential clipping is not used. A guard aborts an invalid attempted solve.
"""
from __future__ import annotations
import argparse,json,sys,time,warnings
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import PPoly
import scipy
ROOT=Path(__file__).resolve().parents[1]

def strict_write(path:Path,data:dict):
    def clean(o):
        if isinstance(o,dict):return {str(k):clean(v) for k,v in o.items()}
        if isinstance(o,(list,tuple)):return [clean(v) for v in o]
        if isinstance(o,np.ndarray):return clean(o.tolist())
        if isinstance(o,np.integer):return int(o)
        if isinstance(o,(float,np.floating)):return float(o) if np.isfinite(o) else None
        if isinstance(o,np.bool_):return bool(o)
        return o
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(clean(data),indent=2,allow_nan=False)+'\n')

def equations(mu,state,lam):
    U,V,z=state
    if not np.all(np.isfinite(state)) or np.max(np.abs(U))>600 or not np.isfinite(lam):
        raise FloatingPointError('state/parameter guard; no clipping was applied')
    return np.array([V,lam*(1-z)-np.sqrt(mu)*V,(np.exp(U)-z)/np.sqrt(mu)])

def load_cycle(path:Path):
    q=np.load(path)
    pp=PPoly(q['coefficients'],q['breaks'],extrapolate=False)
    return SimpleNamespace(sol=lambda t: pp(t).T,p=np.array([float(q['period_tau']),float(q['lam'])])),float(q['h'])

def solve_cycle(mu,h,previous=None,tol=2e-9,max_nodes=22000):
    start=time.time();mesh=np.linspace(0,1,501)
    if previous is None:
        theta=2*np.pi*mesh
        yy=np.array([h*np.cos(theta),-h*np.sin(theta),1+h*(np.cos(theta)+np.sqrt(mu)*np.sin(theta))/(1+mu)])
        par=np.array([2*np.pi*(1+.02*h*h),(1+mu)*(1+.03*h*h)])
    else:
        old,oldh=previous;yy=old.sol(mesh);scale=h/oldh
        yy[0]*=scale;yy[1]*=scale;yy[2]=1+(yy[2]-1)*scale
        par=old.p.copy();par[1]=(1+mu)+(par[1]-(1+mu))*scale**2
    def fun(t,Y,pp):
        T,lam=pp
        if not np.isfinite(T) or abs(T)>2000:raise FloatingPointError('period guard')
        return T*equations(mu,Y,lam)
    def jac(t,Y,pp):
        T,lam=pp;U,V,z=Y;n=len(t);rt=np.sqrt(mu)
        J=np.zeros((3,3,n));J[0,1]=T;J[1,1]=-T*rt;J[1,2]=-T*lam;J[2,0]=T*np.exp(U)/rt;J[2,2]=-T/rt
        jp=np.zeros((3,2,n));jp[:,0]=equations(mu,Y,lam);jp[1,1]=T*(1-z)
        return J,jp
    def bc(a,b,pp):return np.r_[a-b,a[0]-h,a[1]]
    row={'mu':mu,'amplitude':h,'tolerance':tol,'max_nodes':max_nodes,'guard_triggered':False,'accepted':False}
    try:
        with warnings.catch_warnings(record=True) as ww:
            warnings.simplefilter('always')
            sol=solve_bvp(fun,bc,mesh,yy,p=par,tol=tol,max_nodes=max_nodes,fun_jac=jac,bc_tol=tol)
        row.update(solver_success=bool(sol.success),solver_status=int(sol.status),message=sol.message,nodes=len(sol.x),max_rms_residual=float(np.max(sol.rms_residuals)),boundary_residual=float(np.max(np.abs(bc(sol.y[:,0],sol.y[:,-1],sol.p)))),warnings=list(dict.fromkeys(str(w.message) for w in ww)))
        T,lam=map(float,sol.p);grid=np.linspace(0,1,12001);values=sol.sol(grid);U,V,z=values;Y=1-np.sqrt(mu)*V/lam
        checks={'finite':bool(np.all(np.isfinite(values)) and np.all(np.isfinite(sol.p))),'positive_period':T>0,'period_above_bound':T>2*np.pi,'oscillatory_parameter':lam>1+mu,'physical_positive':bool(np.min(z)>0 and np.min(Y)>0),'nonconstant':np.ptp(U)>h,'crest_phase':z[0]>1 and abs(V[0])<50*tol,'residual_small':np.max(sol.rms_residuals)<2*tol and row['boundary_residual']<2*tol}
        row.update(period_tau=T,period_fast=T/np.sqrt(mu),lam=lam,k=mu*lam,rho=lam/(1+mu),u_min=float(np.min(U)),u_max=float(np.max(U)),min_z=float(np.min(z)),min_Y=float(np.min(Y)),acceptance_checks=checks,accepted=bool(sol.success and all(checks.values())))
        row['seconds']=time.time()-start
        return sol,row
    except (FloatingPointError,ValueError,OverflowError,np.linalg.LinAlgError) as exc:
        row.update(solver_success=False,guard_triggered=isinstance(exc,FloatingPointError),message=str(exc),seconds=time.time()-start)
        return None,row

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'evidence/cycle_trials.json');ap.add_argument('--only-mu',type=float);ap.add_argument('--resume',action='store_true');args=ap.parse_args()
    mus=[.1,.3,1.,10/3,10.] if args.only_mu is None else [args.only_mu]
    hs=[.05,.2,.5,1.,2.,3.];rows=[]
    if args.resume and args.out.exists():rows=json.loads(args.out.read_text())['attempts']
    folder=ROOT/'evidence/cycles';folder.mkdir(exist_ok=True)
    data={'schema_version':'1.0','kind':'floating_point_amplitude_fixed_bvp','not_interval_certified':True,'software':{'numpy':np.__version__,'scipy':scipy.__version__},'attempts':rows}
    for mu in mus:
        previous=None
        for h in hs:
            ident=f'mu_{mu:.8g}_h_{h:.8g}'
            old=next((r for r in rows if r['id']==ident),None)
            if old:
                if old['accepted']:previous=load_cycle(ROOT/old['mesh'])
                continue
            sol,row=solve_cycle(mu,h,previous=previous);row['id']=ident
            if row['accepted']:
                meshfile=folder/(ident+'.npz');np.savez_compressed(meshfile,breaks=sol.sol.x,coefficients=sol.sol.c,mu=mu,h=h,period_tau=sol.p[0],lam=sol.p[1],mesh=sol.x,states=sol.y)
                row['mesh']=str(meshfile.relative_to(ROOT));previous=(sol,h)
            rows.append(row);print(ident,'ACCEPTED' if row['accepted'] else 'REJECTED',row.get('rho'),row.get('nodes'),flush=True)
            data.update(number_attempted=len(rows),number_accepted=sum(r['accepted'] for r in rows));strict_write(args.out,data)
    ident='mu_0.3_h_0.01'
    if args.only_mu is None and not any(r['id']==ident for r in rows):
        mu=.3;h=.01;sol,row=solve_cycle(mu,h,tol=2e-11,max_nodes=30000);row['id']=ident
        if row['accepted']:
            meshfile=folder/(ident+'.npz');np.savez_compressed(meshfile,breaks=sol.sol.x,coefficients=sol.sol.c,mu=mu,h=h,period_tau=sol.p[0],lam=sol.p[1],mesh=sol.x,states=sol.y);row['mesh']=str(meshfile.relative_to(ROOT))
        rows.append(row)
    data.update(number_attempted=len(rows),number_accepted=sum(r['accepted'] for r in rows));strict_write(args.out,data)
    print(data['number_accepted'],'/',data['number_attempted'],'accepted',flush=True)
    return 0
if __name__=='__main__':sys.exit(main())
