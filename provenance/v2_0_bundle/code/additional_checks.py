#!/usr/bin/env python3
"""Cross-method, reciprocal-rate, convex-envelope, and negative-control checks.
All numerical checks remain floating-point diagnostics, never certificates.
"""
from __future__ import annotations
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from cycles import ROOT,strict_write
from numerical_audit import audit,load,valid_phase_order,valid_winding,radial_certificate_pass


def run():
    inputs=json.loads((ROOT/'evidence/cycle_trials.json').read_text())['attempts']
    primary=json.loads((ROOT/'evidence/numerical_audit.json').read_text())['records']
    byid={r['id']:r for r in inputs};pa={r['id']:r for r in primary}
    # Independently re-integrate with an implicit solver, preserving primary profiles.
    cross=[]
    for ident in ['mu_0.3_h_1','mu_1_h_3','mu_10_h_2','mu_0.3_h_0.01']:
        try:
            other=audit(byid[ident],method='Radau',save_profiles=False);one=pa[ident]
            ediff=abs(other['lambda_c_projective']-one['lambda_c_projective']);mdiff=float(np.max(np.abs(np.array(other['return_matrix'])-np.array(one['return_matrix']))))
            passed=other['passed'] and ediff<1e-8 and mdiff<1e-7
            cross.append({'id':ident,'passed':bool(passed),'exponent_difference':ediff,'return_matrix_max_difference':mdiff,'Radau_record':other})
            print('cross-method',ident,passed,ediff,mdiff,flush=True)
        except Exception as exc:cross.append({'id':ident,'passed':False,'error':repr(exc)})
    strict_write(ROOT/'evidence/cross_method.json',{'kind':'floating_point_solver_sensitivity','records':cross,'number_checked':len(cross),'number_passed':sum(r['passed'] for r in cross),'passed':all(r['passed'] for r in cross)})
    reciprocal=[]
    for m1,m2 in [(.1,10.),(.3,10/3)]:
        for h in [.05,.2,.5,1.,2.,3.]:
            i1=f'mu_{m1:.8g}_h_{h:.8g}';i2=f'mu_{m2:.8g}_h_{h:.8g}'
            s1,mu1,l1,T1=load(ROOT/byid[i1]['mesh']);s2,mu2,l2,T2=load(ROOT/byid[i2]['mesh']);g=np.linspace(0,1,5001)
            u1,v1,z1=s1(g);u2,v2,z2=s2(g)
            dy=float(np.max(np.abs((1-np.sqrt(mu1)*v1/l1)-(1-np.sqrt(mu2)*v2/l2))))
            du=float(np.max(np.abs(u1-u2)));drho=abs(l1/(1+mu1)-l2/(1+mu2));dt=abs(T1-T2)
            ok=max(du,dy,drho,dt)<2e-7
            reciprocal.append({'mu_pair':[m1,m2],'amplitude':h,'max_difference_U':du,'max_difference_Y':dy,'difference_rho':drho,'difference_normalised_period':dt,'passed':ok})
    strict_write(ROOT/'evidence/reciprocal_rates.json',{'kind':'floating_point_coordinate_symmetry','records':reciprocal,'number_checked':len(reciprocal),'number_passed':sum(r['passed'] for r in reciprocal),'passed':all(r['passed'] for r in reciprocal)})
    # Verify the local envelope Hessian by separately locating the maximising phase.
    ident='mu_1_h_1';state,mu,lam,Ttau=load(ROOT/byid[ident]['mesh']);k=mu*lam;Tf=Ttau/np.sqrt(mu)
    prof=np.load(ROOT/pa[ident]['profile']);eg=prof['eta'].copy();zg=prof['zeta'].copy()
    correction=max(abs(eg[-1]-eg[0]),abs(zg[-1]-zg[0]));eg[-1]=eg[0];zg[-1]=zg[0]
    es=CubicSpline(prof['phase'],eg,bc_type='periodic');zs=CubicSpline(prof['phase'],zg,bc_type='periodic')
    def point(t):
        U,V,z=state(t);e=es(t%1);q=zs(t%1);x=k*np.exp(U);p=np.sqrt(mu)*V;y=k*z
        gp=np.array([x-(1+e)*q,-(q+(mu-1)*e-e*e)])
        velocity=np.array([p,k-mu*p-y])
        return np.array([U,p]),y,np.array([q,-e]),gp,velocity
    se=brentq(lambda e:e*(e+1)*(e+1-mu)-k,max(0.,mu-1),mu+k+2);ze=se*(se+1-mu)
    ts=np.arange(4096)/4096;uv=np.array([point(t)[0] for t in ts]);yy=np.array([point(t)[1] for t in ts]);gg=np.array([point(t)[2] for t in ts])
    def optimum(query,base):
        def f(t):r,y,g,gp,v=point(t);return np.dot(gp,query-r)
        root=brentq(f,base-.01,base+.01,xtol=1e-14,rtol=1e-14)
        r,y,g,gp,v=point(root);height=y+np.dot(g,query-r)
        gridheight=yy+np.sum(gg*(query-uv),axis=1)
        win=height>=max(float(np.max(gridheight)),k+ze*query[0]-se*query[1])-1e-10
        return root,height,g,bool(win)
    envelopes=[]
    for ph in [.08,.35,.55,.8]:
        r,y,g,gp,velocity=point(ph);W=float(np.dot(gp,velocity));H=np.outer(gp,gp)/W
        step=2e-5;Hfd=np.zeros((2,2));wins=[];resid=[]
        for j in range(2):
            v=np.eye(2)[j]*step
            op=optimum(r+v,ph);om=optimum(r-v,ph);Hfd[:,j]=(op[2]-om[2])/(2*step);wins.extend([op[3],om[3]])
            for query,out in [(r+v,op),(r-v,om)]:
                tt,hh,gr,_=out;r0,y0,_,_,_=point(tt);x0=k*np.exp(r0[0]);du=query[0]-r0[0]
                residual=gr[0]*query[1]+gr[1]*(k-mu*query[1]-hh)-(k*np.exp(query[0])-hh)
                resid.append(abs(residual+x0*(np.expm1(du)-du)))
        err=float(np.max(np.abs(Hfd-H)))/max(1.,float(np.max(np.abs(H))))
        ok=W>0 and all(wins) and err<1e-4 and max(resid)<1e-7
        envelopes.append({'phase':ph,'curvature_W':W,'predicted_Hessian':H,'finite_difference_Hessian':Hfd,'relative_scaled_error':err,'defect_identity_max_residual':max(resid),'maximiser_beats_grid_and_equilibrium':all(wins),'passed':bool(ok)})
        print('envelope',ph,ok,err,max(resid),flush=True)
    strict_write(ROOT/'evidence/envelope_checks.json',{'kind':'floating_point_convex_envelope_diagnostic','base_cycle':ident,'periodic_spline_endpoint_adjustment':float(correction),'finite_difference_step':2e-5,'comparison_phase_grid':4096,'records':envelopes,'number_checked':len(envelopes),'number_passed':sum(r['passed'] for r in envelopes),'passed':all(r['passed'] for r in envelopes),'scope':'The finite phase grid does not certify a global maximum or a continuum inequality.'})
    # Failure semantics: the test runner succeeds only if deliberate errors fail.
    wrong=ROOT/'evidence/deliberately_wrong_exact_receipt.json'
    result=subprocess.run([sys.executable,str(ROOT/'code/exact_checks.py'),'--inject-error','--out',str(wrong)],text=True,capture_output=True)
    (ROOT/'logs/deliberate_wrong_exact.log').write_text(result.stdout+result.stderr)
    record=json.loads(wrong.read_text())
    controls=[{'id':'wrong_identity_rejected_by_predicate_and_exit','passed':result.returncode!=0 and not record['passed'],'observed_exit_status':result.returncode,'observed_scientific_pass':record['passed']}]
    rr=pa['mu_0.3_h_1'];seq=list(rr['phase_fractions'].values())
    ordered=valid_phase_order
    controls.append({'id':'reversed_extrema_order_rejected','passed':ordered(seq) and not ordered(seq[::-1]),'correct_order_pass':ordered(seq),'reversed_order_pass':ordered(seq[::-1])})
    controls.append({'id':'wrong_winding_rejected','passed':valid_winding(rr['winding']) and not valid_winding(1.),'observed_winding':rr['winding']})
    small=pa['mu_0.3_h_0.01'];v=np.array(small['square_root_test_ratios'])
    controls.append({'id':'false_radial_certificate_detected','passed':not radial_certificate_pass(v) and float(v[0])>1 and small['weak_multiplier']<1,'radial_ratio_1':float(v[0]),'actual_weak_multiplier_approximation':small['weak_multiplier'],'certificate_pass':radial_certificate_pass(v)})
    strict_write(ROOT/'evidence/negative_controls.json',{'kind':'deliberate_error_regression','controls':controls,'number_checked':len(controls),'number_passed':sum(r['passed'] for r in controls),'passed':all(r['passed'] for r in controls),'interpretation':'A control passes only when the deliberately wrong scientific assertion is rejected.'})
    return all(r['passed'] for r in cross+reciprocal+envelopes+controls)
if __name__=='__main__':sys.exit(0 if run() else 1)
