#!/usr/bin/env python3
"""Exact finite algebra underlying the manuscript. Not a formal proof of analysis.
Run from any directory. --out chooses a receipt; --inject-error is a negative control.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]

def run(inject_error: bool = False) -> dict:
    mu,k,x,p,y,e,z,U = s.symbols('mu k x p y eta zeta U', positive=False)
    A=1+mu
    fp=k-mu*p-y
    ed=z+(mu-1)*e-e**2
    zd=x-(1+e)*z
    states=(x,p,y,e,z,U)
    field=(p*x,fp,x-y,ed,zd,p)
    def D(f): return s.expand(sum(s.diff(f,a)*b for a,b in zip(states,field)))
    rows=[]
    def ck(name, expression):
        r=s.simplify(s.expand(expression))
        passed = r == 0 or (isinstance(r,s.MatrixBase) and r == s.zeros(*r.shape))
        rows.append({'id':name,'passed':bool(passed),'residual':str(r)})
    J=s.Matrix([[0,1,0],[0,-mu,-1],[x,0,-1]])
    S=s.diag(1,-1,-1)
    B=s.Matrix([[0,0,x],[1,mu,0],[0,1,1]])
    ck('adjoint_sign_conjugacy',B+S*J.T*S)
    h1,h2,h3=s.symbols('h1 h2 h3',positive=True)
    H=s.Matrix([h1,h2,h3]); Hd=B*H
    ck('projective_eta',((Hd[1]*h3-h2*Hd[2])/h3**2).subs({h1:z*h3,h2:e*h3})-ed)
    ck('projective_zeta',((Hd[0]*h3-h1*Hd[2])/h3**2).subs({h1:z*h3,h2:e*h3})-zd)
    n=s.Matrix([z,-e,-1]); nd=s.Matrix([zd,-ed,0])
    ck('normal_adjoint_equation',nd+J.T*n+(1+e)*n)
    O=z*p-e*fp-(x-y)
    ck('orthogonality_transport',D(O)+(1+e)*O)
    pdd=D(fp)
    ck('tangent_second_order',pdd-(e-mu)*fp+z*p-O)
    ck('operator_D2',mu-e+1+e-A)
    ck('operator_D1',z-ed+(1+e)*(mu-e)-mu)
    ck('operator_D0',zd+(1+e)*z-x)
    f=e*(e+1)*(e+1-mu)
    ck('gain_primitive',D(z+e+e**2/2)-x+f)
    ck('gain_polynomial',f-((1-mu)*e+(2-mu)*e**2+e**3))
    K=z*U-e*p-y+k
    ck('support_equilibrium_primitive',D(K)+(1+e)*K-(x*U-x+k))
    W=zd*p-ed*fp
    ck('support_curvature',D(W)+(1+e)*W-x*p**2+ed*O)
    primitive=(1+e)*p**2/2+x-k*U
    integrand=((1+e)*(e-mu)+s.Rational(3,2)*ed)*p**2
    ck('weighted_tangent_primitive',D(primitive)-integrand+(1+e)*p*O)
    # The two-trajectory remainder uses x2=x exp(dU); no Taylor truncation.
    du,dp,dy=s.symbols('du dp dy')
    dexp=s.exp(du)
    delta=z*du-e*dp-dy
    deltad=zd*du+z*dp-ed*dp-e*(-mu*dp-dy)-(x*(dexp-1)-dy)
    ck('strict_support_identity',deltad+(1+e)*delta+x*(dexp-1-du))
    # Tangent identities at extrema.
    ck('crest_identity_x', (O-(x-k-(1+e)*(y-k))*(-1)).subs(p,0))
    ck('filter_second_derivative',D(x-y).subs(y,x)-x*p)
    ck('p_extremum_sign',(pdd+z*p-O).subs(y,k-mu*p))
    # Compound matrix and crest wedge.
    C2=s.Matrix([[-mu,-1,0],[0,-1,1],[-x,0,-A]])
    C=s.Matrix([[-mu,1,0],[0,-1,1],[x,0,-A]])
    ck('compound_sign_conjugacy',S*C2*S-C)
    a,b=s.symbols('variation_U variation_y')
    F=s.Matrix([p,fp,x-y]); V=s.Matrix([a,0,b])
    wedge=s.Matrix([F[0]*V[1]-F[1]*V[0],F[0]*V[2]-F[2]*V[0],F[1]*V[2]-F[2]*V[1]])
    ck('crest_wedge',(S*wedge).subs(p,0)-s.Matrix([(y-k)*a,(x-y)*a,(y-k)*b]))
    ms,mt=s.symbols('m_s m_t')
    ck('return_trace_sign',ms+mt-1-ms*mt-(mt-1)*(1-ms))
    ck('neutral_variance_expansion',f.subs(e,mu+s.Symbol('q'))-(mu*A+(1+3*mu+mu**2)*s.Symbol('q')+2*A*s.Symbol('q')**2+s.Symbol('q')**3))
    ck('neutral_variance_lower_integrand',2*A*s.Symbol('q')**2+s.Symbol('q')**3-(2+mu)*s.Symbol('q')**2-s.Symbol('q')**2*(s.Symbol('q')+mu))
    # At equilibrium, identify the stable direction and the positive normal.
    g=s.symbols('gamma',positive=True); kg=g*(g-1)*(g-mu)
    Je=J.subs(x,k)
    vs=s.Matrix([1,-g,k/(1-g)])
    ck('stable_equilibrium_eigenvector',((Je+g*s.eye(3))*vs).subs(k,kg))
    ne=s.Matrix([(g-1)*(g-mu),-(g-1),-1])
    ck('equilibrium_normal',((Je.T+g*s.eye(3))*ne).subs(k,kg))
    ck('equilibrium_characteristic',Je.charpoly(s.Symbol('L')).as_expr()-(s.Symbol('L')**3+A*s.Symbol('L')**2+mu*s.Symbol('L')+k))
    ck('critical_characteristic',(s.Symbol('L')**3+A*s.Symbol('L')**2+mu*s.Symbol('L')+k).subs(k,mu*A)-(s.Symbol('L')+A)*(s.Symbol('L')**2+mu))
    # Explicit periodic ratio identities for the all-cycle state box.
    X,Y,Z,lam=s.symbols('X Y Z lam',positive=True)
    FX=lam*X*(1-Y); FY=mu*(Z-Y); FZ=X-Z
    def PD(q):return s.factor(s.diff(q,X)*FX+s.diff(q,Y)*FY+s.diff(q,Z)*FZ)
    ck('physical_Z_over_X',PD(Z/X)-(1-(1+lam*(1-Y))*Z/X))
    ck('physical_Y_over_X',PD(Y/X)-(mu*Z/X-(mu+lam*(1-Y))*Y/X))
    ck('physical_transform_U',FX/X-(lam*(1-Y)))
    ck('physical_transform_p',(-lam*FY).subs(lam,k/mu)-(k-mu*(k/mu)*(1-Y)-k*Z))
    ck('physical_transform_y',k*FZ-(k*X-k*Z))
    ck('reciprocal_time_D2',A/mu-(1+1/mu))
    ck('reciprocal_time_D1',mu/mu**2-1/mu)
    # Convex envelope: derivative, Hessian, defect and positivity algebra.
    u0,p0,y0,x0,z0,e0,d0,d1=s.symbols('u0 p0 y0 x0 z0 e0 d0 d1')
    ell=y0+z0*d0-e0*d1
    xp=x0*s.exp(d0)
    zbprime=x0-(1+e0)*z0
    ebprime=z0+(mu-1)*e0-e0**2
    baseO=z0*p0-e0*(k-mu*p0-y0)-(x0-y0)
    residual=z0*(p0+d1)-e0*(k-mu*(p0+d1)-ell)-(xp-ell)
    stationarity=zbprime*d0-ebprime*d1
    ck('envelope_graph_defect',residual-baseO+stationarity+x0*(s.exp(d0)-1-d0))
    gu,gp,vU,vP=s.symbols('grad_dot_U grad_dot_p vel_U vel_p')
    curvature=gu*vU+gp*vP
    HH=s.Matrix([gu,gp])*s.Matrix([[gu,gp]])/curvature
    ck('envelope_Hessian_rank_one',s.det(HH))
    ck('envelope_Hessian_tangent',HH*s.Matrix([vU,vP])-s.Matrix([gu,gp]))
    # Cubic Hopf coefficient in the frequency-one convention.
    sig=s.symbols('sigma',positive=True); I=s.I
    J0=s.Matrix([[0,1,0],[0,0,1],[-sig,-1,-sig]])
    q=s.Matrix([1,I,-1]); left=s.Matrix([[I*sig,sig+I,1]])/(-2+2*I*sig)
    ck('Hopf_right_eigenvector',J0*q-I*q)
    ck('Hopf_left_eigenvector',left*J0-I*left)
    ck('Hopf_normalisation',(left*q)[0]-1)
    def QB(v,w):return s.Matrix([0,0,-sig*v[0]*w[0]])
    def QC(v,w,t):return s.Matrix([0,0,-sig*v[0]*w[0]*t[0]])
    qbar=s.conjugate(q)
    G=(left*(QC(q,q,qbar)-2*QB(q,J0.inv()*QB(q,qbar))+QB(qbar,(2*I*s.eye(3)-J0).inv()*QB(q,q))))[0]
    ck('Hopf_G21',G-sig*(sig+3*I)/(3*(sig+2*I)*(-1+I*sig)))
    ck('Hopf_l1',s.re(s.simplify(G))/2+sig/((sig**2+1)*(sig**2+4)))
    ck('Hopf_transversality',-sig/(3*I**2+2*sig*I+1)-sig*(1+I*sig)/(2*(1+sig**2)))
    # The earlier false radial test is reproduced as an algebraic caution.
    Dm=mu**2+6*mu+1; Cm=mu**2+3*mu+1
    c20=mu*(mu**2+4*mu+1)/(2*Dm); c11=-mu*A/Dm; c02=mu/Dm
    Hg=mu*U-mu*p+c20*U**2+c11*U*p+c02*p**2
    ff=s.expand(s.diff(Hg,U)*p+s.diff(Hg,p)*(-mu*p-Hg)-(mu*A*(U+U**2/2)-Hg))
    for ij in [(1,0),(0,1),(2,0),(1,1),(0,2)]:
        ck('centre_graph_'+''.join(map(str,ij)),ff.coeff(U,ij[0]).coeff(p,ij[1]))
    h,qf=s.symbols('h q_f')
    zcurve=1+h/A+c20/(mu*A)*h**2
    v1=2*(1-s.exp(-h/2)); v2=2*(zcurve-s.sqrt(zcurve))
    rv=s.series(v2/v1,h,0,2).removeO(); re=1/A+2*c20/(mu*A)*h
    mismatch=mu*(1-mu)*(mu+3)/(4*A**2*Dm)
    ck('radial_direction_mismatch',s.expand(rv-re).coeff(h,1)-mismatch)
    coeff=(1-qf)*A*mu/Cm*mismatch
    ck('radial_failure_coefficient',coeff-(1-qf)*mu**2*(1-mu)*(mu+3)/(4*A*Cm*Dm))
    # A deliberately false statement must not pass.
    if inject_error:ck('injected_wrong_support_sign',2*x*(s.exp(du)-1-du))
    return {'schema_version':'1.0','kind':'exact_symbolic_equalities','software':{'sympy':s.__version__},'checks':rows,'number_checked':len(rows),'number_passed':sum(r['passed'] for r in rows),'passed':all(r['passed'] for r in rows),'scope':'Finite algebra only. Analysis and topology require the written proofs.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=ROOT/'evidence/exact_checks.json');ap.add_argument('--inject-error',action='store_true');args=ap.parse_args()
    data=run(args.inject_error);args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(f"{data['number_passed']}/{data['number_checked']} exact equalities passed")
    for row in data['checks']:
        if not row['passed']:print(row)
    return 0 if data['passed'] else 1
if __name__=='__main__':sys.exit(main())
