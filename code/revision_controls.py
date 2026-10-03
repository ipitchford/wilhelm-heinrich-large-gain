"""Bounded regression and semantic controls for the v3.1 revisions; not a proof checker."""
import importlib.util
from pathlib import Path
import sympy as s
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('interval_audit',root/'code/rational_interval_audit.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
k=s.symbols('k',positive=True)
actual=s.diff(3*k**s.Rational(2,3)-4*k**s.Rational(1,3)-1,k)
correct=s.Rational(2,3)*k**s.Rational(-2,3)*(3*k**s.Rational(1,3)-2)
wrong=2*k**s.Rational(-2,3)*(k**s.Rational(1,3)-1)
if s.simplify(actual-correct)!=0 or s.simplify(actual-wrong)==0:raise ValueError('derivative regression')
u,v=s.symbols('u v');h=u*u+v*v;F=s.Matrix([v,u])
grad=s.Matrix([s.diff(h,u),s.diff(h,v)])
material=s.diff((grad.dot(F)),u)*F[0]+s.diff((grad.dot(F)),v)*F[1]
quadratic=(F.T*s.hessian(h,(u,v))*F)[0]
if s.simplify(material-quadratic-grad.dot(F.jacobian([u,v])*F))!=0:raise ValueError('chain rule')
if s.simplify(material-quadratic)==0:raise ValueError('missing curvature term undetected')
for bad in [lambda:m.sqrt(-1),lambda:m.log(0),lambda:m.I(1)/m.I(-1,1),lambda:m.chk('wrong contraction bound',m.I('.82717'),'<=' ,m.I('.82'))]:
    try:bad()
    except (ValueError,AssertionError):pass
    else:raise ValueError('semantic negative control accepted corruption')
print('PASS: derivative, material derivative, and four interval-domain/false-bound rejection controls; optimization-safe.')
