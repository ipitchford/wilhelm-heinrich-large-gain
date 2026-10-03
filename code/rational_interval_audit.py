#!/usr/bin/env python3
"""Exact-rational, outward-rounded enclosure of the release's scalar constants.

No third-party arithmetic library is used. All decision arithmetic uses Fraction,
integer isqrt and rational series with explicit remainder bounds. Elementary
results are rounded OUTWARDS onto a 10^-50 rational grid. This audits the finite
endpoint comparisons, NOT the continuum monotonicity arguments in the paper.
The mathematical expressions are replayed from the supplied author checker;
the interval arithmetic and root enclosure are new referee code.
Usage: python rational_interval_audit.py /path/to/lg_constants_author.py
"""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import sys, json, hashlib
from decimal import Decimal, localcontext
D=10**50

def rat(x):
    if isinstance(x,Q): return x
    if isinstance(x,float): return Q(str(x))  # the original script's finite decimal literals
    return Q(x)
def down(x): return Q((x.numerator*D)//x.denominator,D)
def up(x): return -down(-x)
def fmt(x):
    with localcontext() as c:
        c.prec=22
        return str(Decimal(x.numerator)/Decimal(x.denominator))
class I:
    def __init__(self,a,b=None):
        if isinstance(a,I) and b is None: self.lo,self.hi=a.lo,a.hi; return
        a=rat(a); b=a if b is None else rat(b)
        if a>b: raise ValueError('reversed interval')
        self.lo,self.hi=down(a),up(b)
    def __add__(self,y):
        y=I(y); return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,y): return self+-I(y)
    def __rsub__(self,y): return I(y)+-self
    def __mul__(self,y):
        y=I(y); vs=[a*b for a in (self.lo,self.hi) for b in (y.lo,y.hi)]
        return I(min(vs),max(vs))
    __rmul__=__mul__
    def __truediv__(self,y):
        y=I(y)
        if y.lo<=0<=y.hi: raise ValueError('division through zero')
        return self*I(1/y.hi,1/y.lo)
    def __rtruediv__(self,y): return I(y)/self
    def __pow__(self,n):
        n=I(n)
        if n.lo==n.hi==Q(3,2): return sqrt(self)**3
        if n.lo!=n.hi or n.lo.denominator!=1: raise ValueError('unsupported power')
        n=int(n.lo)
        if n<0: return 1/(self**(-n))
        r=I(1); a=self
        while n:
            if n%2:r=r*a
            a=a*a; n//=2
        return r
    def __le__(self,y): return self.hi<=I(y).lo
    def __ge__(self,y): return self.lo>=I(y).hi
    def __lt__(self,y): return self.hi<I(y).lo
    def __gt__(self,y): return self.lo>I(y).hi
    def __str__(self): return '['+fmt(self.lo)+', '+fmt(self.hi)+']'
    __repr__=__str__

def interval_max(*values):
    values=[I(v) for v in values]
    return I(max(v.lo for v in values), max(v.hi for v in values))

def sqrt(x):
    x=I(x)
    if x.lo<0: raise ValueError('negative square root')
    a=isqrt((x.lo.numerator*D*D)//x.lo.denominator)
    b=isqrt((x.hi.numerator*D*D)//x.hi.denominator)
    lo=Q(a,D); hi=Q(b,D)
    if hi*hi!=x.hi: hi+=Q(1,D)
    return I(lo,hi)

def exp_point(x):
    if x<0:return 1/exp_point(-x)
    m=0;r=x
    while r>Q(1,8):r/=2;m+=1
    term=Q(1);total=term
    N=60
    for j in range(1,N+1):term=term*r/j;total+=term
    nxt=term*r/(N+1)
    # For all subsequent terms the ratio is <= r/(N+2).
    ans=I(total,total+nxt/(1-r/(N+2)))
    for _ in range(m):ans=ans*ans
    return ans

def exp(x):
    x=I(x); return I(exp_point(x.lo).lo,exp_point(x.hi).hi)

def log_unit(r):
    if not 1<=r<=2:raise ValueError('log range')
    z=(r-1)/(r+1); total=Q(0); power=z;N=80
    for j in range(N):total+=2*power/(2*j+1);power*=z*z
    tail=2*power/((2*N+1)*(1-z*z))
    return I(total,total+tail)
LN2=log_unit(Q(2))
def log_point(x):
    if x<=0:raise ValueError('nonpositive logarithm')
    n=0;r=x
    while r>=2:r/=2;n+=1
    while r<1:r*=2;n-=1
    return n*LN2+log_unit(r)
def log(x):
    x=I(x);return I(log_point(x.lo).lo,log_point(x.hi).hi)

def atan_small(z):
    # Alternating-series enclosure for 0<z<1; N even means next term positive.
    N=80; s=sum(((-1)**j*z**(2*j+1)/Q(2*j+1) for j in range(N)),Q(0))
    nxt=z**(2*N+1)/Q(2*N+1)
    return I(s,s+nxt)
PI=16*atan_small(Q(1,5))-4*atan_small(Q(1,239))
E=exp(1)
# A certified positive root bracket for v=2(1-exp(-v)).
vl=I('1.59362426004004009'); vh=I('1.59362426004004010')
f=lambda v:v-2*(1-exp(-v))
if not (f(vl)<0 and f(vh)>0 and 1-2*exp(-vl)>0):
    raise ValueError('positive root bracket or derivative certificate failed')
V2=I(vl.lo,vh.hi)
def findroot(fn,guess):return V2
class MP:
    dps=0
    @staticmethod
    def nstr(x,n=8):return str(I(x))
mp=MP(); records=[]
def chk(label,lhs,op,rhs):
    lhs,rhs=I(lhs),I(rhs)
    margin=lhs-rhs if op=='>=' else rhs-lhs
    ok=margin.lo>=0
    row={'name':label,'operator':op,'passed':ok,
         'lhs_interval':[str(lhs.lo),str(lhs.hi)],'rhs_interval':[str(rhs.lo),str(rhs.hi)],
         'margin_lower_rational':str(margin.lo),'margin_lower_decimal':fmt(margin.lo)}
    records.append(row)
    print(('PASS ' if ok else 'FAIL ')+label+'; margin >= '+fmt(margin.lo))
    if not ok:raise AssertionError(label)

def main():
    src=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'v3/lib/lg_constants_author.py'
    text=src.read_text()
    # Preserve all author expressions and comparisons, replace only the arithmetic
    # backend, output helper, numerical root finder, and exact cube roots of 1000.
    start=text.index('# ---------------------------------------------------------------- fixed constants')
    body=text[start:]
    body=body.replace("K0**(mpf(2)/3)","I(100)").replace("K0**(mpf(1)/3)","I(10)")
    ns=dict(I=I,mp=mp,mpf=I,exp=exp,log=log,sqrt=sqrt,pi=PI,findroot=findroot,E=E,
            K0=I(1000),FAILS=[],RECORDS=[],chk=chk,sys=sys,max=interval_max)
    exec(compile(body,str(src),'exec'),ns)
    chk('REFEREE: theorem bound 3/4 + q <= 0.83',I(3)/4+ns['q'],'<=',I('0.83'))
    out={'method':'exact-rational interval arithmetic; explicit series remainder bounds',
         'grid_denominator':str(D),'author_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
         'scope':'Finite endpoint comparisons only; no automatic verification of the analytic continuum reductions.',
         'root_bracket_v2':[str(V2.lo),str(V2.hi)],'root_signs_certified':True,
         'q_interval':[str(ns['q'].lo),str(ns['q'].hi)],'q_interval_display':str(ns['q']),
         'checks':len(records),'passed':all(r['passed'] for r in records),'records':records}
    output = Path(sys.argv[2]) if len(sys.argv)>2 else Path(__file__).resolve().parents[1]/'evidence'/'rational_interval_results.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    print('EXACT-RATIONAL INTERVAL CHECKS PASSED:',len(records))
if __name__=='__main__':main()
