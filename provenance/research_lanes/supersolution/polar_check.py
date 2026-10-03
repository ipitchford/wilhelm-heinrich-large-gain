#!/usr/bin/env python3
"""Check the polar orbit equation used in Proposition T:
   d log r / d theta = sin(th) (m + mu sin(th) - cos(th)) / sigma,  sigma = sin^2 + mu sin cos + m cos = S1/r^2,  m = (y-k)/r."""
import numpy as np
from base import load
D = load()
worst = 0
for key in sorted(D):
    P = D[key]; mu = P['mu']; k = P['k']
    U, p, y = P['U'], P['p'], P['y']
    pd = k - mu*p - y
    r = np.hypot(U, p); s = p/r; c = U/r; m = (y-k)/r
    sig = s*s + mu*s*c + m*c
    lhs = ((U*p + p*pd)/r**2) / (-(p*p - U*pd)/r**2)          # (r'/r)/theta'
    rhs = s*(m + mu*s - c)/sig
    e1 = np.max(np.abs(sig - (p*p-U*pd)/r**2)); e2 = np.max(np.abs(lhs-rhs)/(1+np.abs(lhs)))
    worst = max(worst, e1, e2)
print("max discrepancy of the polar orbit equation over 55 cycles:", worst)
