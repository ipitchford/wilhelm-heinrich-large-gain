#!/usr/bin/env python3
"""Figure 1 of the paper: an illustration, not evidence.

Left: projection to the (U,p) plane of one orbit on the unstable manifold of E (grey) spiralling out to the
cycle (black) at mu = 1, rho = 1.3.  Right: the three concentrations over one period, from the X crest.
Floating point (DOP853).  Self-contained: needs only numpy, scipy and matplotlib.
Run:  python make_fig_attractor.py   (writes fig_attractor.pdf next to this file)
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MU, RHO = 1.0, 1.3
K = RHO * MU * (1 + MU)
LAM = K / MU
INK, MUTED, GRID = "#0b0b0b", "#8a8984", "#d9d8d3"
SERIES = {"X": "#2a78d6", "Z": "#eb6834", "Y": "#1baf7a"}      # validated categorical slots 1-3
STYLE = {"X": "-", "Z": (0, (5, 2)), "Y": (0, (5, 1.5, 1, 1.5))}  # secondary encoding for print


def rhs(t, z):
    U, p, y = z
    return [p, K - MU * p - y, K * np.exp(U) - y]


def crest(t, z):
    return z[1]
crest.direction = -1

J = np.array([[0, 1, 0], [0, -MU, -1], [K, 0, -1.0]])
w, V = np.linalg.eig(J)
i = int(np.argmax(w.imag))
v = np.real(V[:, i]); v /= np.linalg.norm(v)
z0 = np.array([0.0, 0.0, K]) + 2e-2 * v
t_end = (np.log(1 / 2e-2) + 14) / w[i].real
sol = solve_ivp(rhs, [0, t_end], z0, method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True, events=crest)
tc = sol.t_events[0]
T = tc[-1] - tc[-2]
ts = np.linspace(0, tc[-2], 12000)
spiral = sol.sol(ts)
tt = np.linspace(tc[-2], tc[-1], 4000)
U, p, y = sol.sol(tt)
phase = (tt - tt[0]) / T
X, Z, Y = np.exp(U), y / K, 1 - p / LAM

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.linewidth": 0.6, "pdf.fonttype": 42})
fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.9), gridspec_kw={"width_ratios": [1.05, 1]})

a.plot(spiral[0], spiral[1], color=MUTED, lw=0.55, zorder=1)
a.plot(U, p, color=INK, lw=1.5, zorder=3)
a.axhline(0, color=GRID, lw=0.6, zorder=0); a.axvline(0, color=GRID, lw=0.6, zorder=0)
iu, il = int(np.argmax(p)), int(np.argmin(p))
for j, d in ((iu, 1), (il, 1)):
    a.annotate("", xy=(U[j + 25 * d], p[j + 25 * d]), xytext=(U[j], p[j]),
               arrowprops=dict(arrowstyle="-|>", color=INK, lw=0, mutation_scale=9), zorder=4)
a.plot([0], [0], "o", ms=3.2, color=INK, zorder=5)
a.plot([U.max(), U.min()], [0, 0], "o", ms=4.2, mfc="white", mec=INK, mew=1.0, zorder=5)
a.annotate("$E$", (0, 0), xytext=(5, 5), textcoords="offset points", color=INK)
a.annotate("crest", (U.max(), 0), xytext=(-4, 12), textcoords="offset points", ha="right", color=INK,
           bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none"))
a.annotate("trough", (U.min(), 0), xytext=(6, -12), textcoords="offset points", color=INK)
a.text(0.42, 0.93, "rising arc, $p>0$", transform=a.transAxes, ha="center", color=INK)
a.text(0.42, 0.04, "falling arc, $p<0$", transform=a.transAxes, ha="center", color=INK)
a.set_xlabel("$U=\\log(x/k)$"); a.set_ylabel("$p=U'$")
a.set_ylim(p.min() * 1.32, p.max() * 1.32)
for s in ("top", "right"):
    a.spines[s].set_visible(False)

for name, ser in (("X", X), ("Z", Z), ("Y", Y)):
    b.plot(phase, ser, color=SERIES[name], ls=STYLE[name], lw=1.5, label=f"${name}$")
    jmax, jmin = int(np.argmax(ser)), int(np.argmin(ser))
    b.plot(phase[jmax], ser[jmax], "o", ms=4.2, color=SERIES[name], mec="white", mew=0.8, zorder=5)
    b.plot(phase[jmin], ser[jmin], "o", ms=4.2, mfc="white", mec=SERIES[name], mew=1.1, zorder=5)
for name, ser, ph, dx, dy in (("X", X, 0.22, -15, -4), ("Z", Z, 0.52, -15, -11), ("Y", Y, 0.56, 6, 6)):
    j = int(np.argmin(np.abs(phase - ph)))
    b.annotate(f"${name}$", (phase[j], ser[j]), xytext=(dx, dy), textcoords="offset points", color=INK)
b.axhline(1, color=GRID, lw=0.6, zorder=0)
b.set_yscale("log"); b.set_xlim(0, 1)
b.set_xlabel("time from the $X$ crest, in periods"); b.set_ylabel("concentration (log scale)")
b.legend(frameon=False, loc="lower right", handlelength=2.6, labelcolor=INK)
for s in ("top", "right"):
    b.spines[s].set_visible(False)
fig.tight_layout(w_pad=1.6)
out = Path(__file__).with_name("fig_attractor.pdf")
fig.savefig(out); fig.savefig(out.with_suffix(".png"), dpi=170)
order = sorted([(phase[int(np.argmax(s))] % 1, n + "max") for n, s in (("X", X), ("Z", Z), ("Y", Y))] +
               [(phase[int(np.argmin(s))] % 1, n + "min") for n, s in (("X", X), ("Z", Z), ("Y", Y))])
print("T =", round(T, 6), " U range", round(U.min(), 4), round(U.max(), 4), " order from X crest:", [n for _, n in order])
