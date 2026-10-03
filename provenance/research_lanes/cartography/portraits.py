#!/usr/bin/env python3
"""S6 portraits: the four source ratios whose cycle means all equal 1 + <eta> (so lambda_c = <R_i> - A):
   R1 = x/zeta,  R2 = zeta/eta + mu,  R3 = k(U e^U - e^U + 1)/K,  R4 = x p^2 / W.
Saves portrait_mu<mu>_rho<rho>.npz/.png for three parameter sets."""
import sys
sys.path.insert(0, '/Users/admin/claude-sandbox/wilhelm-heinrich-2026-10-03/work/lanes/cartography')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cart as C
from hunt import primary_cycle
from stests import g2

SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']      # validated categorical slots 1-4 (light surface)
INK, INK2, SURF, GRID = '#0b0b0b', '#52514e', '#fcfcfb', '#e4e3df'
NAMES = ['R1 = x/ζ', 'R2 = ζ/η + μ', 'R3 = k(U·eᵁ − eᵁ + 1)/K', 'R4 = x·p²/W']

for mu, rho in [(1.0, 2.0), (0.1, 1.5), (0.01, 3.0)]:
    k = rho * C.kH(mu)
    A = 1 + mu
    cyc, mono, hb, rb = primary_cycle(mu, rho)
    pr = C.profile(mu, k, cyc['U0'], cyc['y0'], cyc['T'], nmin=8000)
    t, U, p, y, eta, zeta, K, W = (pr[n] for n in ('t', 'U', 'p', 'y', 'eta', 'zeta', 'K', 'W'))
    T = cyc['T']
    x = k * np.exp(U)
    with np.errstate(divide='ignore', invalid='ignore'):
        R = [x / zeta, zeta / eta + mu, k * g2(U) / K, x * p * p / W]
    means = [float(np.trapezoid(r, t) / T) for r in R]
    np.savez_compressed(f'portrait_mu{mu:g}_rho{rho:g}.npz', t=t, U=U, p=p, y=y, x=x, eta=eta, zeta=zeta, K=K, W=W,
                        R1=R[0], R2=R[1], R3=R[2], R4=R[3], mu=mu, rho=rho, k=k, T=T, A=A,
                        lambda_c=pr['lambda_c'], means=np.array(means))
    print(f"mu={mu} rho={rho}: T={T:.4f} lambda_c={pr['lambda_c']:+.6f}; <R_i>-A = " + ', '.join(f'{m - A:+.6f}' for m in means)
          + f"; min zeta/eta={np.min(zeta/eta):.4f}; ranges R_i-A: " + ', '.join(f'[{np.min(r)-A:+.3g},{np.max(r)-A:+.3g}]' for r in R))
    fig = plt.figure(figsize=(11, 7.2), facecolor=SURF)
    gs = fig.add_gridspec(4, 2, width_ratios=[1.45, 1], hspace=0.45, wspace=0.3)
    ph = t / T
    for i in range(4):
        ax = fig.add_subplot(gs[i, 0], facecolor=SURF)
        v = R[i] - A
        ax.axhline(0, color=INK2, lw=0.8)
        ax.plot(ph, v, color=SERIES[i], lw=2)
        ax.fill_between(ph, 0, v, where=v < 0, color=SERIES[i], alpha=0.18, lw=0)
        ax.set_yscale('symlog', linthresh=0.1 * A)
        ax.set_xlim(0, 1)
        ax.grid(True, color=GRID, lw=0.6)
        for s in ax.spines.values():
            s.set_color(GRID)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.set_title(f'{NAMES[i]}   minus A = 1+μ      (cycle mean {means[i] - A:+.4f} = λc)', fontsize=9.5, color=INK, loc='left')
        if i == 3:
            ax.set_xlabel('phase t/T (t = 0 at the X crest)', fontsize=9, color=INK2)
        else:
            ax.set_xticklabels([])
    ax = fig.add_subplot(gs[:, 1], facecolor=SURF)
    ax.plot(U, p / np.sqrt(mu), color=INK2, lw=2)
    for f_, lab in ((0.0, 'crest t=0'), (0.25, 't/T=0.25'), (0.5, '0.5'), (0.75, '0.75')):
        j = int(np.argmin(np.abs(ph - f_)))
        ax.plot(U[j], p[j] / np.sqrt(mu), 'o', ms=8, color=INK, mec=SURF, mew=2)
        ax.annotate(lab, (U[j], p[j] / np.sqrt(mu)), textcoords='offset points', xytext=(8, 6), fontsize=8.5, color=INK)
    ax.plot(0, 0, 'x', color=INK2, ms=7)
    ax.annotate('E', (0, 0), textcoords='offset points', xytext=(6, -12), fontsize=9, color=INK2)
    ax.set_xlabel('U = log(x/k)', fontsize=9, color=INK2)
    ax.set_ylabel('p/√μ', fontsize=9, color=INK2)
    ax.grid(True, color=GRID, lw=0.6)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_title('projected cycle (clockwise)', fontsize=9.5, color=INK, loc='left')
    fig.suptitle(f'Wilhelm–Heinrich cycle, μ = {mu:g}, ρ = {rho:g}:  four source ratios around the cycle;  λc = ⟨Rᵢ⟩ − A = {pr["lambda_c"]:+.5f}',
                 fontsize=11, color=INK, x=0.02, ha='left')
    fig.savefig(f'portrait_mu{mu:g}_rho{rho:g}.png', dpi=130, facecolor=SURF, bbox_inches='tight')
    plt.close(fig)
