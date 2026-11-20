#!/usr/bin/env python3
"""
Visualization & analysis for chua.f90 output.

Usage:
    ./chua input.dat                              # writes trajectory.txt, poincare.txt
    python plot_chua.py attractor                 # double-scroll projections + 3D
    python plot_chua.py timeseries                 # x(t), y(t), z(t)
    python plot_chua.py poincare                   # section at y = 0 (x vs z)
    python plot_chua.py lyapunov                    # d(t) from trajectory vs results_pert.txt
    python plot_chua.py bifurcation bifur.dat      # param vs local maxima of x

File formats (all with '#' header lines):
    trajectory.txt  : tau x y z
    results_pert.txt: tau x y z
    poincare.txt    : x z
    bifur.dat       : param x_max   (build it with the sweep script, Ex. 3.2)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# keep in sync with the namelist defaults in chua.f90
ALPHA, BETA, GAMMA = 9.0, 14.286, 0.05
A, B = -1.143, -0.714

# Literature benchmark (Matsumoto-Chua-Komuro 1985, canonical parameters
# alpha=9, beta=100/7, gamma=0, a=-8/7, b=-5/7):
#   lambda ~ (0.23, 0, -1.77), Kaplan-Yorke dimension ~ 2.13
# Renormalized two-trajectory estimates for the same system give
# lambda1 ~ 0.33 +- 0.02 (see Exercise 2.2: why do published values differ?).


# ---------------------------------------------------------------- attractor
def plot_attractor(fname="trajectory.txt", skip=0):
    data = np.loadtxt(fname)
    t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    x, y, z = x[skip:], y[skip:], z[skip:]      # drop the transient

    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(x, y, lw=0.3, alpha=0.6)
    ax1.set_xlabel("x"); ax1.set_ylabel("y")
    ax1.set_title("x-y projection")
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(x, z, lw=0.3, alpha=0.6)
    ax2.set_xlabel("x"); ax2.set_ylabel("z")
    ax2.set_title("x-z projection")
    ax3 = fig.add_subplot(1, 3, 3, projection="3d")
    ax3.plot(x, y, z, lw=0.2, alpha=0.5)
    ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z")
    ax3.set_title("double scroll")
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- timeseries
def plot_timeseries(fname="trajectory.txt"):
    data = np.loadtxt(fname)
    t, x, y, z = data.T
    fig, ax = plt.subplots(3, 1, figsize=(9, 6), sharex=True)
    for a, v, name in zip(ax, (x, y, z), "xyz"):
        a.plot(t, v, lw=0.6)
        a.set_ylabel(name)
    ax[-1].set_xlabel("tau")
    fig.suptitle("Chua oscillator: time series")
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- poincare
def plot_poincare(fname="poincare.txt"):
    data = np.loadtxt(fname)
    x, z = data[:, 0], data[:, 1]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x, z, ".", ms=1.5, alpha=0.7)
    ax.set_xlabel("x"); ax.set_ylabel("z")
    ax.set_title(r"Poincaré section at $y = 0$ (upward crossings)")
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- lyapunov
def plot_lyapunov(fa="trajectory.txt", fb="results_pert.txt"):
    """Plain two-trajectory separation (no renormalization!):
    exponential growth -> saturation. Slope of the growth window = lambda1.
    Only trustworthy if results_pert.txt was sampled at the SAME times."""
    A = np.loadtxt(fa); Bv = np.loadtxt(fb)
    n = min(len(A), len(Bv))
    t = A[:n, 0]
    d = np.sqrt(((A[:n, 1:] - Bv[:n, 1:])**2).sum(axis=1))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(t, d, lw=0.8)
    ax.set_xlabel(r"$\tau$"); ax.set_ylabel("phase-space distance $d(\\tau)$")
    ax.set_title("trajectory separation (growth window slope = lambda1)")
    mask = (d > 100 * d[0]) & (d < 0.5 * d.max())
    if mask.sum() > 5:
        lam = np.polyfit(t[mask], np.log(d[mask]), 1)[0]
        print(f"two-trajectory estimate: lambda1 ~ {lam:.3f}"
              "  (saturates once d ~ attractor size)")
        ax.semilogy(t[mask], np.exp(np.polyval(
            np.polyfit(t[mask], np.log(d[mask]), 1), t[mask])), "r--", lw=1.5,
            label=f"slope = {lam:.2f}")
        ax.legend()
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- bifurcation
def plot_bifurcation(fname="bifur.dat"):
    data = np.loadtxt(fname)
    p, xmax = data[:, 0], data[:, 1]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(p, xmax, ".", ms=0.8, alpha=0.6)
    ax.set_xlabel("parameter (e.g. alpha)")
    ax.set_ylabel("local maxima of x")
    ax.set_title("bifurcation diagram (period doubling -> chaos)")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "attractor"
    if mode == "attractor":
        plot_attractor(*(sys.argv[2:3] or []))
    elif mode == "timeseries":
        plot_timeseries(*(sys.argv[2:3] or []))
    elif mode == "poincare":
        plot_poincare(*(sys.argv[2:3] or []))
    elif mode == "lyapunov":
        plot_lyapunov(*(sys.argv[2:4] or []))
    elif mode == "bifurcation":
        plot_bifurcation(*(sys.argv[2:3] or []))
    else:
        print(__doc__)
