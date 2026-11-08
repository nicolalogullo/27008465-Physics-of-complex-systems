#!/usr/bin/env python3
"""
Visualization & analysis for double_pendulum.f90 output.

Usage:
    ./double_pendulum                      # writes results.dat
    python plot_double_pendulum.py angles       # theta1(t), theta2(t)
    python plot_double_pendulum.py energy       # E(t) drift (RK4 convergence)
    python plot_double_pendulum.py phase        # phase portraits
    python plot_double_pendulum.py poincare     # Poincare section (th1 crossing 0)
    python plot_double_pendulum.py lyapunov       # results.dat vs results_pert.dat
    python plot_double_pendulum.py animate

Expected file format (one row per step):  t th1 th2 w1 w2
If a 6th column is present it is assumed to be the energy E(t).
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

G = 9.81
M1, M2 = 0.01, 1.0   # keep in sync with the namelist defaults in double_pendulum.f90
L1, L2 = 1.0, 1.0


def load(fname):
    data = np.loadtxt(fname)
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]


def energy(th1, th2, w1, w2):
    """Total energy (theory.md, Ex 1.1)."""
    v = (0.5 * M1 * L1**2 * w1**2
         + 0.5 * M2 * (L1**2 * w1**2 + L2**2 * w2**2
                       + 2 * L1 * L2 * w1 * w2 * np.cos(th2 - th1))
         - (M1 + M2) * G * L1 * np.cos(th1)
         - M2 * G * L2 * np.cos(th2))
    return v


def positions(th1, th2):
    x1, y1 = L1 * np.sin(th1), -L1 * np.cos(th1)
    x2, y2 = x1 + L2 * np.sin(th2), y1 - L2 * np.cos(th2)
    return x1, y1, x2, y2


# ---------------------------------------------------------------- angles
def plot_angles(fname="results.dat"):
    t, th1, th2, w1, w2 = load(fname)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(t, th1, lw=0.8, label=r"$\theta_1$")
    ax.plot(t, th2, lw=0.8, alpha=0.8, label=r"$\theta_2$")
    ax.set_xlabel("t"); ax.set_ylabel("angle [rad]")
    ax.set_title("Double pendulum: angles vs time")
    ax.legend(); fig.tight_layout(); plt.show()


# ---------------------------------------------------------------- energy
def plot_energy(files=("results_dt0.01.dat", "results_dt0.005.dat",
                       "results_dt0.001.dat", "results_dt0.0005.dat")):
    """Relative energy drift vs t for several dt; slope ~ -4 confirms RK4."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4))
    dts = []
    for f in files:
        try:
            t, th1, th2, w1, w2 = load(f)
        except OSError:
            continue
        E = energy(th1, th2, w1, w2)
        drift = np.abs(E - E[0]) / np.abs(E[0])
        dt = np.median(np.diff(t))
        dts.append((dt, drift[-1]))
        a.loglog(t[1:], drift[1:], lw=1, label=f"dt={dt:g}")
        b.plot(t, E, lw=0.8)
    a.set_xlabel("t"); a.set_ylabel(r"$|E(t)-E(0)|/E(0)$")
    a.set_title("Energy drift"); a.legend()
    b.set_xlabel("t"); b.set_ylabel("E(t)"); b.set_title("Total energy")
    if len(dts) >= 2:
        d = np.array(dts)
        slope = np.polyfit(np.log(d[:, 0]), np.log(d[:, 1]), 1)[0]
        print(f"convergence slope: {slope:.2f}  (RK4 expectation: -4)")
    fig.tight_layout(); plt.show()


# ---------------------------------------------------------------- phase
def plot_phase(fname="results.dat"):
    t, th1, th2, w1, w2 = load(fname)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for a, th, w, name in ((ax[0], th1, w1, "1"), (ax[1], th2, w2, "2")):
        a.plot(np.mod(th + np.pi, 2 * np.pi) - np.pi, w, ".", ms=0.5,
               alpha=0.4, rasterized=True)
        a.set_xlabel(rf"$\theta_{name}$ (wrapped)")
        a.set_title(rf"phase portrait: $\omega_{name}$ vs $\theta_{name}$")
    ax[0].set_ylabel(r"$\omega$ [rad/s]")
    fig.tight_layout(); plt.show()


# ---------------------------------------------------------------- Poincare
def plot_poincare(fname="results.dat"):
    """Section: states where theta1 crosses 0 upward."""
    t, th1, th2, w1, w2 = load(fname)
    s = th1[:-1] < 0
    crossing = (th1[:-1] < 0) & (th1[1:] >= 0)
    th2c, w2c = th2[:-1][crossing], w2[:-1][crossing]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(np.mod(th2c + np.pi, 2 * np.pi) - np.pi, w2c, ".", ms=2)
    ax.set_xlabel(r"$\theta_2$"); ax.set_ylabel(r"$\omega_2$")
    ax.set_title(r"Poincaré section at $\theta_1 = 0$")
    fig.tight_layout(); plt.show()


# ---------------------------------------------------------------- Lyapunov
def plot_lyapunov(fa="results.dat", fb="results_pert.dat"):
    ta, a1, a2, a3, a4 = load(fa)
    tb, b1, b2, b3, b4 = load(fb)
    n = min(len(ta), len(tb))
    d = np.sqrt((a1[:n]-b1[:n])**2 + (a2[:n]-b2[:n])**2
                + (a3[:n]-b3[:n])**2 + (a4[:n]-b4[:n])**2)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(ta[:n], d, lw=0.8)
    ax.set_xlabel("t"); ax.set_ylabel(r"phase-space distance $d(t)$")
    ax.set_title("Trajectory separation (slope in the linear regime = λ)")
    # crude automatic slope estimate on the growth window
    mask = (d > 10*d[0]) & (d < 0.5*d.max()) & (d > 0)
    if mask.sum() > 10:
        lam = np.polyfit(ta[:n][mask], np.log(d[mask]), 1)[0]
        print(f"largest Lyapunov estimate: λ ≈ {lam:.3f} 1/s")
        ax.semilogy(ta[:n][mask], d[mask], "r-", lw=0,
                    label=f"fit region")  # highlight
        ax.semilogy(ta[:n][mask], np.exp(np.polyval(
            np.polyfit(ta[:n][mask], np.log(d[mask]), 1), ta[:n][mask])),
            "r--", lw=1.5, label=rf"$\lambda \approx {lam:.2f}$")
        ax.legend()
    fig.tight_layout(); plt.show()


# ---------------------------------------------------------------- animate
def animate(fname="results.dat", stride=5):
    from matplotlib.animation import FuncAnimation
    t, th1, th2, w1, w2 = load(fname)
    x1, y1, x2, y2 = positions(th1, th2)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(-(L1+L2)-0.2, (L1+L2)+0.2)
    ax.set_ylim(-(L1+L2)-0.2, (L1+L2)+0.2)
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    trace, = ax.plot([], [], "r-", lw=0.7, alpha=0.5)
    rod, = ax.plot([], [], "o-", lw=2, ms=8, color="C0")

    def update(i):
        j = i * stride
        trace.set_data(x2[:j], y2[:j])
        rod.set_data([0, x1[j], x2[j]], [0, y1[j], y2[j]])
        return rod, trace

    ani = FuncAnimation(fig, update, frames=len(t)//stride,
                        interval=20, blit=True)
    plt.show()
    return ani


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "angles"
    if mode == "angles":
        plot_angles(*(sys.argv[2:3] or []))
    elif mode == "energy":
        plot_energy()
    elif mode == "phase":
        plot_phase(*(sys.argv[2:3] or []))
    elif mode == "poincare":
        plot_poincare(*(sys.argv[2:3] or []))
    elif mode == "lyapunov":
        plot_lyapunov(*(sys.argv[2:4] or []))
    elif mode == "animate":
        animate(*(sys.argv[2:3] or []))
    else:
        print(__doc__)
