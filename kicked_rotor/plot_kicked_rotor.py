
#!/usr/bin/env python3
"""
Visualization & analysis for kicked_rotor.f90 output.

Usage:
    ./kicked_rotor input.dat                    # writes quantum.txt, classical.txt,
                                               #   momentum.txt, map.txt
    python plot_kicked_rotor.py diffusion       # THE figure: <p^2> quantum vs classical
    python plot_kicked_rotor.py momentum        # |psi_m|^2 + exponential-tail fit
    python plot_kicked_rotor.py map             # standard-map phase portrait
    python plot_kicked_rotor.py norm            # norm drift (unitarity check)

File formats (headers start with '#'):
    quantum.txt   : n p2 norm
    classical.txt  : n p2_avg
    momentum.txt  : m prob
    map.txt       : n x p
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# keep in sync with the namelist defaults in kicked_rotor.f90
K, HBAR = 5.0, 1.0
D = K**2 / 2          # classical diffusion coefficient (p^2 per kick)

# Reference values for K = 5 (measured with this code, hbar = 1, N = 128):
#   <p^2>_sat ~ 100-140 (fluctuates with the initial state / seed)
#   break time  t* = <p^2>_sat / D ~ 8-12 kicks
#   localization length from the exponential fit: l ~ D/hbar^2 within a
#   factor ~ 2 (fluctuates: single-exponential fits are window-dependent)


# ---------------------------------------------------------------- diffusion
def plot_diffusion():
    q = np.loadtxt("quantum.txt")
    nq, p2q = q[:, 0], q[:, 1]
    norm_drift = np.abs(q[:, 2] - 1.0).max()
    print(f"max |norm - 1| over the run: {norm_drift:.3e}  (should be ~1e-15)")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(nq, p2q, lw=1.2, label="quantum $\\langle p^2\\rangle$")
    if norm_drift > 1e-8:
        print("WARNING: norm not conserved — check the DFT / phases!")

    try:
        c = np.loadtxt("classical.txt")
        ax.plot(c[:, 0], c[:, 1], lw=1.2, alpha=0.7,
                label="classical ensemble $\\langle p^2\\rangle$")
    except OSError:
        pass
    nmax = nq.max()
    ax.plot(nq, D * nq, "k--", lw=1,
            label=f"classical diffusion $Dn$, $D=K^2/2={D}$")

    # break time: t* = <p^2>_sat / D, saturation = mean of the last quarter
    psat = p2q[len(p2q) // 4 * 3:].mean()
    tstar = psat / D
    print(f"<p^2>_sat = {psat:.1f}   t* = <p^2>_sat/D = {tstar:.1f} kicks")
    print(f"theory: t* ~ D/(2 hbar^2) = {D/(2*HBAR**2):.1f}  "
          "(order of magnitude: prefactor fluctuates)")

    ax.axhline(psat, color="C0", ls=":", lw=1)
    ax.set_xlabel("kick number $n$")
    ax.set_ylabel("$\\langle p^2\\rangle$")
    ax.set_title("Dynamical localization: quantum diffusion stops, classical does not")
    ax.legend()
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- momentum
def plot_momentum():
    data = np.loadtxt("momentum.txt")
    m, prob = data[:, 0], data[:, 1]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogy(m[np.abs(m) > 0], prob[np.abs(m) > 0], ".", ms=3,
                label="$|\\psi_m|^2$")
    # exponential-tail fit on the positive side: P(m) ~ exp(-2 m / l)
    lhat = D / HBAR**2
    lo, hi = max(2, int(lhat / 2)), int(3 * lhat)
    sel = (m >= lo) & (m <= hi) & (prob > 0)
    if sel.sum() >= 5:
        slope = np.polyfit(m[sel], np.log(prob[sel]), 1)[0]
        ell = -2.0 / slope
        mm = np.linspace(lo, hi, 50)
        ax.semilogy(mm, np.exp(slope * mm + np.polyfit(m[sel], np.log(prob[sel]), 1)[1]),
                    "r--", lw=1.5,
                    label=f"fit $e^{{-2|m|/\\ell}}$, $\\ell={ell:.1f}$")
        print(f"localization length (fit window m in [{lo},{hi}]): "
              f"l = {ell:.1f},  theory ~ D/hbar^2 = {lhat:.1f}")
        print("(expect agreement only within a factor ~2 — discuss why)")
    ax.set_xlabel("$m$   (p = m hbar)")
    ax.set_ylabel("$|\\psi_m|^2$")
    ax.set_title("Momentum distribution: exponentially localized")
    ax.legend()
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- map
def plot_map():
    data = np.loadtxt("map.txt")
    n, x, p = data[:, 0], data[:, 1], data[:, 2]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.plot(x, p, ".", ms=1.5, alpha=0.6)
    ax.set_xlabel("$x$ (wrapped)")
    ax.set_ylabel("$p$")
    ax.set_title(f"Standard map, K = {K}: stroboscopic phase portrait")
    fig.tight_layout()
    plt.show()


# ---------------------------------------------------------------- norm
def plot_norm():
    q = np.loadtxt("quantum.txt")
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.semilogy(q[:, 0], np.abs(q[:, 2] - 1.0), lw=1)
    ax.set_xlabel("kick number $n$")
    ax.set_ylabel("|norm - 1|")
    ax.set_title("Unitarity check (exact split-step: machine precision)")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "diffusion"
    if mode == "diffusion":
        plot_diffusion()
    elif mode == "momentum":
        plot_momentum()
    elif mode == "map":
        plot_map()
    elif mode == "norm":
        plot_norm()
    else:
        print(__doc__)
