# Kicked Rotor — Exercises

The Fortran program `kicked_rotor.f90` evolves both the **classical** standard map and the **quantum** kicked rotor. It reads three namelists from `input.dat`:

```fortran
&rotor      ! K, hbar, N          (kick strength, effective hbar, basis size)
&initial    ! m0, x0, p0, seed    (quantum initial state |m0>, classical IC, RNG seed)
&simulation ! nkicks, stride, ensemble, quantum, classical
```

and writes `quantum.txt` (`# n p2 norm`), `momentum.txt` (`# m prob`), `classical.txt` (`# n p2_avg`), and `map.txt` (`# n x p`). The quantum evolution is an **exact** split-step Floquet iteration — no integrator error. Defaults ($K=5$, $\hbar=1$, $N=128$, 2000 kicks) reproduce the theory figure.

Python plotting is in `plot_kicked_rotor.py` (modes: `diffusion`, `momentum`, `map`, `norm`).

---

## Part 0 — Getting started

**Exercise 0.1 (Build and run).**
Write `input.dat`, compile with `gfortran -O2 -o kicked_rotor kicked_rotor.f90`, run, then:

```bash
python plot_kicked_rotor.py diffusion
python plot_kicked_rotor.py map
```

You should see: (a) the classical ensemble diffusing linearly, the quantum $\langle p^2\rangle$ following it briefly and then flattening; (b) the phase portrait of the standard map. This one figure — *classical diffusion versus quantum localization* — is the summary of the whole section.

**Exercise 0.2 (Read the code).**
Find where: (1) the standard map is iterated; (2) the DFT matrix $F$ is built and verify from the code that $F^\dagger F = I$ (why is unitarity essential here?); (3) the free-evolution phase $e^{-im^2\hbar/2}$ and kick phase $e^{-iK\cos x_j/\hbar}$ are applied. Why is the split-step **exact**, while RK4 in the pendulum/Chua sections had a step-size error? (Hint: what does a delta-kick do to the ordering of the two exponential terms?)

**Exercise 0.3 (Effective $\hbar$).**
In this dimensionless formulation $\hbar$ is a free parameter. What physical quantities of the real experiment (kick period $T$, rotor moment of inertia $I$, laser wavelength) set it, and what do the limits $\hbar\to0$ and $\hbar\to\infty$ correspond to physically?

---

## Part 1 — Validation

**Exercise 1.1 (Unitarity: the built-in check).**
The third column of `quantum.txt` is $\|\psi_n\|^2$. Plot it with `plot_kicked_rotor.py norm` and state the drift: it should be at **machine precision** ($\sim10^{-15}$), not at RK4 level ($\sim10^{-8}$). Compare with the energy drift of the pendulum exercises and explain the difference in terms of the structure of the two algorithms (unitary matrix products vs. approximate time integration). Then break it on purpose: modify a kick phase by hand (remove the minus sign in one exponent) — and check that the norm still stays 1. Why? What would you have to break to violate the norm?

**Exercise 1.2 (The classical limit: one kick at a time).**
Starting from $p_0=0$, compute by hand $\langle p^2\rangle$ after one kick for a uniform ensemble of $x_0$ (answer: $K^2/2$). Check that both `classical.txt` (first row after one kick) and `quantum.txt` reproduce it. This is the cleanest one-line validation that the quantum and classical codes describe the same physics.

**Exercise 1.3 (Area preservation).**
Compute the Jacobian determinant of the standard map analytically and show it equals 1 for all $(x,p)$. Verify numerically: iterate the corners of a small square of initial conditions (say $10^{-3}$-sized) for 100 kicks and estimate the area of the resulting cloud (e.g. by the determinant of the covariance matrix times $\pi$). Compare with the pendulum section: which system there played the same "area preservation" role, and in which variables?

---

## Part 2 — Classical chaos and diffusion

**Exercise 2.1 (Phase portraits: the KAM scenario again).**
Produce the standard-map portrait (`map.txt`, wrapped $x\in[-\pi,\pi)$) for $K = 0.5,\ 1.0,\ 1.5,\ 5.0$. Describe the transition: invariant circles and islands at small $K$, the appearance of a connected chaotic sea near $K_c\approx1$, global chaos at large $K$. Compare side by side with the double pendulum Poincaré sections of the first section: same phenomenology in a discrete-time system, at a fraction of the numerical cost.

**Exercise 2.2 (Diffusion coefficient).**
With $K=5$ and `ensemble = 5000`, fit $\langle p^2\rangle$ vs $n$ over $n\in[10,100]$ and extract the slope. Compare with $D=K^2/2=12.5$. Repeat for $K=3$: is the deviation from $K^2/2$ larger? Why? (Hint: "effectively uncorrelated kicks" requires $K\gg1$; at moderate $K$ correlations $\langle\sin x_n\sin x_m\rangle\neq0$ renormalize $D$ — the same physics as the decorrelation assumption behind the Lyapunov estimate in the pendulum.)

**Exercise 2.3 (Lyapunov exponent, analytically and numerically).**
For $K\gg1$ the map locally stretches by $\approx K/2$ per kick, giving $\lambda\approx\ln(K/2)$. Verify numerically for $K=5$ using two classical trajectories separated by $10^{-9}$ in initial $x$ (add a `perturb`-style option to the Fortran code, or post-process two `map.txt` runs). You should find $\lambda \approx 0.9$: positive Lyapunov exponent, yet the *quantum* evolution of the same $K$ shows no exponential growth of anything. Reconcile the two statements — this is the exercise the whole section exists for.

---

## Part 3 — Dynamical localization

**Exercise 3.1 (The break time).**
With the defaults ($K=5$, $\hbar=1$): from `quantum.txt` measure the saturation value $\langle p^2\rangle_\mathrm{sat}$ (average of the last quarter of the run) and define $t^\ast = \langle p^2\rangle_\mathrm{sat}/D$. Expected: $\langle p^2\rangle_\mathrm{sat}\approx100$–$140$, $t^\ast\approx8$–$12$ kicks (the range is genuine run-to-run fluctuation, not a bug: vary `m0` and `seed` and watch it move). Compare with the order-of-magnitude prediction $t^\ast\sim D/(2\hbar^2)$.

**Exercise 3.2 (Localization length).**
Plot `momentum.txt` on a semilog scale and fit the envelope $|\psi_m|^2\sim e^{-2|m|/\ell}$ on the positive side (the Python `momentum` mode does this; the fit window is printed). For $K=5$, $\hbar=1$ expect $\ell\approx10$–$17$, to be compared with the theory $\ell\approx D/\hbar^2=12.5$. Then repeat for $K=7$: does $\ell$ scale with $D\propto K^2$? Note honestly: the fit is window-dependent, because the distribution is exponential only in its tail. Which window choices make the comparison fair?

**Exercise 3.3 (The correspondence principle in a chaotic system).**
Repeat Exercise 3.1 with $\hbar=0.5$ (`N = 256`, keep $K=5$; run at least 150 kicks). Predict first, then measure: $t^\ast$ should grow roughly like $\hbar^{-2}$ (theory: factor 4 when halving $\hbar$). You will find a factor nearer 2 than 4. Discuss the discrepancy — prefactors of the localization length fluctuate strongly, and single-parameter scaling is only asymptotic. What would you do to establish the exponent convincingly (hint: many $\hbar$ values, logarithmic fit, and *averaging over initial states* — the same statistical honesty lesson as the Chua Lyapunov exponent).

**Exercise 3.4 (Sensitivity of localization to $\hbar$: quantum resonances).**
Try $\hbar = 2\pi$ and $\hbar = 4\pi/5$. At special rational values of $\hbar/4\pi$ the free-evolution phases rephase coherently and the rotor **ballistically accelerates**: $\langle p^2\rangle\propto n^2$ instead of saturating. Verify this "anti-localization", and explain why it disappears for generic irrational $\hbar/4\pi$.

---

## Part 4 — Going further (optional / project-level)

- **4.1 Fidelity decay.** Run two quantum evolutions differing by a tiny perturbation $\delta K$ of the kick strength and track the fidelity $F(n)=|\langle\psi_n(K)|\psi_n(K+\delta K)\rangle|^2$. In the chaotic regime $F(n)\sim e^{-\lambda_{cl} n}$ (classical Lyapunov rate) up to a $\hbar$-dependent break, then faster — the Loschmidt-echo signature of quantum chaos. Compare the observed decay rate with your Exercise 2.3 Lyapunov exponent.
- **4.2 FFT.** Replace the $O(N^2)$ DFT with a Cooley–Tukey FFT ($O(N\log N)$) and benchmark both. With the FFT, push to $N=1024$, $\hbar=0.25$, and re-test the $t^\ast\propto\hbar^{-2}$ scaling with much smaller finite-size effects.
- **4.3 Husimi sections.** Compute the Husimi distribution of $\psi_n$ (coherent states on the $(x,p)$ cylinder) at several $n$ and watch the localized packet sit *inside* the classical chaotic sea of Exercise 2.1 — the quantum-classical phase-space comparison, made direct.
- **4.4 Literature search.** As in the previous sections: find one real application of the kicked rotor / dynamical localization (cold-atom experiments of Moore–Raizen, atom-optics realization of Anderson localization, chaos-assisted tunneling, quantum simulation of localization in photonic lattices, kicked tops in quantum information) and write the standard half-page report (reference, mapping, simplification, role of chaos/localization).
- **4.5 The course in one figure.** Design a single slide comparing the three systems of the course — double pendulum (conservative classical chaos), Chua oscillator (dissipative classical chaos), kicked rotor (quantum chaos) — in terms of: what "unpredictability" means, what the conserved/checkable quantity is, what the asymptotic state of the dynamics is, and which numerical method and why.

---

## Suggested checks built into the reference code

- The **norm column** of `quantum.txt` (unitarity to machine precision).
- **One-kick check**: $\langle p^2\rangle = K^2/2$ after the first kick, identical in `quantum.txt` and `classical.txt`.
- **Diffusion line** $Dn$ with $D=K^2/2$ drawn by the `diffusion` mode, against both ensembles.
- Reference values for $K=5$, $\hbar=1$ (measured, fluctuating): $\langle p^2\rangle_\mathrm{sat}\approx100$–$140$, $t^\ast\approx8$–$12$, $\ell\approx10$–$17$ vs. $D/\hbar^2=12.5$.
