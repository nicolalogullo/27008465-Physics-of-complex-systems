# Double Pendulum — Exercises

The Fortran program `double_pendulum.f90` integrates the equations of motion with a 4th-order Runge–Kutta scheme. It reads all parameters from an input file (default `input.dat`), organized as **three namelists**, in this order:

```fortran
&pendulum1                  ! first arm:  m1, l1, th1, w1
&pendulum2                  ! second arm: m2, l2, th2, w2
&simulation                 ! dt, steps, poincare, perturb
```

and writes `results.dat` with one row per time step:

```
# t  th1  th2  w1  w2  E
```

The two extra options in `&simulation`:
- `poincare = .true.` — write only the states where $\theta_1$ crosses zero upward;
- `perturb = δ` — additionally integrate a second trajectory with $\theta_2 \to \theta_2 + \delta$ into `results_pert.dat`.

Exercises are grouped by theme; each builds on the previous one. Python plotting code is provided in `plot_double_pendulum.py`.

---

## Part 0 — Getting started

**Exercise 0.1 (Build and run).**
Compile and run the simulation:

```bash
gfortran -O2 -o double_pendulum double_pendulum.f90
./double_pendulum input.dat
```

Write your first `input.dat` (see the header comment at the top of the source for the exact format). Plot $\theta_1(t)$, $\theta_2(t)$ with the provided Python script. Describe qualitatively the difference between the two angles' behavior.

**Exercise 0.2 (Read the code).**
Identify in `double_pendulum.f90`:
1. where the first-order system (§4 of the theory) is implemented (`rhs`);
2. where the RK4 update (Eq. for the 4 slopes $k_1,\dots,k_4$) is implemented (`rk4_step`);
3. how the input file is read (the three namelists in the main program) and how the values are transferred into the module variable `P`.

**Exercise 0.3 (Learning to love namelists).**
Experiment with the input file and discuss:
1. what happens when a key is **misspelled** (try `step = 50000`): why is the loud runtime error a good thing here?
2. what happens when a key is **omitted**: it silently keeps its default value — convenient or dangerous? How would you check, from the program output alone, that the parameters actually used are the ones you intended?
3. why the three groups must appear **in a fixed order** in the file;
4. what happens when the file is missing or a group is malformed.
Argue, in a few sentences, which of these behaviors you would change for a scientific code used by other people.

---

## Part 1 — Validation: does the integrator tell the truth?

**Exercise 1.1 (Energy conservation).**
The double pendulum is a conservative system, so the total energy

$$
E = \tfrac12 m_1 L_1^2\dot\theta_1^2 + \tfrac12 m_2\left(L_1^2\dot\theta_1^2 + L_2^2\dot\theta_2^2 + 2L_1L_2\dot\theta_1\dot\theta_2\cos(\theta_2-\theta_1)\right) - (m_1+m_2)gL_1\cos\theta_1 - m_2gL_2\cos\theta_2
$$

must be constant in time. The reference code already writes $E(t)$ as a 6th column — but do not trust it blindly: **recompute $E$ independently in Python** from $(t, \theta_1, \theta_2, \omega_1, \omega_2)$ and check that your values coincide with the file's. Then:
1. Plot $E(t)$ vs $t$ and $|E(t)-E(0)|/E(0)$. What is the relative energy drift over the full simulation?
2. Repeat for `dt = 0.01, 0.005, 0.001, 0.0005` and plot the drift vs `dt` on a log–log scale. RK4 is 4th order: does the slope match $-4$?

**Exercise 1.2 (Small-amplitude limit: normal modes).**
Set $m_1=m_2=m$, $L_1=L_2=L$, and initial conditions $\theta_1=\theta_2=\epsilon$ (in-phase) and $\theta_1=\epsilon,\ \theta_2=-\epsilon$ (anti-phase) with small $\epsilon \approx 0.01$ rad. In the small-angle limit the system decouples into two normal modes with frequencies

$$
\omega_\pm = \omega_0\sqrt{2 \pm \sqrt{2}}, \qquad \omega_0=\sqrt{g/L}.
$$

Measure the period of $\theta_1(t)$ in both runs and compare with the prediction. This validates the equations of motion themselves, not just the integrator.

**Exercise 1.3 (Add a new output quantity — Fortran practice).**
Modify `double_pendulum.f90` to add a 7th column with the angular momentum of the second bob,
$$p_{\theta_2} = m_2 L_2^2\,\dot\theta_2 + m_2 L_1 L_2\,\dot\theta_1\cos(\theta_2-\theta_1),$$
computed in a new module function alongside `total_energy`. Convince yourself it is **not** conserved (unlike $E$), plot it over time, and explain physically why momentum exchange with the first arm breaks its conservation.

---

## Part 2 — Chaos: sensitive dependence on initial conditions

**Exercise 2.1 (The butterfly effect).**
Run a simulation with `perturb = 1.0e-3` in `&simulation`: the program produces both trajectories (`results.dat` and `results_pert.dat`) in one run. Plot $\theta_2(t)$ for both on the same axes. Estimate visually the time $t^{\ast}$ at which the two trajectories become clearly distinct, and repeat with $\delta = 10^{-6}$: roughly how much larger is $t^{\ast}$? What does this imply about long-term prediction?

**Exercise 2.2 (Divergence rate and the Lyapunov exponent).**
Define the phase-space distance

$$
d(t) = \sqrt{(\Delta\theta_1)^2 + (\Delta\theta_2)^2 + (\Delta\omega_1)^2 + (\Delta\omega_2)^2}
$$

between the two trajectories from `results.dat` and `results_pert.dat` (the Python script's `lyapunov` mode does this automatically). Plot $\log d(t)$ vs $t$: you should see (i) a transient, (ii) an exponential growth regime $d(t)\sim d(0)e^{\lambda t}$, (iii) saturation when the trajectory fills the accessible energy shell. Estimate the largest Lyapunov exponent $\lambda$ from the slope of region (ii). Repeat for a **low-energy** initial condition (small angles) — what do you find, and why?

**Exercise 2.3 (Regularity map).**
Systematically scan the initial angle: run simulations for $\theta_2(0) \in [-\pi, \pi]$ (e.g. 100 values) with $\omega_1=\omega_2=0$, $\theta_1=0$, all at the same energy scale — you will need a small script (shell, Python, or Makefile) that rewrites `th2` in the input file, runs the program, and collects the result. For each run compute a simple "chaos indicator": the energy drift is *not* a good one — use instead the maximal Lyapunov estimate from Ex. 2.2, or the RMS of $\theta_2$'s high-frequency content. Plot the indicator vs $\theta_2(0)$ and identify the intervals of regular vs chaotic initial conditions. The boundary regions correspond to the chaotic separatrix discussed in the theory.

---

## Part 3 — Phase space and Poincaré sections

**Exercise 3.1 (Phase portraits).**
Plot $\omega_1$ vs $\theta_1$ and $\omega_2$ vs $\theta_2$ (use $\theta \to \theta \bmod 2\pi$ for the angles) for:
1. a low-energy run (small angles) — quasi-periodic curves on invariant tori;
2. a high-energy run (large initial $\omega$) — chaotic sea.

**Exercise 3.2 (Poincaré section).**
Since the energy surface is 3D, a 2D section reveals the structure. Set `poincare = .true.` in `&simulation`: the program records the state **only when $\theta_1$ crosses zero** upward (this is what `results.dat` then contains — use large `steps`, e.g. $10^6$–$10^7$, to populate the section). Plot $(\theta_2, \omega_2)$ at those crossings, superposing sections from ~10 different initial conditions at the same energy:
- low energy: smooth closed curves (KAM tori);
- intermediate energy: coexisting islands of regular motion and chaotic layers;
- high energy: mostly scattered points (chaotic sea).

This is the most instructive figure of the whole set — budget time for it.

**Exercise 3.3 (Energy as the control parameter).**
Repeat the Poincaré section for $E/(m g L) \approx 0.5,\ 1,\ 2,\ 4$ (tune the initial velocities to hit each energy) and describe the order→chaos transition: at which energy scale does a macroscopic chaotic component first appear?

---

## Part 4 — Going further (optional / project-level)

- **4.1 Momentum–energy coordinates.** Replot sections in terms of the conserved quantities: $p_{\theta_2}$ (angular momentum is *not* conserved; verify this numerically by plotting $p_{\theta_2}$ over time) and $E$.
- **4.2 Fractal basin / final-state map.** For fixed energy, color-code each initial condition by some late-time observable (e.g. sign of $\theta_2(t_{end})$) and produce a 2D map: near the chaotic boundary the pattern is fractal.
- **4.3 Animation.** Use the provided Python script to render an animation of the pendulum (positions $x_i, y_i$ from the angles), and a "double-pendulum art" figure obtained by fading the trace of $m_2$.
- **4.4 Stroboscopic vs driving.** Replace RK4 with a symplectic integrator (e.g. leapfrog/velocity-Verlet adapted to this non-separable system is tricky — discuss why, and test the energy drift of RK4 vs. an implicit midpoint rule).
- **4.5 Literature search — the double pendulum in the wild.**
  Find **one real system in your own field** (physics, engineering, biology, economics, or anything else) that has been modeled in the scientific literature as a double pendulum or as two coupled nonlinear oscillators. Good starting points for a search: "double pendulum model" combined with *robot arm / crane / biomechanics / walking / energy harvesting / machine learning benchmark / SINDy*.
  Write a half-page report containing:
  1. the reference (paper, book, or thesis) and the system studied;
  2. what plays the role of $m_1, m_2, L_1, L_2$ in that system, and what the two angles represent physically;
  3. one simplification the authors made compared to the ideal double pendulum (damping, driving, small-angle limit, …) and why it was justified;
  4. whether chaos plays a role in their analysis (e.g. do they estimate a Lyapunov exponent, as in Ex. 2.2?).
  Be ready to present your findings in 3 minutes: this exercise connects the numerical methods of this section to how double-pendulum-type models are actually used in research.

