
# Chua Oscillator — Exercises

The Fortran program `chua.f90` integrates the dimensionless Chua equations with a 4th-order Runge–Kutta scheme. It reads three namelists from `input.dat`:

```fortran
&circuit     ! alpha, beta, gamma, a, b
&initial     ! x0, y0, z0
&simulation  ! dt, steps, stride, poincare, perturb, analysis
```

and writes `trajectory.txt` (`# tau x y z`, one row every `stride` steps), `poincare.txt` (`# x z`, states where $y$ crosses 0 upward), and, if `perturb ≠ 0`, `results_pert.txt` (a second trajectory with $x_0 \to x_0 + \delta$). With `analysis = .true.` it does not integrate: it prints the equilibria and the Jacobian eigenvalues (Part 1.3). Default parameters are the canonical double-scroll set: $\alpha=9$, $\beta=100/7$, $\gamma \approx 0$, $a=-8/7$, $b=-5/7$.

Python plotting is in `plot_chua.py`. Exercises are grouped by theme; each builds on the previous one.


## Part 0 — Getting started

**Exercise 0.1 (Build and run).**
Compile, write `input.dat` (see the header comment of the source), run, and produce the double-scroll attractor:

```bash
gfortran -O2 -o chua chua.f90
./chua input.dat
python plot_chua.py attractor
```

Compare your $x$–$y$ projection with the double-scroll figures of the reference paper (Chua 1993 / Matsumoto–Chua–Komuro 1985).

**Exercise 0.2 (Read the code).**
Find in `chua.f90` where (1) the piecewise nonlinearity $f(x)$ is implemented, (2) the RK4 update happens, (3) the Poincaré crossing condition is detected. Why is the crossing condition evaluated **after** each step, using the sign of $y$ at the previous step?

**Exercise 0.3 (Compare with the double pendulum).**
The Chua oscillator is a **dissipative** system, the pendulum a **conservative** one. Predict two qualitative consequences for the numerical output (hint: what happens to a cloud of initial conditions?), then verify both by running several trajectories from different initial conditions.


## Part 1 — Validation

**Exercise 1.1 (Integrator convergence).**
The Chua system has no conserved energy, so the energy-drift test of the pendulum section does not apply. Instead: run the same initial condition with `dt = 0.01, 0.005, 0.001` over a short time window and compare $x(\tau)$ curve against curve. At which $\tau$ does the coarsest run first deviate visibly? Repeat at a longer horizon and explain why *no* step size gives you a reliable trajectory forever (this is physics, not numerics).

**Exercise 1.2 (Symmetry).**
Show analytically that the system is invariant under $(x,y,z) \to (-x,-y,-z)$. Verify numerically: run from $(0.1,0,0)$ and from $(-0.1,0,0)$, and check that the two attractors are mirror images. Why does this imply the equilibria come in symmetric pairs $P_\pm$?

**Exercise 1.3 (Equilibria and stability — the built-in check).**
Run with `analysis = .true.`. The program prints the three equilibria and their Jacobian eigenvalues, plus a self-check (sum of eigenvalues must equal $-\mathrm{tr}\,J$). Now do the theory side:
1. derive the equilibrium coordinates by solving $\dot x=\dot y=\dot z=0$:
$$y^\ast = \frac{\gamma x^\ast}{\beta+\gamma},\qquad z^\ast = -\frac{\beta x^\ast}{\beta+\gamma},\qquad x^{\ast} = \frac{(b-a)(\gamma+\beta)}{\beta + b(\gamma+\beta)}$$
(note $y^\ast, z^\ast$ are **proportional to** $x^\ast$ — beware of "reference" notes that state them as constants, which is only valid if $x^\ast = \pm 1$) — for the canonical parameters you should get exactly $x^\ast = \pm 1.5$ and hence $P_\pm = (\pm 1.5,\ 0,\ \mp 1.5)$;
2. derive the characteristic polynomial $\lambda^3 + c_2\lambda^2 + c_1\lambda + c_0$ with
$$c_2 = \gamma + 1 + \alpha(1+m),\qquad c_1 = \gamma+\beta+\alpha(1+m)(\gamma+1)-\alpha,\qquad c_0 = \alpha\big[(1+m)(\gamma+\beta)-\gamma\big],$$
where $m = f'(x^\ast) = a$ or $b$ depending on the region;
3. confirm that for the canonical parameters $P_0$ is an **index-1 saddle** (one positive real eigenvalue, one stable complex pair) while $P_\pm$ are **unstable spirals** (one stable real eigenvalue $\approx -3.9$, one unstable pair with real part $\approx +0.2$). The double scroll winds around $P_\pm$ — check on your attractor plot that the two lobes are centered on $(\pm 1.5,\, y^\ast,\, z^\ast) = (\pm 1.5, 0, \mp 1.5)$.

**Exercise 1.4 (The nonlinearity).**
Plot $f(x)$ over $x \in [-3,3]$ (either by hand or with a small program) and verify its three segments and slopes. Then show that the one-line `abs()` form used in the code is equivalent to the explicit piecewise form
$$f(x)=\begin{cases} bx + (a-b) & x>1\\ ax & |x|\le 1\\ bx - (a-b) & x<-1.\end{cases}$$
What would change physically if the slopes were both positive?


## Part 2 — Chaos

**Exercise 2.1 (The double scroll).**
With the canonical parameters, produce: (a) the three 2D projections, (b) the 3D attractor, (c) the time series $x(\tau)$. Describe in words how the trajectory alternates between the two lobes, and relate the lobe centers to Exercise 1.3.

**Exercise 2.2 (The Lyapunov exponent — with a twist).**
Run with `perturb = 1.0e-8` and use `plot_chua.py lyapunov`: you should see the three regimes (transient, exponential growth, saturation). A more accurate **renormalized** estimate (restart the perturbation at size $d_0$ every $m$ steps and accumulate $\ln(d/d_0)$) gives, for the canonical parameters,
$$\lambda_1 \approx 0.33 \pm 0.02 .$$
Some references report $\lambda_1 \approx 0.23$ for "the same" circuit. **Check, and explain, the discrepancy:** where do published Lyapunov values legitimately differ (time-step accuracy, transient handling, renormalization, actual parameter set, e.g. $\gamma=0$ vs $\gamma=0.05$)? This is the most honest lesson in numerical chaos measurement: always state your method together with your number.

**Exercise 2.3 (Basin of the attractor).**
From 10 widely different initial conditions (including large ones), verify that after the transient **all** trajectories settle onto the same attractor (dissipation!). Estimate how long the transient lasts as a function of the distance of the initial condition from the attractor.


## Part 3 — Poincaré sections and routes to chaos

**Exercise 3.1 (Poincaré section).**
With `poincare = .true.`, plot the section at $y=0$. Instead of a point-cloud like in the pendulum, you should see the attractor's section as a set of **spiraling arcs** — quasi-1D curves. Explain why the section of the Chua attractor is so thin (hint: dissipation again: the attractor has fractal dimension barely above 2).

**Exercise 3.2 (Bifurcation diagram — the classic route to chaos).**
Fix $\beta = 100/7$, $a=-8/7$, $b=-5/7$, $\gamma=0$, and sweep $\alpha$ from about 6 to 10 (say 400 values). For each $\alpha$: integrate, **discard the transient**, record the local maxima of $x(\tau)$, and print pairs $(\alpha, x_{max})$ to a file (a small shell or Python driver that rewrites `input.dat` and runs `./chua` is expected — same workflow as the pendulum regularity map). Plot with `plot_chua.py bifurcation`. Identify:
1. the stable-equilibrium regime at small $\alpha$ (a single branch);
2. a **Hopf bifurcation** where the equilibrium becomes a limit cycle;
3. the **period-doubling cascade** (1 → 2 → 4 → 8 branches) ending in chaos;
4. periodic windows inside the chaotic band.
This figure, the "rosetta stone" of nonlinear dynamics, is fully computable in one afternoon.

**Exercise 3.3 (Hopf, analytically).**
Locate the Hopf bifurcation of Exercise 3.2 on the $x$ axis analytically: a Hopf occurs when a complex eigenvalue pair crosses the imaginary axis, i.e. when $c_1 c_2 = c_0$ (Routh–Hurwitz condition) for the appropriate fixed point. Solve for $\alpha$ and compare with the numerically observed onset in the bifurcation diagram within your parameter resolution.


## Part 4 — Going further (optional / project-level)

- **4.1 Full Lyapunov spectrum.** Implement the Benettin algorithm with three perturbation vectors and Gram–Schmidt re-orthonormalization to get $(\lambda_1, \lambda_2, \lambda_3)$. Two independent checks: (i) $\lambda_2 \approx 0$ (flow on the attractor); (ii) the **sum must equal the time-averaged trace** $\langle \nabla\cdot F\rangle = \langle -\alpha(1+f'(x)) - 1 - \gamma\rangle$ — average the trace along your trajectory and compare.
- **4.2 Lyapunov dimension.** From the spectrum, compute the Kaplan–Yorke dimension $D_{KY} = k + \sum_{i\le k}\lambda_i/|\lambda_{k+1}|$ and compare with the literature value $D \approx 2.13$ for the canonical double scroll (Matsumoto–Chua–Komuro 1985). Why is it "barely above 2", and how does that show in the Poincaré section?
- **4.3 From equations to hardware.** Look up the op-amp realization of the Chua diode and find the component values that realize $a \approx -8/7$, $b \approx -5/7$. Which real-world effects (saturation, noise) make a hardware attractor look different from the simulated one?
- **4.4 Literature search.** As for the double pendulum: find one research application built on the Chua oscillator (secure communication via chaos synchronization, random number generation, memristor-based Chua variants, neuromorphic computing) and write a half-page report with the standard four points (reference, mapping, simplification, role of chaos).
- **4.5 Chaos with and without dissipation.** Write one page comparing the Chua oscillator and the double pendulum as two complementary paradigms: attractor vs. energy shell, fractal dimension vs. KAM tori, bifurcations vs. transitions. Supporting figures: your Poincaré sections from both sections of the course side by side.


## Suggested fixes for `theory.md`

1. **Fixed-point table (§5) — a real physics error, not a typo:** the outer equilibria are $(x^{\pm},\ \tfrac{\gamma x^\pm}{\beta+\gamma},\ -\tfrac{\beta x^\pm}{\beta+\gamma})$. The table's entries $y=\gamma/(\beta+\gamma)$, $z=-\beta/(\beta+\gamma)$ hold **only when $x^\ast=\pm 1$**, which is why the x-column shows "1". For the canonical parameters the correct equilibria are $P_\pm = (\pm 1.5,\, 0,\, \mp 1.5)$; the values printed in the table ($(1.5,0,-1)$) leave a residual of $0.5$ in the equation $\dot y = x-y+z$.
2. **Normalization table (§4):** `b=R_Gb` should read $b = R\,G_b$; also note that the three namelist parameters $a, b$ are related to the diode slopes by $a = R\,G_a$ (inner) and $b = R\,G_b$ (outer).
3. **Fixed-point derivation (§5):** "Set $dx/d\tau=0, dy/d\tau=0, dx/d\tau=0$" — the third should be $dz/d\tau = 0$.
4. **Add the canonical parameter set explicitly** ($\alpha=9$, $\beta=100/7$, $a=-8/7$, $b=-5/7$, $\gamma \approx 0$): every numerical exercise refers to it as "the canonical double scroll".
5. **Add the closed-form characteristic coefficients** $c_2, c_1, c_0$ (as in Exercise 1.3) and the Routh–Hurwitz criterion $c_1 c_2 > c_0$ for stability — this grounds the Hopf-onset analysis of Exercise 3.3.
6. **Add a short Lyapunov-section:** definition of the exponents, the renormalized (Benettin) estimator, the trace-sum rule $\sum_i \lambda_i = \langle\nabla\cdot F\rangle$, and the Kaplan–Yorke dimension — the tools used in Exercises 2.2, 4.1, 4.2.
7. Minor: in §3 the second Kirchhoff rule application mentions "the loop containing $L$, $R_0$ and $C_2$" — a small sketch of which node/loop is which would help students reproduce the derivation unaided.
