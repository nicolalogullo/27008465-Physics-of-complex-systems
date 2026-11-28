# The Kicked Rotor: Quantum Chaos and Dynamical Localization

**References:**
- G. Casati, B. V. Chirikov, F. M. Izrailev, J. Ford, in *Stochastic Behavior in Classical and Quantum Systems* (Lecture Notes in Physics 93, Springer, 1979) — the discovery of dynamical localization.
- S. Fishman, D. R. Grempel, R. E. Prange, *Phys. Rev. A* 29, 1639 (1984) — mapping to the Anderson model.
- B. V. Chirikov, F. M. Izrailev, D. L. Shepelyansky, *Sov. Sci. Rev. C* 2, 209 (1981).
- F. L. Moore, J. C. Robinson, C. F. Bharucha, B. Sundaram, M. G. Raizen, *Phys. Rev. Lett.* 75, 4594 (1995) — realization with cold atoms.

---

## Part I — Theory

### 1. Motivation: what does "quantum chaos" mean?

The double pendulum and the Chua oscillator showed that classical systems can be unpredictable in principle: nearby trajectories separate exponentially. Quantum mechanics, however, poses an immediate problem:

- the Schrödinger equation is **linear**, so two initially close wavefunctions stay close;
- the evolution is **unitary**: norms (probabilities) are exactly conserved;
- Heisenberg's principle forbids the arbitrarily sharp initial conditions needed to trigger classical sensitivity.

"Quantum chaos" is therefore not exponential divergence of wavefunctions, but the study of **signatures**: how does quantization of a classically chaotic system differ from quantization of a regular one? The kicked rotor is the cleanest laboratory for this question — as central to quantum chaos as the double pendulum is to classical chaos:

- its classical limit is a textbook chaotic map;
- the quantum evolution is a **one-line exact algorithm** (one Floquet operator);
- it displays the paradigmatic phenomenon of **dynamical localization**: the quantum suppression of classical diffusion, the temporal analogue of Anderson localization in disordered solids;
- it has been realized experimentally with laser-cooled atoms (Moore et al. 1995).

---

### 2. The classical kicked rotor and the standard map

Consider a free rotor (angle $x \in [0,2\pi)$, momentum $p$) kicked periodically with period $T=1$ by a potential $K\cos x$:

$$
H(x,p,t) = \frac{p^2}{2} + K\cos(x)\sum_{n=-\infty}^{\infty}\delta(t-n).
$$

Integrating across each kick (the free part gives $x \to x + p$ between kicks, the delta kick gives $p \to p + K\sin x$) yields the **standard map** (Chirikov–Taylor map), one of the most studied objects in nonlinear dynamics:

$$
\begin{aligned}
p_{n+1} &= p_n + K\sin x_n,\\
x_{n+1} &= x_n + p_{n+1} \pmod{2\pi}.
\end{aligned}
$$

Properties to be verified numerically in the exercises:

1. **Area preservation.** The Jacobian determinant is exactly 1: the map is symplectic, i.e. the discrete analogue of a Hamiltonian flow.
2. **Order-to-chaos transition.** For $K \lesssim 1$ phase space is mostly regular (invariant circles block momentum transport — the KAM scenario of the pendulum section); for $K \gtrsim K_c \approx 1$ a connected chaotic sea appears and momentum diffuses without bound.
3. **Diffusion.** For $K \gg 1$ successive kicks are effectively uncorrelated, $\langle \sin x_n \sin x_m\rangle \approx \delta_{nm}/2$, so an ensemble started at rest obeys
$$\langle p^2\rangle_n = D\,n, \qquad D = \frac{K^2}{2}.$$
4. **Lyapunov exponent.** For $K\gg 1$, $\lambda \approx \ln(K/2)$, provided $\lambda > 0$ (e.g. $\lambda \approx 0.9$ for $K=5$).

---

### 3. Quantization: the Floquet operator

On the rotor, wavefunctions are periodic in $x$, so momentum is quantized:

$$\hat p |m\rangle = m\hbar\, |m\rangle, \qquad m \in \mathbb{Z},$$

where $\hbar$ here is the *effective, dimensionless* Planck constant set by the physical kick period and strengths (in the code it is simply a parameter). Because the Hamiltonian is periodic in time, the evolution over one period is a single **Floquet operator**:

$$
\hat U = e^{-i\hat p^2/(2\hbar)}\; e^{-iK\cos\hat x/\hbar},
\qquad \psi_{n+1} = \hat U \psi_n .
$$

Crucially, for a delta-kicked Hamiltonian this factorized form is **exact**, not a Trotter approximation: the kick acts instantaneously, the free rotation acts between kicks. This is the quantum analogue of the "exact discrete map" property of the classical limit.

Each factor is diagonal in one basis:

- $e^{-i\hat p^2/2\hbar}$ is diagonal in the **momentum basis** $|m\rangle$: it multiplies each amplitude by $e^{-i m^2\hbar/2}$;
- $e^{-iK\cos\hat x/\hbar}$ is diagonal in the **position basis**: it multiplies each $x_j = 2\pi j/N$ component by $e^{-iK\cos x_j/\hbar}$.

The two bases are connected by a discrete Fourier transform, which on a finite basis of $N$ momentum states is an $N\times N$ **unitary matrix** $F$. One Floquet step is therefore three matrix operations,

$$\psi \;\longmapsto\; e^{-iK\cos\hat x/\hbar}\, F\,\psi \;\longmapsto\; F^\dagger\!\left(\cdot\right) \;\longmapsto\; e^{-i\hat p^2/2\hbar}\left(\cdot\right),$$

with cost $O(N^2)$ per kick (an FFT would give $O(N\log N)$ — see the projects). No RK4, no step-size error, no external libraries: the only numerical parameter is the basis size $N$, which must be large enough that the localized state never touches its edges.

---

### 4. Dynamical localization

The central result. Run the classical map and the quantum evolution with the same $K$:

- **classical**: $\langle p^2\rangle = D\,n$ grows forever (diffusion);
- **quantum**: $\langle p^2\rangle$ follows the classical line only up to a **break time** $t^\ast$, then **saturates**: the probability distribution in momentum stops spreading and settles into an exponentially localized profile

$$|\psi_m|^2 \sim e^{-2|m|/\ell}, \qquad \ell \approx \frac{D}{\hbar^2}$$

(up to a fluctuating order-unity prefactor — see the exercises: single-exponential fits are window-dependent, and the localization length fluctuates from run to run).

Why does this happen? Interference between the many kick-to-kick paths to the same momentum state, which is invisible in the classical random-walk picture, progressively dephases the diffusion. Fishman, Grempel and Prange showed the mechanism precisely maps onto **Anderson localization** in a 1D disordered lattice: momentum states play the role of lattice sites, the quasi-energy disorder comes from the kick phases. Quantum mechanics does not "remove" chaos here — it removes *transport*.

Two quantitative consequences used as checks in the exercises:

- break time: $t^\ast \sim \dfrac{\langle p^2\rangle_\mathrm{sat}}{D} \sim \dfrac{D}{2\hbar^2}$ (order of magnitude), i.e. quantum mechanics follows classical mechanics for a time diverging as $\hbar\to0$ — this **is** the correspondence principle in a chaotic system;
- saturation: $\langle p^2\rangle_\mathrm{sat} \sim \dfrac{D^2}{2\hbar^2}$.

A final contrast worth meditating on, completing the course's arc: in the double pendulum (conservative, classical) chaos *spreads* trajectories over the energy shell; in Chua (dissipative, classical) chaos concentrates them onto a fractal attractor; in the kicked rotor (unitary, quantum) chaos is *muffled* into a localized, non-transporting state.

---

### 5. Numerical implementation and diagnostics

The reference Fortran code `kicked_rotor.f90` reads three namelists (`&rotor`, `&initial`, `&simulation`) and produces:

- `quantum.txt` — $\langle p^2\rangle$ and the norm per kick: the norm column is the built-in check (exact unitarity demands it stays 1 to machine precision, ~$10^{-15}$, unlike the RK4 energy drift of the pendulum);
- `momentum.txt` — the final $|\psi_m|^2$, for the localization-length fit;
- `classical.txt` — ensemble $\langle p^2\rangle$: the classical diffusion line $Dn$;
- `map.txt` — one classical trajectory, i.e. the standard-map phase portrait.

Defaults $K=5$, $\hbar=1$, $N=128$, 2000 kicks reproduce the theory figure: classical diffusion with $D=12.5$, quantum break at $t^\ast \approx 10$ kicks, saturation $\langle p^2\rangle_\mathrm{sat}\approx 100$–$140$ (fluctuating), exponentially localized momentum distribution with $\ell \approx 10$–$17$.
