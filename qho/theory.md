# Quantum Chaos and Energy Level Statistics

## 1. Introduction

In classical mechanics, chaos is characterized by an extreme sensitivity to initial conditions.
In quantum mechanics, however, the notion of a trajectory does not exist and the Schrödinger
equation is linear. As a result, chaos cannot be defined in terms of trajectories or Lyapunov
exponents.

Instead, **quantum chaos** is identified through *statistical properties* of quantum systems,
in particular through the **statistics of their energy levels**.

---

## 2. Integrability and Chaos in Quantum Systems

A quantum system is said to be **integrable** if it possesses a complete set of commuting
conserved quantities. In such systems, the energy levels are largely uncorrelated.

A **quantum chaotic system** is one whose classical counterpart is chaotic and whose quantum
spectral properties follow universal predictions from **random matrix theory (RMT)**.

The key observable distinguishing these two regimes is the **distribution of energy level
spacings**.

---

## 3. Energy Levels and Level Spacings

Let
\[
E_1 \le E_2 \le E_3 \le \cdots
\]
be the ordered energy eigenvalues of a quantum Hamiltonian.

The **nearest-neighbor level spacing** is defined as:
\[
s_n = E_{n+1} - E_n
\]

To compare different systems, the spectrum is usually **unfolded** so that the mean spacing
is equal to one:
\[
\langle s \rangle = 1
\]

The statistical quantity of interest is the **probability distribution \(P(s)\)** of the
spacings \(s_n\).

---

## 4. Poisson Distribution: Integrable Systems

For **integrable quantum systems**, the energy levels are essentially uncorrelated.
As a consequence, level crossings are allowed and small spacings occur frequently.

The level spacing distribution is given by the **Poisson distribution**:
\[
P_{\text{Poisson}}(s) = e^{-s}
\]

### Key features:
- No level repulsion
- Finite probability at \(s = 0\)
- Typical of integrable systems
- Example: particle in a box, harmonic oscillator

---

## 5. Wigner Distribution: Quantum Chaotic Systems

For **quantum chaotic systems**, energy levels strongly repel each other.
Exact degeneracies are avoided unless enforced by symmetries.

The level spacing distribution follows the **Wigner–Dyson distribution**.
For systems with time-reversal symmetry, the relevant case is the **Gaussian Orthogonal Ensemble (GOE)**:
\[
P_{\text{Wigner}}(s) = \frac{\pi}{2} s \, e^{-\pi s^2 / 4}
\]

### Key features:
- Strong level repulsion
- \(P(s) \sim s\) for small \(s\)
- Universal behavior independent of microscopic details
- Signature of quantum chaos

---

## 6. Physical Interpretation: Level Repulsion

The difference between Poisson and Wigner statistics reflects the underlying dynamics:

- **Integrable systems** have independent quantum numbers → uncorrelated levels
- **Chaotic systems** have strongly mixed eigenstates → correlated levels

Level repulsion is a direct manifestation of the nonintegrability of the system.

---

## 7. Bohigas–Giannoni–Schmit Conjecture

The **Bohigas–Giannoni–Schmit (BGS) conjecture** states:

> Quantum systems whose classical counterparts are chaotic exhibit energy level statistics
> described by random matrix theory.

This conjecture has been verified numerically and experimentally for a wide class of systems.

---

## 8. Summary

- Quantum chaos cannot be defined via trajectories
- Energy level statistics provide a powerful diagnostic
- **Poisson distribution** → integrable systems
- **Wigner distribution** → quantum chaotic systems
- Level repulsion is the key signature of chaos

---

## 9. Outlook

Energy level statistics form the foundation of modern studies of:
- Floquet (periodically driven) systems
- Many-body quantum chaos
- Eigenstate thermalization
- Quantum information scrambling

They provide a universal bridge between classical chaos and quantum mechanics.

