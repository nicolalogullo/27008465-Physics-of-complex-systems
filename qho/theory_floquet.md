# Floquet Theory: A Simple Introduction

## 1. Why Floquet Theory?

Many physical systems are subject to **periodic driving**, such as:
- Atoms in oscillating electromagnetic fields
- Trapped ions driven by lasers
- Cold atoms in shaken optical lattices
- Driven anharmonic oscillators

In these systems, the Hamiltonian depends explicitly on time:
$H(t) = H(t + T)$

where $T$ is the period of the drive.

Because energy is **not conserved**, stationary energy eigenstates no longer exist.
Floquet theory provides the appropriate framework to describe such systems.

---

## 2. Analogy with Bloch’s Theorem

Floquet theory is the **time-domain analogue of Bloch’s theorem** for spatially
periodic systems.

- Bloch theorem: periodicity in space → crystal momentum
- Floquet theorem: periodicity in time → quasienergy

Just as momentum is defined modulo a reciprocal lattice vector,
**quasienergy is defined modulo $2\pi/T$**.

---

## 3. Floquet Theorem (Quantum Version)

Consider the time-dependent Schrödinger equation:
$i\hbar \frac{d}{dt} |\psi(t)\rangle = H(t) |\psi(t)\rangle$

If $H(t)$ is periodic in time, then solutions can be written as:
$|\psi_n(t)\rangle = e^{-i \varepsilon_n t / \hbar} |u_n(t)\rangle$

where:
- $|u_n(t)\rangle = |u_n(t+T)\rangle$ is periodic in time
- $\varepsilon_n$ is the **quasienergy**

The quasienergy plays a role analogous to energy in time-independent systems.

---

## 4. The Floquet Operator

Instead of diagonalizing the Hamiltonian, one studies the **time-evolution operator
over one period**:
$U(T) = \mathcal{T} \exp\!\left( -\frac{i}{\hbar} \int_0^T H(t)\, dt \right)$

The Floquet states are eigenstates of $U(T)$:
$U(T) |\phi_n\rangle = e^{-i \varepsilon_n T / \hbar} |\phi_n\rangle$

Key points:
- $U(T)$ is unitary
- Eigenvalues lie on the unit circle
- Quasienergies are extracted from the eigenphases

---

## 5. Quasienergies and Their Properties

### Periodicity
Quasienergies are defined modulo:
$\varepsilon \equiv \varepsilon + \frac{2\pi \hbar}{T}$

This means:
- The quasienergy spectrum is **compact**
- Only differences of quasienergies matter

### Physical interpretation
- Quasienergy determines phase accumulation per driving period
- Long-time dynamics are governed by Floquet eigenstates

---

## 6. Floquet Systems and Quantum Chaos

Floquet systems are especially important for quantum chaos because:

- Even a one-dimensional system becomes effectively **two-dimensional**
  (phase space + time)
- Periodic driving can destroy integrability
- Classical chaos emerges naturally

As a result, **Floquet quasienergy spectra** often show:
- Poisson statistics in regular regimes
- Wigner–Dyson statistics in chaotic regimes

---

## 7. Level Spacing Statistics in Floquet Systems

To study chaos in Floquet systems:

1. Compute quasienergies $\varepsilon_n$
2. Order them within one Floquet zone
3. Compute spacings:
   $s_n = \varepsilon_{n+1} - \varepsilon_n$
4. Normalize to unit mean spacing
5. Analyze the distribution $P(s)$

### Expected results:
- Regular dynamics → Poisson distribution
- Chaotic dynamics → Wigner distribution

---

## 8. Why Floquet Theory Is Ideal for Studying Chaos

Floque

