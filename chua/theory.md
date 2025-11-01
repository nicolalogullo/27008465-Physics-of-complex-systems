# Lecture: Chua’s Oscillator — A Paradigm of Low-Dimensional Chaos
**References:**  
- Leon O. Chua, *IEEE Transactions on Circuits and Systems*, 1983–1986  
- Matsumoto T., *IEEE Trans. Circuits Syst.*, 1984  

---

## 🧮 Part I — Theoretical Lecture

### 1. Motivation and Context

- **Goal:** Understand how simple deterministic circuits can exhibit complex, aperiodic, and chaotic dynamics.  
- **Why Chua’s circuit?**  
  - It’s the *simplest physical chaotic oscillator*.  
  - Realizable with off-the-shelf electronic components.  
  - Ideal for exploring bifurcations, attractors, Lyapunov exponents, and fractal geometry.  

**Key chaotic oscillator families:**
Lorenz, Rössler, Duffing, **Chua** — each is a low-dimensional autonomous system.

---

### 2. Circuit Description

**Components**
| Symbol | Element | Role |
|---------|----------|------|
| C₁ | Capacitor | voltage \(v_1\) |
| C₂ | Capacitor | voltage \(v_2\) |
| L  | Inductor  | current \(i_L\) |
| R  | Resistor  | linear damping |
| N_R | Nonlinear resistor (“Chua diode”) | source of nonlinearity |

---

### 3. Governing Equations

Apply Kirchhoff’s laws:

$$
\begin{aligned}
C_1 \frac{dv_1}{dt} &= \frac{1}{R}(v_2 - v_1) - i_{NR}(v_1), \\
C_2 \frac{dv_2}{dt} &= \frac{1}{R}(v_1 - v_2) + i_L, \\
L \frac{di_L}{dt} &= -v_2.
\end{aligned}
$$

---

### 4. Nonlinear Resistor Characteristic

Piecewise-linear function:

$$ i_{NR}(v_1) = m_1 v_1 + \tfrac{1}{2}(m_0 - m_1)
  \big(|v_1 + B_p| - |v_1 - B_p|\big) $$

- $m_0$: inner (negative) slope  
- $m_1$: outer slope  
- $B_p$: break voltage  

This defines a **three-segment** $I–V$ curve.

---

### 5. Normalization and Dimensionless Form

Let
- $x = \frac{v_1}{B_p}$
- $y = \frac{v_2}{B_p}$
- $z = \frac{R i_L}{B_p}$
- $\tau = \frac{t}{R C_2}.$

Then
$$
\begin{aligned}
\dot{x} = \alpha (y - x - f(x)),
\dot{y} = x - y + z,\quad
\dot{z} = -\beta y,
\end{aligned}
$$

with
$$ f(x) = m_1 x + \tfrac{1}{2}(m_0 - m_1)(|x + 1| - |x - 1|), $$
and parameters
$$
\begin{aligned}
\alpha = \frac{C_2}{C_1},\\
\beta  = \frac{R^2 C_2}{L}.
\end{aligned}
$$

---

### 6. Equilibria and Linear Stability

Set $\dot{x}=\dot{y}=\dot{z}=0$:

$$
f(x^*) = -x^*, \quad y^* = 0, \quad z^* = -x^*.
$$

Hence there are **three fixed points**:
$$
E_0=(0,0,0), \quad E_{\pm}=(x_{\pm},0,-x_{\pm}).
$$

**Jacobian:**
$$
J =
\begin{pmatrix}
-\alpha(1+f'(x^*)) & \alpha & 0 \\
1 & -1 & 1 \\
0 & -\beta & 0
\end{pmatrix}
$$

Eigenvalues determine local type:
- Inner point → unstable focus or saddle-focus  
- Outer points → saddle or node  

---

### 7. Bifurcations and Routes to Chaos

Typical parameter variations cause:
- **Hopf bifurcation** → oscillations  
- **Period-doubling cascade** → chaos  
- **Crisis bifurcation** → sudden attractor enlargement  
- **Shilnikov homoclinic orbit** → guarantee of chaos for saddle-focus equilibrium

**Common parameter set:**
$$
\alpha = 9.0,\quad \beta = 14.286,\quad m_0=-1.143,\quad m_1=-0.714.
$$

Produces the famous **double-scroll attractor**.

---

### 8. Chaotic Dynamics and Observables

- **Double-scroll attractor:** trajectory alternates between two lobes.  
- **Sensitivity to initial conditions:** exponential divergence of nearby states.  
- **Lyapunov spectrum:** one positive, one zero, one negative exponent.  
- **Fractal dimension:** \(D_2 \approx 2.0\).

---

### 9. Quantitative Measures

**Largest Lyapunov exponent**
$$
\lambda_1 = \lim_{t \to \infty}
 \frac{1}{t} \ln \frac{||\delta X(t)||}{||\delta X(0)||}.
$$

**Kaplan–Yorke dimension**
$$
D_{KY} = j + \frac{\sum_{i=1}^{j}\lambda_i}{|\lambda_{j+1}|}.
$$

**Poincaré map**
- Plane $y=0,\ \dot{y}>0$
- Projection $(x,z)$ shows fractal cross-section of attractor.

---

### 10. Discussion

- Chua’s circuit demonstrates *deterministic chaos* with minimal components.  
- Connects nonlinear electronics, bifurcation theory, and complex systems.  
- Rich phenomena: multistability, intermittency, crises, chaos control.

---
