# Time Series Analysis — Comprehensive Theory Reference

**Course slides cross-referenced with:**
*Time Series Analysis*, Henrik Madsen, Chapman & Hall/CRC, 2008.

---

## Table of Contents

1. [Regression-Based Methods](#1-regression-based-methods)
2. [Exponential Smoothing and Local Trend Models](#2-exponential-smoothing-and-local-trend-models)
3. [Time Series Operators](#3-time-series-operators)
4. [Linear Dynamic Systems](#4-linear-dynamic-systems)
5. [Stochastic Processes](#5-stochastic-processes)
6. [ARMA Model Families](#6-arma-model-families)
7. [Identification: ACF and PACF](#7-identification-acf-and-pacf)
8. [Multivariate (Vector) Time Series](#8-multivariate-vector-time-series)
9. [Estimation of ARMA/ARMAX Models](#9-estimation-of-armaarmax-models)
10. [State Space Models and the Kalman Filter](#10-state-space-models-and-the-kalman-filter)
11. [Named Theorems, Definitions, and Rules — Master Index](#11-named-theorems-definitions-and-rules--master-index)

---

## 1. Regression-Based Methods

**Textbook:** Chapter 3 (p. 31–68), especially Sections 3.1–3.4.

### 1.1 The General Linear Model (GLM)

The model is:

```
Y = x θ + ε
```

where `Y` is the (N×1) response vector, `x` is the (N×p) design matrix, `θ` is the (p×1) parameter vector, and `ε ~ N(0, σ² I)`.

**Textbook:** Section 3.2, p. 33.

### 1.2 Ordinary Least Squares (OLS)

**Definition:** The OLS estimator minimises the sum of squared residuals:

```
S(θ) = (Y - x θ)ᵀ(Y - x θ)
```

**OLS Estimator:**

```
θ̂ = (xᵀx)⁻¹ xᵀ Y  =  F⁻¹ h
```

where `F = xᵀx` and `h = xᵀY`.

**Textbook:** Section 3.2.1, p. 34.

**Properties:** BLUE (Best Linear Unbiased Estimator) under Gauss-Markov conditions.

### 1.3 Weighted Least Squares (WLS)

**When to use:** When observations have unequal variances, i.e. `Var(ε) = σ² Σ` with `Σ ≠ I`.

**WLS Criterion:**

```
S(θ) = (Y - x θ)ᵀ Σ⁻¹ (Y - x θ)
```

**WLS Estimator:**

```
θ̂ = (xᵀ Σ⁻¹ x)⁻¹ xᵀ Σ⁻¹ Y
```

or equivalently:

```
θ̂_N = F_N⁻¹ h_N
F_N = Σ_{j=0}^{N-1} λʲ f(-j) fᵀ(-j)
h_N = Σ_{j=0}^{N-1} λʲ f(-j) Y_{N-j}
```

In the local (forgetting factor) formulation: `Σ = diag[1/λ^{N-1}, ..., 1/λ, 1]`, giving exponentially decreasing weights.

**Textbook:** Section 3.2.1 (p. 34), local formulation in Section 3.4 (p. 47).

### 1.4 Maximum Likelihood (ML) Estimation

Under normality `ε ~ N(0, σ² I)`, the ML estimator coincides with OLS. The ML estimator for `σ²` is biased:

```
σ̂²_ML = εᵀε / N
```

while the unbiased estimator is:

```
σ̂² = εᵀε / (N - p)
```

**Textbook:** Section 3.2.2, p. 40.

### 1.5 Prediction in the General Linear Model

**l-step-ahead prediction:**

```
Ŷ_{N+l|N} = fᵀ(l) θ̂_N
```

**Variance of the prediction error:**

```
V[Y_{N+l} - Ŷ_{N+l|N}] = σ² [1 + fᵀ(l) F_N⁻¹ f(l)]
```

**100(1-α)% prediction interval:**

```
Ŷ_{N+l|N} ± t_{α/2}(N-p) σ̂ √(1 + fᵀ(l) F_N⁻¹ f(l))
```

where `σ̂² = εᵀε/(N-p)`.

**Textbook:** Section 3.3 (p. 44), Section 3.3.1 (p. 45).

### 1.6 Trend Models

**General trend model:**

```
Y_{N+j} = fᵀ(j) θ + ε_{N+j}
```

The 2×2 transition matrix **L** satisfies `f(j+1) = L f(j)`.

**Standard trend models** (Textbook Section 3.4, p. 47):

| Model | Formula | f(j) |
|---|---|---|
| Constant mean | `Y_{N+j} = θ₀ + ε_{N+j}` | `(1)` |
| Linear trend | `Y_{N+j} = θ₀ + θ₁j + ε_{N+j}` | `(1, j)ᵀ` |
| Quadratic trend | `Y_{N+j} = θ₀ + θ₁j + θ₂ j²/2 + ε_{N+j}` | `(1, j, j²/2)ᵀ` |
| k-th order polynomial | `Y_{N+j} = θ₀ + θ₁j + ... + θₖ jᵏ/k! + ε_{N+j}` | polynomial basis |
| Harmonic (period p) | `Y_{N+j} = θ₀ + θ₁ sin(2π/p · j) + θ₂ cos(2π/p · j) + ε_{N+j}` | trigonometric |

For linear trend: `L = [[1,0],[1,1]]`.

**Textbook:** Section 3.4, p. 47–58.

### 1.7 Iterative (Recursive) Parameter Updating — Global Model

When a new observation `Y_{N+1}` arrives, update without recomputing from scratch:

```
F_{N+1} = F_N + f(-N) fᵀ(-N)
h_{N+1} = L⁻¹ h_N + f(0) Y_{N+1}
θ̂_{N+1} = F_{N+1}⁻¹ h_{N+1}
```

**Textbook:** Section 3.2.1 (recursive form), Section 3.6 (example, p. 62).

---

## 2. Exponential Smoothing and Local Trend Models

**Textbook:** Section 3.4 (p. 47–58), Section 3.4.2 (p. 50), Section 3.4.4 (p. 56).

### 2.1 Simple Exponential Smoothing (SES)

**Definition (Simple Exponential Smoothing):**
Given forgetting factor `λ ∈ ]0,1[`, define the sequence `S_N` by:

```
S_N = (1 - λ) Y_N + λ S_{N-1}
```

where `α = 1 - λ` is the **smoothing constant**. The weighted mean estimator is:

```
μ̂_N = c Σ_{j=0}^{N-1} λʲ Y_{N-j}
```

with `c = (1-λ)/(1-λᴺ)`. For large N: `c ≈ 1-λ` so:

```
μ̂_N ≈ (1-λ) Y_N + λ μ̂_{N-1}
```

**Used as a prediction model:** `Ŷ_{N+l|N} = μ̂_N` for all horizons l (constant forecast).

**Updating predictions with new observations:**

```
Ŷ_{N+l+1|N+1} = (1-λ) Y_{N+1} + λ Ŷ_{N+l|N}
```

**Textbook:** Section 3.4.2, p. 50.

### 2.2 Choosing the Smoothing Constant α = 1 - λ

Minimise the sum of squared one-step-ahead prediction errors:

```
S(α) = Σ_{t=1}^{N} (Y_t - Ŷ_{t|t-1}(α))²
```

The optimal `α` depends on the forecast horizon. Longer horizons → smaller `α` (more smoothing).

**Textbook:** Section 3.4.2 (p. 50).

### 2.3 Local Trend Models (WLS with Forgetting Factor)

Forget old observations exponentially:

```
θ̂_N = argmin_θ S(θ; N)
```

where:

```
S(θ; N) = Σ_{j=0}^{N-1} λʲ [Y_{N-j} - fᵀ(-j) θ]²
```

This is a WLS criterion with `Σ = diag[1/λ^{N-1}, ..., 1/λ, 1]`.

**WLS solution:**

```
θ̂_N = F_N⁻¹ h_N
F_N = Σ_{j=0}^{N-1} λʲ f(-j) fᵀ(-j)
h_N = Σ_{j=0}^{N-1} λʲ f(-j) Y_{N-j}
```

**Textbook:** Section 3.4.4, p. 56.

### 2.4 Recursive Updating for Local Model

When `Y_{N+1}` arrives:

```
F_{N+1} = F_N + λᴺ f(-N) fᵀ(-N)
h_{N+1} = λ L⁻¹ h_N + f(0) Y_{N+1}
θ̂_{N+1} = F_{N+1}⁻¹ h_{N+1}
```

Initial values: `h_0 = 0`, `F_0 = 0`.

**Stationary limit** (for many standard functions `f`): as `N → ∞`, `λᴺ f(-N) fᵀ(-N) → 0`, so `F_{N+1} → F_N = F`. The update simplifies to:

```
θ̂_{N+1} = Lᵀ θ̂_N + F⁻¹ f(0) [Y_{N+1} - Ŷ_{N+1|N}]
```

**Textbook:** Section 3.4.4, p. 56.

### 2.5 Total Memory and Variance Estimation

**Total memory:**

```
T = Σ_{j=0}^{N-1} λʲ = (1 - λᴺ)/(1 - λ)
```

T measures how many observations the estimation is essentially based upon.

**Variance estimator in local trend models:**

```
σ̂² = (Y - x_N θ̂_N)ᵀ Σ⁻¹ (Y - x_N θ̂_N) / (T - p),    T > p
```

Note: This estimator is not in the textbook.

### 2.6 Holt-Winters Procedure (Seasonal Exponential Smoothing)

For time series with seasonal variation, the Holt-Winters procedure extends exponential smoothing to handle both trend and seasonal components.

**Textbook:** Section 3.5.2, p. 61.

---

## 3. Time Series Operators

**Textbook:** Section 4.5, p. 87–89.

### 3.1 Backshift Operator B

```
B Xₜ = X_{t-1}
Bᵏ Xₜ = X_{t-k}
```

### 3.2 Forward Shift Operator F

```
F Xₜ = X_{t+1}
Fᵏ Xₜ = X_{t+k}
```

Note: `F = B⁻¹`.

### 3.3 Difference Operator ∇

```
∇ Xₜ = Xₜ - X_{t-1} = (1 - B) Xₜ
∇ᵈ Xₜ = (1 - B)ᵈ Xₜ   (d-th order differencing)
```

### 3.4 Summation Operator S

```
S Xₜ = Σ_{j=0}^{∞} Xₜ₋ⱼ = (1/(1-B)) Xₜ
```

S is the inverse of ∇: `S = ∇⁻¹`.

### 3.5 Seasonal Operators

```
∇_s Xₜ = (1 - Bˢ) Xₜ = Xₜ - X_{t-s}     (seasonal differencing)
```

**Textbook:** Section 4.5, p. 87.

---

## 4. Linear Dynamic Systems

**Textbook:** Chapter 4 (p. 69–96).

### 4.1 Linear Systems in the Time Domain

A linear time-invariant (LTI) system is characterised by its **impulse response** `{h_k}`:

```
Y_t = Σ_{k=0}^{∞} h_k u_{t-k}  =  H(B) u_t
```

where `H(B) = Σ_{k=0}^{∞} h_k Bᵏ` is the transfer function operator.

**Textbook:** Section 4.1, p. 70.

### 4.2 Step Response

The **step response** `{g_k}` is the output when the input is a unit step (u_t = 1 for t ≥ 0):

```
g_k = Σ_{j=0}^{k} h_j
```

The impulse response is recovered as: `h_k = g_k - g_{k-1}`.

**Textbook:** Section 4.1, p. 70.

### 4.3 Stability

A LTI system is **stable** (BIBO stable) if the impulse response is absolutely summable:

```
Σ_{k=0}^{∞} |h_k| < ∞
```

Equivalently, for rational systems `H(z) = B(z)/A(z)`: all poles must lie inside the unit circle.

**Textbook:** Section 4.1, p. 70.

### 4.4 The z-Transform

The z-transform of a sequence `{x_k}`:

```
X(z) = Σ_{k=-∞}^{∞} x_k z⁻ᵏ
```

Properties:
- Convolution in time domain → multiplication in z-domain
- Shift: `Z{x_{k-d}} = z⁻ᵈ X(z)`
- The backshift operator `B` corresponds to `z⁻¹`

**Textbook:** Section 4.4, p. 80.

### 4.5 Frequency Domain (Linear Systems)

The **frequency response** is `H(e^{iω})`, the z-transform evaluated on the unit circle.

The **power spectral density** of the output `Y_t = H(B) X_t` is:

```
S_Y(ω) = |H(e^{iω})|² S_X(ω)
```

**Textbook:** Section 4.2, p. 73.

### 4.6 Convolution

For two sequences `{a_k}` and `{b_k}`, their convolution is:

```
(a * b)_k = Σ_j a_j b_{k-j}
```

In the z-domain: `(A * B)(z) = A(z) · B(z)`.

**Textbook:** Section 4.1, p. 70.

---

## 5. Stochastic Processes

**Textbook:** Chapter 5 (p. 97–143).

### 5.1 Definition of a Stochastic Process

A **stochastic process** `{X_t, t ∈ T}` is a collection of random variables indexed by time. A specific observed sequence is a **realisation** of the process.

**Characterisation:** Fully described by all finite-dimensional distributions `F(x_{t₁}, ..., x_{tₙ})` for all n and all time indices.

**Textbook:** Section 5.1, p. 97; Section 5.2, p. 97.

### 5.2 Second-Order Moment Representation

For practical purposes, a stochastic process is often characterised by its second-order moments:
- **Mean function:** `μ(t) = E[X_t]`
- **Variance function:** `σ²(t) = Var(X_t) = E[(X_t - μ(t))²]`
- **Autocovariance function:** `γ(t₁, t₂) = Cov(X_{t₁}, X_{t₂}) = E[(X_{t₁} - μ(t₁))(X_{t₂} - μ(t₂))]`
- **Autocorrelation function:** `ρ(t₁, t₂) = γ(t₁, t₂) / √(σ²(t₁) σ²(t₂))`

**Textbook:** Section 5.2, p. 97; Section 5.2.2, p. 103.

### 5.3 Stationarity

#### Strong (Strict) Stationarity

A process is **strongly stationary** if all finite-dimensional distributions are invariant to time shifts:

```
(X_{t₁}, ..., X_{tₙ}) =ᵈ (X_{t₁+τ}, ..., X_{tₙ+τ})  for all τ, n, t₁,...,tₙ
```

**Textbook:** Section 5.2.1.1, p. 99.

#### Weak (Second-Order) Stationarity

A process is **weakly stationary** (or covariance-stationary) if:
1. `μ(t) = μ` (constant mean)
2. `σ²(t) = σ²` (constant variance)
3. `γ(t₁, t₂) = γ(t₁ - t₂)` (autocovariance depends only on the lag `k = t₁ - t₂`)

For a weakly stationary process, write `γ(k) = Cov(X_t, X_{t+k})` and `ρ(k) = γ(k)/γ(0)`.

**Textbook:** Section 5.2.1.1, p. 99.

**Properties of ACF for stationary processes:**
- `γ(0) = σ² ≥ 0`
- `γ(-k) = γ(k)` (symmetry)
- `|γ(k)| ≤ γ(0)` (Cauchy-Schwarz)
- The autocovariance function is **positive semi-definite**

**Textbook:** Section 5.2.2, p. 103.

### 5.4 Ergodicity

A stationary process is **mean-ergodic** if the time average converges to the ensemble mean:

```
(1/N) Σ_{t=1}^{N} X_t  →  E[X_t] = μ   as N → ∞
```

**Mean-ergodic condition:**

```
Σ_{k=-∞}^{∞} |γ(k)| < ∞
```

Ergodicity allows estimation of population moments from a single realisation.

**Textbook:** Section 5.2.1, p. 99.

### 5.5 Classes of Stochastic Processes

**Textbook:** Section 5.2, p. 97; Section 5.3, p. 107.

#### Gaussian (Normal) Process

A process `{X_t}` is **Gaussian** if all finite-dimensional distributions are multivariate normal. For Gaussian processes, weak stationarity implies strong stationarity.

#### Markov Process

A process is a **Markov process** if:

```
P(X_{t+1} ≤ x | X_t, X_{t-1}, ...) = P(X_{t+1} ≤ x | X_t)
```

(the future depends on the past only through the present state).

#### Decomposition Theorem

Any stochastic process can be decomposed as:

```
X_t = S_t + D_t
```

where `S_t` is a purely stochastic component and `D_t` is a deterministic component.

**Textbook:** Section 5.3, p. 107.

#### Pure Stochastic Process (White Noise)

A sequence `{ε_t}` is **white noise** if:
- `E[ε_t] = 0`
- `Var(ε_t) = σ²`
- `Cov(ε_t, ε_s) = 0` for `t ≠ s`

Written: `ε_t ~ WN(0, σ²)`. If also normally distributed: `ε_t ~ N(0, σ²)`.

### 5.6 Correlation and R²

The **coefficient of determination R²** measures the proportion of variance explained by the model:

```
R² = 1 - SS_res / SS_tot = 1 - Σ(Y_t - Ŷ_t)² / Σ(Y_t - Ȳ)²
```

**Textbook:** Section 3.2 (regression context), p. 33.

---

## 6. ARMA Model Families

**Textbook:** Chapter 5, Sections 5.5–5.6 (p. 117–134).

### 6.1 Moving Average Process MA(q)

**Definition:**

```
X_t = ε_t + θ₁ε_{t-1} + ... + θ_q ε_{t-q} = θ(B) ε_t
```

where `θ(B) = 1 + θ₁B + ... + θ_q Bᵠ` and `ε_t ~ WN(0, σ²)`.

**Properties:**
- Always weakly stationary (for any θ coefficients)
- Mean: `E[X_t] = 0`
- Variance: `γ(0) = σ²(1 + θ₁² + ... + θ_q²)`
- ACF cuts off after lag q: `ρ(k) = 0` for `|k| > q`
- PACF decays exponentially (or oscillates)

**Autocovariance (eq. 5.65 in textbook):**

```
γ(k) = σ² Σ_{j=0}^{q-k} θⱼ θ_{j+k},   k = 0, 1, ..., q
γ(k) = 0,                                 k > q
```

**Textbook:** Section 5.5.1, p. 117; autocovariance eq. (5.65).

### 6.2 Invertibility of MA(q)

An MA(q) process is **invertible** if the roots of `θ(z⁻¹) = 0` all lie inside the unit circle (equivalently, roots of `θ(z) = 0` lie outside the unit circle).

An invertible MA can be written as an infinite AR:

```
π(B) X_t = ε_t    where    π(B) = θ⁻¹(B)
```

**When to check:** Invertibility ensures a unique representation and is needed for the Wold decomposition.

**Textbook:** Section 5.5.1, p. 117.

### 6.3 Autoregressive Process AR(p)

**Definition:**

```
X_t = φ₁ X_{t-1} + ... + φ_p X_{t-p} + ε_t
φ(B) X_t = ε_t
```

where `φ(B) = 1 - φ₁B - ... - φ_p Bᵖ`.

**Properties:**
- Stationary if roots of `φ(z⁻¹) = 0` are inside the unit circle
- ACF decays exponentially (or with damped oscillations)
- PACF cuts off after lag p: `φ_{kk} = 0` for `k > p`

**AR(1) example:**
- `X_t = φ X_{t-1} + ε_t`
- Stationary iff `|φ| < 1`
- `γ(k) = φᵏ σ²/(1-φ²)`, `ρ(k) = φᵏ`

**Yule-Walker Equations for AR(p):**

```
γ(k) = φ₁ γ(k-1) + φ₂ γ(k-2) + ... + φ_p γ(k-p),    k ≥ 1
```

In matrix form (Theorem 5.10):

```
Γ φ = γ
```

where `Γ` is the (p×p) autocovariance matrix, `φ = (φ₁,...,φ_p)ᵀ`, `γ = (γ(1),...,γ(p))ᵀ`.

**Textbook:** Section 5.5.2, p. 119; Theorem 5.10 (Yule-Walker).

### 6.4 Stationarity Condition for AR(p)

The AR(p) process is **stationary** if and only if all roots of:

```
φ(z⁻¹) = 0
```

lie strictly inside the unit circle (i.e., `|z| < 1`).

Equivalently, all roots of the characteristic polynomial `1 - φ₁z - ... - φ_p zᵖ = 0` lie outside the unit circle.

**Textbook:** Section 5.5.2, p. 119.

### 6.5 ARMA(p,q) Process

**Definition:**

```
φ(B) X_t = θ(B) ε_t
(1 - φ₁B - ... - φ_p Bᵖ) X_t = (1 + θ₁B + ... + θ_q Bᵠ) ε_t
```

**Stationarity condition:** roots of `φ(z⁻¹) = 0` inside the unit circle.

**Invertibility condition:** roots of `θ(z⁻¹) = 0` inside the unit circle.

**ACF/PACF behaviour:**
- ACF: decays exponentially after lag (q-p) if q ≥ p, or purely exponentially
- PACF: decays exponentially after lag (p-q) if p ≥ q

**Autocovariance of ARMA(p,q) (eqs. 5.100, 5.101 in textbook):**

For `k > q`:
```
γ(k) = φ₁ γ(k-1) + ... + φ_p γ(k-p)       (homogeneous recursion)
```

For `k = 0, 1, ..., q` (initial conditions involve MA terms, Textbook eq. 5.101).

**Textbook:** Section 5.5.3, p. 125; autocovariance eqs. (5.100) and (5.101).

### 6.6 ARIMA(p,d,q) — Non-stationary Models

**Definition:**

```
φ(B) ∇ᵈ X_t = θ(B) ε_t
```

where `∇ = 1 - B` and d is the degree of differencing. The differenced series `W_t = ∇ᵈ X_t` follows a stationary ARMA(p,q) model.

**When to use:** When the original series is non-stationary (integrated of order d). Identified by ACF that decays very slowly.

**Textbook:** Section 5.6.1, p. 130.

### 6.7 Seasonal ARIMA (SARIMA)

**Definition:** `ARIMA(p,d,q)×(P,D,Q)_s`

```
Φ(Bˢ) φ(B) ∇ᴰ_s ∇ᵈ X_t = Θ(Bˢ) θ(B) ε_t
```

**Textbook:** Section 5.6.2, p. 132.

### 6.8 ARX and ARMAX Models (with Exogenous Inputs)

**ARX (AutoRegressive with eXogenous input):**

```
φ(B) Y_t = ω(B) u_t + ε_t
```

where `u_t` is an observed exogenous input.

**ARMAX (ARX + Moving Average):**

```
φ(B) Y_t = ω(B) u_t + θ(B) ε_t
```

**Estimation of ARX:** Ordinary Least Squares (OLS) directly applicable since the model is linear in parameters.

**Textbook:** Section 5.6.3 (p. 134), Section 6.4.2 (p. 159).

---

## 7. Identification: ACF and PACF

**Textbook:** Chapter 6, Sections 6.2–6.3 (p. 146–156).

### 7.1 Sample Autocovariance Function

```
ĉ(k) = (1/N) Σ_{t=k+1}^{N} (Y_t - Ȳ)(Y_{t-k} - Ȳ)
```

### 7.2 Sample Autocorrelation Function (ACF)

```
r(k) = ĉ(k) / ĉ(0)
```

**95% confidence bands for white noise:** `±1.96/√N`

**Textbook:** Section 6.2.1, p. 146.

### 7.3 Partial Autocorrelation Function (PACF)

The **partial autocorrelation at lag k** is the correlation between `X_t` and `X_{t-k}` after removing the linear dependence on `X_{t-1}, ..., X_{t-k+1}`.

Computed recursively using the **Durbin-Levinson algorithm** (or Yule-Walker equations):

```
φ_{kk} = [ρ(k) - Σ_{j=1}^{k-1} φ_{k-1,j} ρ(k-j)] / [1 - Σ_{j=1}^{k-1} φ_{k-1,j} ρ(j)]
```

**Textbook:** Section 6.5.1, p. 171; Appendix B (Partial autocorrelations), p. 357.

### 7.4 ACF/PACF Patterns for Model Identification

| Model | ACF | PACF |
|---|---|---|
| MA(q) | Cuts off after lag q | Decays exponentially |
| AR(p) | Decays exponentially | Cuts off after lag p |
| ARMA(p,q) | Decays exponentially after lag q-p | Decays exponentially after lag p-q |
| ARIMA(p,d,q) | Very slow decay (non-stationary) | — |

**Textbook:** Section 6.3.2, p. 154.

### 7.5 Information Criteria for Model Order Selection

**AIC (Akaike Information Criterion):**

```
AIC = -2 ln L̂ + 2k
```

**BIC / SBC (Bayesian Information Criterion):**

```
BIC = -2 ln L̂ + k ln N
```

where `k` is the number of estimated parameters and `L̂` is the maximised likelihood.

**Textbook:** Section 6.5.3, p. 174.

---

## 8. Multivariate (Vector) Time Series

**Textbook:** Chapter 9 (p. 247–281).

### 8.1 Vector ARMA (VARMA) Process

A **p-dimensional VARMA(P,Q)** process:

```
Φ(B) Yₜ = Θ(B) εₜ
```

where `Yₜ = (Y₁ₜ, ..., Y_pₜ)ᵀ`, `εₜ ~ WN(0, Σ_ε)` (p-dimensional white noise), and:

```
Φ(B) = I - Φ₁B - ... - Φ_P Bᴾ    (p×p matrix polynomials)
Θ(B) = I + Θ₁B + ... + Θ_Q Bᴼ
```

**Stationarity condition:** all roots of `det(Φ(z⁻¹)) = 0` lie strictly inside the unit circle.

**Invertibility condition:** all roots of `det(Θ(z⁻¹)) = 0` lie strictly inside the unit circle.

**Textbook:** Section 9.3, p. 254.

### 8.2 Vector AR (VAR) Process

**VAR(p):**

```
Yₜ = Φ₁ Y_{t-1} + ... + Φ_P Y_{t-P} + εₜ
```

**Bivariate VAR(1) example** (NO/NO₂ air quality):

```
[Y₁ₜ]   [φ₁₁  φ₁₂] [Y₁,t₋₁]   [ε₁ₜ]
[Y₂ₜ] = [φ₂₁  φ₂₂] [Y₂,t₋₁] + [ε₂ₜ]
```

The cross-terms `φ₁₂`, `φ₂₁` capture Granger causality between the series.

**Textbook:** Section 9.3.4, p. 260.

### 8.3 Auto Covariance Matrix Function (ACMF)

For a stationary multivariate process `{Yₜ}` with mean `μ_Y`:

```
Γ_k = E[(Y_{t-k} - μ_Y)(Yₜ - μ_Y)ᵀ]
```

**Properties:**
- `Γ_k = Γ_{-k}ᵀ` (not symmetric in general, unlike the univariate case)
- `Γ₀ = E[(Yₜ - μ_Y)(Yₜ - μ_Y)ᵀ]` is the variance-covariance matrix

**Theoretical autocovariance for VARMA** is derived via:
- **Theorem 5.10** (Yule-Walker in matrix form) — for the AR part
- **Eq. (5.65)** — pure MA autocovariance formula
- **Eqs. (5.100) and (5.101)** — ARMA autocovariance recursion

These univariate results extend to the multivariate case in Chapter 9.

**Textbook:** Section 9.1, p. 249; Section 9.3.1, p. 255.

### 8.4 Sample Correlation Matrix Function (SCMF)

```
R_k = D⁻¹ Γ̂_k D⁻¹
```

where `D = diag(σ̂₁, ..., σ̂_p)` is the diagonal matrix of sample standard deviations.

**Textbook:** Section 9.6, p. 267; Section 9.3.2, p. 259.

### 8.5 Sample Partial Correlation Matrix Function (SPCMF)

Denoted `S_k`. Used in identification of multivariate models analogously to how PACF is used for univariate AR order identification.

**Textbook:** Section 9.3.2, p. 259; Section 9.3.3 (q-conditioned partial correlation), p. 260.

### 8.6 Multivariate Yule-Walker Equations

For a VAR(P) model, the Yule-Walker equations in matrix form:

```
[Γ₀    Γ₁   ... Γ_{P-1}] [Φ₁ᵀ]   [Γ₁ᵀ]
[Γ₁ᵀ   Γ₀   ... Γ_{P-2}] [Φ₂ᵀ] = [Γ₂ᵀ]
[...   ...  ...  ...   ] [...] = [...]
[Γ_{P-1}ᵀ ... Γ₀       ] [Φ_Pᵀ]   [Γ_Pᵀ]
```

These are solved to obtain the AR coefficient matrices `Φ₁, ..., Φ_P`.

**Textbook:** Section 9.3 (theoretical covariance matrix functions), p. 254.

### 8.7 Multivariate Prewhitening

For identifying the impulse response of a transfer function model (relating input `u_t` to output `Y_t`):

1. Fit a univariate ARIMA model to the input `u_t`: `π(B) u_t = α_t`
2. Apply the same filter to the output: `π(B) Y_t = β_t`
3. Cross-correlate the prewhitened series `{α_t}` and `{β_t}` to identify the impulse response structure

This removes autocorrelation in the input that would otherwise confound identification.

**Textbook:** Section 9.6.1, p. 269.

---

## 9. Estimation of ARMA/ARMAX Models

**Textbook:** Chapter 6 (p. 145–185), Chapter 9 (p. 247–281).

### 9.1 Moment Estimation (Method of Moments)

Equate sample autocovariances to theoretical autocovariances. For AR(p): solve the Yule-Walker equations using sample ACF values.

**Textbook:** Section 6.4.1, p. 157; Section 9.7.1, p. 270.

### 9.2 Least Squares (LS) for ARX

For the ARX model `φ(B) Y_t = ω(B) u_t + ε_t`, write in the linear regression form:

```
Y_t = φ₁ Y_{t-1} + ... + φ_p Y_{t-p} + ω₀ u_t + ... + ω_r u_{t-r} + ε_t
```

OLS directly applies. LS estimator is consistent under the assumption that `ε_t` is white noise.

**Textbook:** Section 6.4.2, p. 159.

### 9.3 Prediction Error Method

Minimise the sum of squared one-step-ahead prediction errors:

```
S(θ) = Σ_t ε_t²(θ)
```

where `ε_t(θ)` are the innovations computed recursively from the model. This is the standard approach for ARMA estimation.

**Textbook:** Section 6.4.3, p. 163.

### 9.4 Maximum Likelihood (ML) for ARMA

Under Gaussian innovations, ML estimation of ARMA parameters maximises:

```
L(θ, σ²) = -(N/2) ln(2πσ²) - (1/2σ²) Σ_t ε_t²(θ)
```

For fixed `σ²`, this reduces to the prediction error method.

**Textbook:** Section 6.4.4, p. 166.

### 9.5 Spliid Method (1983) for Multivariate ARMAX

The **Spliid method** (Spliid, 1983) is an extended LS algorithm for multivariate ARMAX models. It provides consistent and asymptotically efficient estimates without requiring full ML optimisation.

**Procedure:**
1. Fit preliminary AR models to each output
2. Use the residuals as proxies for the MA innovations
3. Apply LS to the augmented regression

**Textbook:** Section 9.7.2, p. 271 ("An extended LS method for multivariate ARMAX models (the Spliid method)").

### 9.6 MARIMA (Multivariate ARIMA)

Software implementation for fitting multivariate ARIMA (VARIMA) models, used in conjunction with the Spliid method and ML estimation.

**Textbook:** Section 9.7, p. 269.

### 9.7 Model Checking

After fitting, verify:
1. **Residual analysis:** residuals should be white noise (check ACF/PACF of residuals)
2. **Cross-validation:** evaluate forecasting performance on held-out data
3. **Information criteria:** use AIC/BIC for model comparison

**Textbook:** Section 6.6, p. 174.

---

## 10. State Space Models and the Kalman Filter

**Textbook:** Chapter 10 (p. 283–311).

### 10.1 Linear Stochastic State Space Model

The linear stochastic state space model consists of two equations:

**System (state transition) equation:**

```
X_t = A X_{t-1} + B u_{t-1} + e_{1,t}
```

**Observation equation:**

```
Y_t = C X_t + e_{2,t}
```

where:
- `X_t` is the (n×1) state vector (latent/unobserved)
- `Y_t` is the (m×1) observation vector
- `u_t` is the (r×1) input vector (exogenous, known)
- `A` (n×n): state transition matrix
- `B` (n×r): input matrix
- `C` (m×n): observation matrix
- `e_{1,t} ~ WN(0, Q)`: system noise (process noise)
- `e_{2,t} ~ WN(0, R)`: observation noise (measurement noise)
- `e_{1,t}` and `e_{2,t}` are uncorrelated

**Textbook:** Section 10.1, p. 284.

### 10.2 Relationship to ARMA

Any stationary ARMA(p,q) model can be written in state space form. This provides a connection between the ARMA and state space representations.

**Textbook:** Section 10.2, p. 286 (Transfer function and state space formulations).

### 10.3 The Kalman Filter

The **Kalman filter** provides optimal (minimum variance) linear estimates of the state `X_t` given observations up to time t.

**Kalman Filter Recursion:**

**Prediction step** (time update):
```
X̂_{t|t-1} = A X̂_{t-1|t-1} + B u_{t-1}
P_{t|t-1} = A P_{t-1|t-1} Aᵀ + Q
```

**Update step** (measurement update):
```
K_t = P_{t|t-1} Cᵀ (C P_{t|t-1} Cᵀ + R)⁻¹        (Kalman gain)
X̂_{t|t} = X̂_{t|t-1} + K_t (Y_t - C X̂_{t|t-1})    (filtered estimate)
P_{t|t} = (I - K_t C) P_{t|t-1}                      (error covariance)
```

The innovation `ν_t = Y_t - C X̂_{t|t-1}` is the one-step-ahead prediction error.

**Textbook:** Section 10.3.1, p. 289.

### 10.4 k-step Predictions in State Space

The k-step-ahead prediction from the state space model:

```
X̂_{t+k|t} = Aᵏ X̂_{t|t} + Σ_{j=0}^{k-1} Aʲ B u_{t+k-1-j}
Ŷ_{t+k|t} = C X̂_{t+k|t}
```

**Textbook:** Section 10.3.2, p. 296.

### 10.5 Signal Extraction and Common State Space Models

Examples of time series in state space form:
- **Local level model:** state = level, system noise = level innovations
- **Local linear trend model:** state = (level, slope)ᵀ
- **Seasonal models:** state includes seasonal harmonics

**Textbook:** Section 10.4, p. 299; Section 10.4.1, p. 301.

### 10.6 ML Estimation of State Space Models

The Kalman filter innovations `ν_t` are used to construct the likelihood function:

```
-2 ln L = Σ_t [ln |F_t| + νₜᵀ F_t⁻¹ ν_t]
```

where `F_t = C P_{t|t-1} Cᵀ + R` is the innovation variance.

**Textbook:** Section 10.6, p. 307.

---

## 11. Named Theorems, Definitions, and Rules — Master Index

Below is a comprehensive list of all named results from the course slides, with textbook cross-references.

---

### Theorem 5.10 — Yule-Walker Equations

**Statement:** For a stationary AR(p) process `φ(B) X_t = ε_t`, the autocovariances satisfy the system of linear equations:

```
γ(k) = φ₁ γ(k-1) + φ₂ γ(k-2) + ... + φ_p γ(k-p),    k ≥ 1
```

In matrix form: `Γ φ = γ` where `Γ_{ij} = γ(i-j)`.

This theorem is used for: (a) computing the theoretical ACF given AR parameters; (b) estimating AR parameters from the sample ACF (Yule-Walker estimator).

**Textbook:** Theorem 5.10, Section 5.5.2, p. 119–124.

**Slides:** Week 8.

---

### Definition: Strong Stationarity

A stochastic process `{X_t}` is strongly (strictly) stationary if all finite-dimensional distributions are time-invariant.

**Textbook:** Section 5.2.1.1, p. 99.

**Slides:** Week 4.

---

### Definition: Weak (Second-Order) Stationarity

A process is weakly stationary if mean, variance, and autocovariance function are time-invariant (autocovariance depends only on lag).

**Textbook:** Section 5.2.1.1, p. 99.

**Slides:** Week 4.

---

### Definition: Ergodicity (Mean-Ergodicity)

A stationary process is mean-ergodic if its time average converges to its ensemble mean. Sufficient condition: `Σ|γ(k)| < ∞`.

**Textbook:** Section 5.2.1, p. 99.

**Slides:** Week 4.

---

### Definition: Simple Exponential Smoothing (SES)

`S_N = (1-λ) Y_N + λ S_{N-1}`, smoothing constant `α = 1-λ`.

**Textbook:** Section 3.4.2, p. 50.

**Slides:** lect04.

---

### Definition: Total Memory T

`T = Σ_{j=0}^{N-1} λʲ = (1-λᴺ)/(1-λ)`. Measures effective number of observations used in WLS estimation.

**Textbook:** Not in textbook (lecture-specific definition). Mentioned in slides: lect04.

---

### Equation (5.65) — MA Autocovariance

```
γ(k) = σ² Σ_{j=0}^{q-k} θⱼ θ_{j+k},   0 ≤ k ≤ q;   γ(k) = 0 for k > q
```

**Textbook:** Eq. (5.65), Section 5.5.1, p. 117.

**Slides:** Week 8 (referenced as the formula for pure MA autocovariance).

---

### Equations (5.100) and (5.101) — ARMA Autocovariance

The autocovariance recursion for ARMA(p,q):
- For `k > q`: `γ(k) = φ₁ γ(k-1) + ... + φ_p γ(k-p)` — Eq. (5.100)
- For `k = 0, 1, ..., q`: initial condition equations involving MA parameters — Eq. (5.101)

**Textbook:** Eqs. (5.100) and (5.101), Section 5.5.3, p. 125.

**Slides:** Week 8.

---

### Stationarity Condition for ARMA

An ARMA(p,q) model `φ(B) X_t = θ(B) ε_t` is stationary iff all roots of `φ(z⁻¹) = 0` lie inside the unit circle.

**Textbook:** Section 5.5.3, p. 125.

**Slides:** Week 5.

---

### Invertibility Condition for ARMA

An ARMA(p,q) model is invertible iff all roots of `θ(z⁻¹) = 0` lie inside the unit circle.

**Textbook:** Section 5.5.1 (for MA), p. 117; Section 5.5.3 (for ARMA), p. 125.

**Slides:** Week 5.

---

### Multivariate Stationarity Condition (VARMA)

A VARMA process `Φ(B) Yₜ = Θ(B) εₜ` is stationary iff all roots of `det(Φ(z⁻¹)) = 0` lie strictly inside the unit circle.

**Textbook:** Section 9.3, p. 254.

**Slides:** Week 8.

---

### Multivariate Invertibility Condition (VARMA)

A VARMA process is invertible iff all roots of `det(Θ(z⁻¹)) = 0` lie strictly inside the unit circle.

**Textbook:** Section 9.3, p. 254.

**Slides:** Week 8.

---

### Spliid Method (1983)

Extended LS procedure for consistent estimation of multivariate ARMAX models.

**Reference:** Spliid, H. (1983). A fast estimation method for the vector autoregressive moving average model with exogenous variables. *Journal of the American Statistical Association*, 78(384), 843–849.

**Textbook:** Section 9.7.2, p. 271.

**Slides:** Week 8.

---

### Kalman Filter (Recursion)

Optimal linear state estimator for linear Gaussian state space models. Consists of prediction and update steps (see Section 10.3 above).

**Textbook:** Section 10.3.1, p. 289.

**Slides:** Week 9.

---

### Wold Decomposition Theorem

Every covariance-stationary process `{X_t}` can be written as:

```
X_t = Σ_{k=0}^{∞} ψ_k ε_{t-k} + V_t
```

where `ε_t` is white noise, `Σ ψ_k² < ∞`, `ψ₀ = 1`, and `V_t` is a linearly deterministic component uncorrelated with `ε_t`.

**Textbook:** Section 5.3.1, p. 107.

---

## Appendix: Textbook Section–Page Reference Table

| Topic | Textbook Section | Page |
|---|---|---|
| Regression model | 3.1 | 31 |
| General Linear Model (GLM) | 3.2 | 33 |
| OLS estimates | 3.2.1 | 34 |
| ML estimates | 3.2.2 | 40 |
| Prediction in GLM | 3.3.1 | 45 |
| Regression and exponential smoothing | 3.4 | 47 |
| Constant mean model (predictions) | 3.4.1 | 48 |
| Simple exponential smoothing | 3.4.2 | 50 |
| Prediction in trend models | 3.4.3 | 52 |
| Local trend and exponential smoothing | 3.4.4 | 56 |
| Seasonal time series / classical decomposition | 3.5.1 | 60 |
| Holt-Winters procedure | 3.5.2 | 61 |
| Global and local trend model (example) | 3.6 | 62 |
| Linear systems (time domain) | 4.1 | 70 |
| Linear systems (frequency domain) | 4.2 | 73 |
| z-transform | 4.4 | 80 |
| Frequently used operators (B, F, ∇, S) | 4.5 | 87 |
| Stochastic processes — introduction | 5.1 | 97 |
| Stochastic processes and their moments | 5.2 | 97 |
| Strong and weak stationarity | 5.2.1.1 | 99 |
| Covariance and correlation functions | 5.2.2 | 103 |
| Linear processes (discrete time) | 5.3.1 | 107 |
| Stationary processes in frequency domain | 5.4 | 113 |
| MA process | 5.5.1 | 117 |
| AR process | 5.5.2 | 119 |
| ARMA process | 5.5.3 | 125 |
| ARIMA process | 5.6.1 | 130 |
| Seasonal models | 5.6.2 | 132 |
| ARX / ARMAX models | 5.6.3 | 134 |
| Optimal prediction of stochastic processes | 5.7 | 135 |
| Estimation of covariance functions (ACF) | 6.2.1 | 146 |
| Identification of degree of differencing | 6.3.1 | 153 |
| Identification of ARMA part | 6.3.2 | 154 |
| Moment estimates (Yule-Walker) | 6.4.1 | 157 |
| LS estimator for linear dynamic models | 6.4.2 | 159 |
| Prediction error method | 6.4.3 | 163 |
| ML method for dynamic models | 6.4.4 | 166 |
| AIC/BIC model selection | 6.5.3 | 174 |
| Residual analysis | 6.6.2 | 175 |
| Transfer function models | 8.3.1 | 222 |
| Identification of transfer function models | 8.4 | 223 |
| Multivariate stationary processes | 9.1 | 249 |
| Multivariate ARMA process | 9.3 | 254 |
| Theoretical covariance matrix functions | 9.3.1 | 255 |
| Partial correlation matrix | 9.3.2 | 259 |
| VAR representation | 9.3.4 | 260 |
| Identification using pre-whitening | 9.6.1 | 269 |
| LS estimation (multivariate) | 9.7.1 | 270 |
| Spliid method (multivariate ARMAX) | 9.7.2 | 271 |
| ML estimates (multivariate) | 9.7.3 | 271 |
| Linear stochastic state space model | 10.1 | 284 |
| Transfer function / state space formulations | 10.2 | 286 |
| Kalman filter | 10.3.1 | 289 |
| k-step predictions (state space) | 10.3.2 | 296 |
| Common models in state space form | 10.4 | 299 |
| ML estimates of state space models | 10.6 | 307 |
| Recursive LS | 11.1 | 313 |
| Recursive LS with forgetting | 11.1.1 | 316 |
| Recursive pseudo-linear regression (RPLR) | 11.2 | 319 |
| Partial autocorrelations (appendix) | App. B | 357 |

---

*Document compiled from course slides: week2.pdf, week4.pdf, week5.pdf, week7.pdf, week8.pdf, week9.pdf, lect02_present.pdf, lect03_present.pdf, lect04_present.pdf. Textbook: Henrik Madsen, Time Series Analysis, Chapman & Hall/CRC, 2008.*
