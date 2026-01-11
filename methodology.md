# Methodology

Detailed algorithm descriptions for Paper V analysis.

## 1. Goldbach Count Computation

### Algorithm

For each even integer N ≥ 4:

```
G(N) = |{(p, q) : p + q = N, p ≤ q, p and q are prime}|
```

**Implementation:**
1. Generate all primes up to N using Sieve of Eratosthenes
2. For each prime p ≤ N/2, check if (N - p) is also prime
3. Count valid pairs

**Complexity:** O(N / ln N) per N

### Hardy-Littlewood Prediction

```
HL(N) = 2 · C₂ · S(N) · Li₂(N)
```

Where:
- C₂ ≈ 0.6602 is the twin prime constant
- S(N) = ∏_{p|N, p>2} (p-1)/(p-2) is the singular series
- Li₂(N) = ∫₂^N dt/(ln t)² is the logarithmic integral

---

## 2. Riemann Radar (FFT Analysis)

### Key Insight

L-function zeros produce oscillations of the form:

```
cos(γ · ln N)
```

Transforming to t = ln N converts these to pure frequencies γ.

### Algorithm

1. **Uniform Sampling in Log-Domain**
   - Create uniform grid: t_i = ln(N_min) + i · Δt
   - Interpolate bias values using cubic spline

2. **Detrending**
   - Remove decay trend: signal - C₂/(2t)

3. **FFT**
   - Apply numpy.fft.fft to detrended signal
   - Convert frequency to γ: γ = 2π · f

4. **Peak Detection**
   - Find peaks above threshold (mean + 2σ)
   - Match to known L-function zeros

### Resolution

The frequency resolution is:

```
Δγ = 2π / (ln N_max - ln N_min)
```

For N ∈ [10³, 1.5×10⁶]: Δγ ≈ 1.25

---

## 3. Hurst Exponent (R/S Analysis)

### Definition

For a time series {x_i}, the Hurst exponent H is defined by:

```
R/S(n) ~ n^H
```

Where R/S is the rescaled range.

### Algorithm

For window size n:

1. **Partition** series into non-overlapping windows of size n
2. For each window:
   - Compute mean μ and standard deviation σ
   - Compute cumulative deviations: Y_k = Σᵢ₌₁ᵏ (x_i - μ)
   - Compute range: R = max(Y) - min(Y)
   - Compute rescaled range: R/S = R / σ
3. **Average** R/S over all windows
4. **Repeat** for multiple n values (log-spaced)
5. **Fit** linear regression: log(R/S) = H · log(n) + c

### Interpretation

- H > 0.5: Persistent (positive correlations)
- H = 0.5: Random walk
- H < 0.5: Anti-persistent

---

## 4. κ/C₂ Statistical Test

### Decay Law

The Dirichlet correction coefficients follow:

```
c_p(N) = κ / (φ(p) · ln N) + O(1/ln²N)
```

### Estimation

1. Use sliding window regression to estimate c_3, c_5, c_7
2. Compute κ_p = c_p · φ(p) · ln(N) for each window
3. Aggregate all κ estimates

### Hypothesis Test

- H₀: κ/C₂ = 1 (Constant Consistency Principle)
- H₁: κ/C₂ ≠ 1

Test statistic:

```
z = (κ̄/C₂ - 1) / SE(κ̄/C₂)
```

p-value from standard normal distribution (two-tailed).

---

## 5. Null Model

### Poisson-Hardy-Littlewood Model

The null hypothesis assumes:
- G(N) follows Hardy-Littlewood prediction exactly
- Deviations are i.i.d. Poisson noise

### Generation

1. Generate Poisson increments with λ = E[G(N)]
2. Compute Hurst exponent for null series
3. Repeat 100+ times to establish null distribution

### Expected Values

- H_null ≈ 0.5 (random walk)
- No spectral peaks at L-function zeros

---

## References

1. Hurst, H.E. (1951). Long-term storage capacity of reservoirs.
2. Mandelbrot, B.B. & Van Ness, J.W. (1968). Fractional Brownian motions.
3. Hardy, G.H. & Littlewood, J.E. (1923). Some problems of 'Partitio Numerorum' III.
