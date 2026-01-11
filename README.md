# Goldbach Fractal Analysis

**Self-Similar Structure in Goldbach Deviations: L-Function Zeros and the Twin Prime Signature**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

This repository contains the code, data, and manuscript for **Paper V** in the Goldbach Conjecture Analysis Series. We present numerical evidence that deviations from the Hardy-Littlewood formula exhibit self-similar (fractal) structure governed by L-function zeros.

### Key Findings

| Result | Value | Significance |
|--------|-------|--------------|
| Hurst Exponent | H = 0.85 ± 0.02 | Long-range persistence (11σ above null) |
| κ/C₂ Ratio | 1.06 ± 0.06 | Consistent with unity (p = 0.32) |
| FFT Peak (ζ) | γ = 13.8 | 2.3% offset from γ₁ = 14.13 |

### The Constant Consistency Principle

We conjecture that the twin prime constant C₂ ≈ 0.6602 governs both:
- The **main term** amplitude in Hardy-Littlewood
- The **fluctuation envelope** of Goldbach deviations

## Repository Structure

```
goldbach-fractal-analysis/
├── README.md                 # This file
├── LICENSE                   # CC-BY-4.0 License
├── requirements.txt          # Python dependencies
├── paper/
│   ├── Chen_2026_Goldbach_Structure_L_Function_Zeros.pdf
│   └── paper_v_v10.tex       # LaTeX source
├── scripts/
│   ├── goldbach_count.py     # Exact G(N) computation
│   ├── riemann_radar.py      # FFT spectral analysis
│   ├── hurst_analysis.py     # R/S Hurst exponent
│   ├── kappa_estimation.py   # κ/C₂ ratio analysis
│   └── generate_figures.py   # Reproduce all figures
├── data/
│   ├── goldbach_counts.csv   # Precomputed G(N) values
│   ├── fft_results.csv       # Spectral analysis results
│   └── hurst_results.csv     # Hurst exponent data
├── figures/
│   ├── fig1_trinity.pdf
│   ├── fig2_riemann_radar.pdf
│   ├── fig3_hurst_improved.pdf
│   ├── fig4_kappa.pdf
│   ├── fig5_kappa_improved.pdf
│   └── fig6_fft_enhanced.pdf
└── docs/
    └── methodology.md        # Detailed algorithm descriptions
```

## Quick Start

### Installation

```bash
git clone https://github.com/Ruqing1963/goldbach-fractal-analysis.git
cd goldbach-fractal-analysis
pip install -r requirements.txt
```

### Run Analysis

```bash
# Compute Goldbach counts
python scripts/goldbach_count.py --nmax 1000000

# Riemann Radar FFT analysis
python scripts/riemann_radar.py --input data/goldbach_counts.csv

# Hurst exponent calculation
python scripts/hurst_analysis.py --input data/goldbach_counts.csv

# Reproduce all figures
python scripts/generate_figures.py
```

## Methodology

### 1. Riemann Radar (FFT in Log-Domain)

The key innovation is performing FFT in the logarithmic domain:

```
t = ln(N)  →  FFT  →  peaks at γ (L-function zeros)
```

Since L-function zeros produce oscillations cos(γ ln N), transforming to t = ln N converts these to pure frequencies detectable by FFT.

### 2. Hurst Exponent (R/S Analysis)

We use rescaled range (R/S) analysis to measure long-range dependence:

```
R/S(n) ~ n^H

H > 0.5  →  Persistent (long memory)
H = 0.5  →  Random walk
H < 0.5  →  Anti-persistent
```

Our result H = 0.85 indicates strong persistence in Goldbach deviations.

### 3. Constant Consistency Test

We test whether κ = C₂ using:

```
H₀: κ/C₂ = 1
Observed: 1.06 ± 0.06
z = 1.0, p = 0.32 → Cannot reject H₀
```

## Paper Series

This is Paper V in a series:

| Paper | Title | DOI |
|-------|-------|-----|
| I | Hardy-Littlewood Validated to N=10¹² | [10.5281/zenodo.18113330](https://doi.org/10.5281/zenodo.18113330) |
| II | Crossover Phenomenon | [10.5281/zenodo.18123132](https://doi.org/10.5281/zenodo.18123132) |
| III | Spectral Rigidity | [10.5281/zenodo.18148544](https://doi.org/10.5281/zenodo.18148544) |
| IV | Second Main Term | [10.5281/zenodo.18149305](https://doi.org/10.5281/zenodo.18149305) |
| **V** | **Self-Similar Structure** | **(this repository)** |

## Citation

```bibtex
@article{chen2026goldbach_v,
  title={Self-Similar Structure in Goldbach Deviations: 
         L-Function Zeros and the Twin Prime Signature},
  author={Chen, Ruqing},
  journal={Zenodo},
  year={2026},
  doi={10.5281/zenodo.XXXXXXXX},
  url={https://github.com/Ruqing1963/goldbach-fractal-analysis}
}
```

## License

This work is licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

- **Author:** Ruqing Chen
- **Email:** ruqing@hotmail.com
- **Institution:** GUT Geoservice Inc., Montreal, Canada

## Acknowledgments

The numerical computations were performed using Python with NumPy and SciPy. We thank the developers of these open-source tools.
