#!/usr/bin/env python3
"""
hurst_analysis.py - Hurst Exponent (R/S) Analysis

This script computes the Hurst exponent H for Goldbach deviation sequences
using rescaled range (R/S) analysis.

Interpretation:
  H > 0.5: Persistent (long-range positive correlations)
  H = 0.5: Random walk (no memory)
  H < 0.5: Anti-persistent (mean-reverting)

Author: Ruqing Chen
Paper: Self-Similar Structure in Goldbach Deviations (Paper V)
Repository: https://github.com/Ruqing1963/goldbach-fractal-analysis
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import argparse

# Twin prime constant
C2 = 0.6601618158468695739278121100145557784326233602847334133194484233354056

def rescaled_range(x, n):
    """
    Compute the rescaled range R/S for a time series x with window size n.
    
    R/S = (max(Y) - min(Y)) / S
    
    where Y is the cumulative deviation from mean, S is standard deviation.
    """
    if len(x) < n:
        return np.nan
    
    # Split into non-overlapping windows
    n_windows = len(x) // n
    rs_values = []
    
    for i in range(n_windows):
        window = x[i*n:(i+1)*n]
        
        # Mean and std
        mean = np.mean(window)
        std = np.std(window, ddof=1)
        
        if std == 0:
            continue
        
        # Cumulative deviations from mean
        deviations = window - mean
        cumsum = np.cumsum(deviations)
        
        # Range
        R = np.max(cumsum) - np.min(cumsum)
        
        # Rescaled range
        rs = R / std
        rs_values.append(rs)
    
    if len(rs_values) == 0:
        return np.nan
    
    return np.mean(rs_values)

def compute_hurst_exponent(x, min_window=10, max_window=None, n_windows=20):
    """
    Compute Hurst exponent using R/S analysis.
    
    H is estimated from: log(R/S) = H * log(n) + c
    """
    if max_window is None:
        max_window = len(x) // 4
    
    # Window sizes (logarithmically spaced)
    window_sizes = np.unique(np.logspace(
        np.log10(min_window), 
        np.log10(max_window), 
        n_windows
    ).astype(int))
    
    rs_values = []
    valid_sizes = []
    
    for n in window_sizes:
        rs = rescaled_range(x, n)
        if not np.isnan(rs) and rs > 0:
            rs_values.append(rs)
            valid_sizes.append(n)
    
    if len(valid_sizes) < 3:
        return np.nan, np.nan, [], []
    
    # Linear regression in log-log space
    log_n = np.log(valid_sizes)
    log_rs = np.log(rs_values)
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_n, log_rs)
    
    return slope, std_err, valid_sizes, rs_values

def generate_null_model(n, n_samples=100):
    """
    Generate null model: Poisson-distributed increments.
    
    This represents the null hypothesis of no long-range correlations.
    """
    hurst_values = []
    
    for _ in range(n_samples):
        # Random Poisson increments
        increments = np.random.poisson(10, n) - 10  # Zero-mean
        cumsum = np.cumsum(increments)
        
        h, _, _, _ = compute_hurst_exponent(cumsum)
        if not np.isnan(h):
            hurst_values.append(h)
    
    return np.mean(hurst_values), np.std(hurst_values)

def detrend_series(x, ln_N):
    """
    Remove the 1/ln(N) decay trend from the series.
    """
    trend = C2 / (2 * ln_N)
    return x - trend

def plot_hurst_analysis(sizes, rs_values, H, H_null, H_null_std, output_path):
    """Create Hurst analysis visualization."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: R/S plot
    ax = axes[0]
    log_n = np.log(sizes)
    log_rs = np.log(rs_values)
    
    ax.scatter(log_n, log_rs, s=50, alpha=0.7, label='Data')
    
    # Fit line
    slope, intercept = np.polyfit(log_n, log_rs, 1)
    fit_line = slope * log_n + intercept
    ax.plot(log_n, fit_line, 'r-', linewidth=2, label=f'Fit: H = {slope:.3f}')
    
    # Reference lines
    ax.plot(log_n, 0.5 * log_n + intercept, 'g--', alpha=0.5, label='H = 0.5 (Random)')
    
    ax.set_xlabel('log(n)')
    ax.set_ylabel('log(R/S)')
    ax.set_title('(A) Rescaled Range Analysis')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel B: H comparison
    ax = axes[1]
    categories = ['Observed\n(Detrended)', 'Null Model\n(Poisson)', 'Random Walk']
    values = [H, H_null, 0.5]
    errors = [0.02, H_null_std, 0]  # Approximate error for observed
    colors = ['steelblue', 'gray', 'lightgray']
    
    bars = ax.bar(categories, values, yerr=errors, capsize=5, color=colors, alpha=0.8)
    ax.axhline(0.5, color='green', linestyle='--', alpha=0.5, label='H = 0.5')
    ax.set_ylabel('Hurst Exponent H')
    ax.set_title('(B) Comparison with Null Model')
    ax.set_ylim(0, 1)
    
    # Add significance annotation
    if H > H_null + 3 * H_null_std:
        ax.annotate(f'{(H - H_null)/H_null_std:.0f}σ above null', 
                    xy=(0, H), xytext=(0.5, H + 0.1),
                    fontsize=12, fontweight='bold', color='red')
    
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Hurst Exponent Analysis')
    parser.add_argument('--input', type=str, default='data/goldbach_counts.csv',
                        help='Input CSV file')
    parser.add_argument('--output', type=str, default='data/hurst_results.csv',
                        help='Output results file')
    parser.add_argument('--figure', type=str, default='figures/hurst_analysis.pdf',
                        help='Output figure file')
    parser.add_argument('--null_samples', type=int, default=100,
                        help='Number of null model samples')
    args = parser.parse_args()
    
    print("Hurst Exponent (R/S) Analysis")
    print("=" * 50)
    
    # Load data
    print(f"Loading data from {args.input}...")
    df = pd.read_csv(args.input)
    print(f"  Loaded {len(df)} data points")
    
    # Extract bias series
    bias = df['bias'].values
    ln_N = df['ln_N'].values
    
    # Detrend
    print("Detrending signal...")
    detrended = detrend_series(bias, ln_N)
    
    # Compute Hurst exponent
    print("Computing Hurst exponent...")
    H, H_err, sizes, rs_values = compute_hurst_exponent(detrended)
    print(f"  H = {H:.4f} ± {H_err:.4f}")
    
    # Generate null model
    print(f"Generating null model ({args.null_samples} samples)...")
    H_null, H_null_std = generate_null_model(len(detrended), args.null_samples)
    print(f"  H_null = {H_null:.4f} ± {H_null_std:.4f}")
    
    # Significance test
    z_score = (H - H_null) / H_null_std
    print()
    print("Significance Test:")
    print(f"  z = (H - H_null) / σ_null = ({H:.3f} - {H_null:.3f}) / {H_null_std:.3f} = {z_score:.1f}")
    print(f"  → {abs(z_score):.0f}σ {'above' if z_score > 0 else 'below'} null hypothesis")
    
    # Save results
    results = pd.DataFrame({
        'metric': ['H_observed', 'H_error', 'H_null', 'H_null_std', 'z_score'],
        'value': [H, H_err, H_null, H_null_std, z_score]
    })
    results.to_csv(args.output, index=False)
    print(f"\nResults saved to {args.output}")
    
    # Generate figure
    print("Generating figure...")
    plot_hurst_analysis(sizes, rs_values, H, H_null, H_null_std, args.figure)
    
    print("\nInterpretation:")
    if H > 0.5:
        print("  H > 0.5 → PERSISTENT long-range correlations")
        print("  The Goldbach deviations exhibit fractal memory structure")
    else:
        print("  H ≤ 0.5 → No significant long-range persistence")
    
    print("\nDone!")

if __name__ == '__main__':
    main()
