#!/usr/bin/env python3
"""
kappa_estimation.py - Estimate κ and test κ/C₂ = 1

This script estimates the universal amplitude κ from the decay law:
  c_p(N) = κ / (φ(p) · ln N)

and tests the Constant Consistency Principle: κ = C₂

Author: Ruqing Chen
Paper: Self-Similar Structure in Goldbach Deviations (Paper V)
Repository: https://github.com/Ruqing1963/goldbach-fractal-analysis
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import argparse

# Twin prime constant
C2 = 0.6601618158468695739278121100145557784326233602847334133194484233354056

def estimate_coefficients(df, window_size=10000, step=5000):
    """
    Estimate c_3, c_5, c_7 using sliding window regression.
    
    Model: bias = c_0 + c_3·χ_3 + c_5·χ_5 + c_7·χ_7
    """
    results = []
    
    n_points = len(df)
    n_windows = (n_points - window_size) // step + 1
    
    for i in range(n_windows):
        start = i * step
        end = start + window_size
        
        if end > n_points:
            break
        
        window_df = df.iloc[start:end]
        
        # Regression
        X = window_df[['chi3', 'chi5', 'chi7']].values
        y = window_df['bias'].values
        
        model = LinearRegression()
        model.fit(X, y)
        
        c3, c5, c7 = model.coef_
        
        # Mean N and ln(N) for this window
        mean_N = window_df['N'].mean()
        mean_ln_N = window_df['ln_N'].mean()
        
        # Compute κ = c_p · φ(p) · ln(N)
        kappa3 = c3 * 2 * mean_ln_N  # φ(3) = 2
        kappa5 = c5 * 4 * mean_ln_N  # φ(5) = 4
        kappa7 = c7 * 6 * mean_ln_N  # φ(7) = 6
        
        results.append({
            'N_center': mean_N,
            'ln_N': mean_ln_N,
            'c3': c3,
            'c5': c5,
            'c7': c7,
            'kappa3': kappa3,
            'kappa5': kappa5,
            'kappa7': kappa7,
            'kappa_mean': (kappa3 + kappa5 + kappa7) / 3
        })
    
    return pd.DataFrame(results)

def test_kappa_equals_c2(kappa_values, c2=C2):
    """
    Statistical test: H₀: κ = C₂
    
    Uses z-test assuming approximate normality (CLT justified).
    """
    mean_kappa = np.mean(kappa_values)
    std_kappa = np.std(kappa_values, ddof=1)
    se_kappa = std_kappa / np.sqrt(len(kappa_values))
    
    # Ratio
    ratio = mean_kappa / c2
    ratio_error = se_kappa / c2
    
    # z-test
    z = (ratio - 1.0) / ratio_error
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))  # Two-tailed
    
    return {
        'kappa_mean': mean_kappa,
        'kappa_std': std_kappa,
        'kappa_se': se_kappa,
        'C2': c2,
        'ratio': ratio,
        'ratio_error': ratio_error,
        'z_statistic': z,
        'p_value': p_value
    }

def plot_kappa_analysis(coef_df, test_result, output_path):
    """Create κ analysis visualization."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: κ estimates over scale
    ax = axes[0]
    ax.plot(coef_df['N_center'], coef_df['kappa3'], 'o-', alpha=0.7, 
            label=r'$\kappa_3 = c_3 \cdot \phi(3) \cdot \ln N$')
    ax.plot(coef_df['N_center'], coef_df['kappa5'], 's-', alpha=0.7,
            label=r'$\kappa_5 = c_5 \cdot \phi(5) \cdot \ln N$')
    ax.plot(coef_df['N_center'], coef_df['kappa7'], '^-', alpha=0.7,
            label=r'$\kappa_7 = c_7 \cdot \phi(7) \cdot \ln N$')
    
    ax.axhline(C2, color='red', linestyle='--', linewidth=2, 
               label=f'$C_2 = {C2:.4f}$')
    ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
    
    ax.set_xlabel('N')
    ax.set_ylabel(r'$\kappa$')
    ax.set_title(r'(A) Local $\kappa$ Estimates')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel B: κ/C₂ ratio test
    ax = axes[1]
    
    ratio = test_result['ratio']
    ratio_err = test_result['ratio_error']
    
    ax.errorbar([1], [ratio], yerr=[ratio_err], fmt='o', markersize=12,
                capsize=5, capthick=2, color='steelblue', linewidth=2)
    ax.axhline(1.0, color='green', linestyle='--', linewidth=2, label='H₀: κ/C₂ = 1')
    
    # Confidence interval shading
    ax.axhspan(1.0 - 1.96 * ratio_err, 1.0 + 1.96 * ratio_err, 
               alpha=0.2, color='green', label='95% CI for H₀')
    
    ax.set_xlim(0.5, 1.5)
    ax.set_ylim(0.8, 1.3)
    ax.set_xticks([1])
    ax.set_xticklabels([r'$\kappa / C_2$'])
    ax.set_ylabel('Ratio')
    ax.set_title(f'(B) Constant Consistency Test\n' + 
                 f'κ/C₂ = {ratio:.3f} ± {ratio_err:.3f}, p = {test_result["p_value"]:.3f}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Kappa Estimation and Testing')
    parser.add_argument('--input', type=str, default='data/goldbach_counts.csv',
                        help='Input CSV file')
    parser.add_argument('--output', type=str, default='data/kappa_results.csv',
                        help='Output results file')
    parser.add_argument('--figure', type=str, default='figures/kappa_analysis.pdf',
                        help='Output figure file')
    parser.add_argument('--window', type=int, default=10000,
                        help='Window size for sliding regression')
    parser.add_argument('--step', type=int, default=5000,
                        help='Step size for sliding window')
    args = parser.parse_args()
    
    print("Kappa Estimation and Constant Consistency Test")
    print("=" * 50)
    
    # Load data
    print(f"Loading data from {args.input}...")
    df = pd.read_csv(args.input)
    print(f"  Loaded {len(df)} data points")
    
    # Estimate coefficients
    print(f"Estimating c_p coefficients (window={args.window}, step={args.step})...")
    coef_df = estimate_coefficients(df, window_size=args.window, step=args.step)
    print(f"  Generated {len(coef_df)} sliding window estimates")
    
    # Aggregate κ values
    all_kappas = np.concatenate([
        coef_df['kappa3'].values,
        coef_df['kappa5'].values,
        coef_df['kappa7'].values
    ])
    
    # Remove outliers (values too close to zero or negative)
    valid_kappas = all_kappas[all_kappas > 0.01]
    
    print()
    print("κ Summary Statistics:")
    print(f"  κ₃ mean: {coef_df['kappa3'].mean():.4f} ± {coef_df['kappa3'].std():.4f}")
    print(f"  κ₅ mean: {coef_df['kappa5'].mean():.4f} ± {coef_df['kappa5'].std():.4f}")
    print(f"  κ₇ mean: {coef_df['kappa7'].mean():.4f} ± {coef_df['kappa7'].std():.4f}")
    print(f"  Combined: {np.mean(valid_kappas):.4f} ± {np.std(valid_kappas):.4f}")
    
    # Statistical test
    print()
    print("Testing H₀: κ = C₂...")
    test_result = test_kappa_equals_c2(valid_kappas)
    
    print(f"  κ/C₂ = {test_result['ratio']:.4f} ± {test_result['ratio_error']:.4f}")
    print(f"  z = {test_result['z_statistic']:.2f}")
    print(f"  p-value = {test_result['p_value']:.4f}")
    print()
    
    if test_result['p_value'] > 0.05:
        print("  Result: Cannot reject H₀ at α = 0.05")
        print("  → Consistent with κ = C₂ (Constant Consistency Principle)")
    else:
        print("  Result: Reject H₀ at α = 0.05")
        print("  → Evidence against κ = C₂")
    
    # Save results
    coef_df.to_csv(args.output, index=False)
    
    # Save test summary
    test_summary = pd.DataFrame([test_result])
    test_summary.to_csv(args.output.replace('.csv', '_test.csv'), index=False)
    print(f"\nResults saved to {args.output}")
    
    # Generate figure
    print("Generating figure...")
    plot_kappa_analysis(coef_df, test_result, args.figure)
    
    print("\nDone!")

if __name__ == '__main__':
    main()
