#!/usr/bin/env python3
"""
riemann_radar.py - Riemann Radar: FFT Analysis in Log-Domain

This script performs spectral analysis of Goldbach deviations in the
logarithmic domain to detect L-function zero signatures.

Key Insight: Since L-function zeros produce oscillations cos(γ ln N),
transforming to t = ln N converts these to pure frequencies γ.

Author: Ruqing Chen
Paper: Self-Similar Structure in Goldbach Deviations (Paper V)
Repository: https://github.com/Ruqing1963/goldbach-fractal-analysis
"""

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import argparse

# Twin prime constant
C2 = 0.6601618158468695739278121100145557784326233602847334133194484233354056

# Known L-function zeros (first imaginary part γ₁)
KNOWN_ZEROS = {
    'ζ(s)': 14.1347,      # Riemann zeta
    'L(s,χ₃)': 8.0398,    # Dirichlet L-function mod 3
    'L(s,χ₅)': 6.0209,    # Dirichlet L-function mod 5
    'L(s,χ₇)': 5.1981,    # Dirichlet L-function mod 7
}

def load_data(filepath):
    """Load Goldbach count data from CSV."""
    df = pd.read_csv(filepath)
    return df

def uniform_log_sampling(df, n_points=1024):
    """
    Resample data uniformly in log-domain.
    
    This is crucial for FFT to correctly detect γ frequencies.
    """
    ln_N = df['ln_N'].values
    bias = df['bias'].values
    
    # Create uniform grid in ln(N)
    t_uniform = np.linspace(ln_N.min(), ln_N.max(), n_points)
    
    # Cubic spline interpolation
    spline = CubicSpline(ln_N, bias)
    bias_uniform = spline(t_uniform)
    
    return t_uniform, bias_uniform

def detrend_signal(t, signal):
    """
    Remove the 1/ln(N) decay trend.
    
    The decay law is: c_p(N) ~ κ/(φ(p) ln N)
    We subtract C₂/(2t) as the expected envelope.
    """
    trend = C2 / (2 * t)
    detrended = signal - trend
    return detrended

def compute_fft(t, signal):
    """
    Compute FFT and return frequency and power spectrum.
    
    The frequencies correspond to γ values of L-function zeros.
    """
    n = len(signal)
    dt = t[1] - t[0]  # Uniform spacing in ln(N)
    
    # FFT
    fft_vals = np.fft.fft(signal)
    power = np.abs(fft_vals[:n//2])**2
    
    # Frequency axis (corresponds to γ)
    freq = np.fft.fftfreq(n, dt)[:n//2]
    gamma = 2 * np.pi * freq  # Convert to γ
    
    return gamma, power

def find_spectral_peaks(gamma, power, height_threshold=None):
    """Find peaks in the power spectrum."""
    if height_threshold is None:
        height_threshold = np.mean(power) + 2 * np.std(power)
    
    peaks, properties = find_peaks(power, height=height_threshold)
    
    return gamma[peaks], power[peaks]

def match_zeros(detected_peaks, known_zeros, resolution=1.25):
    """
    Match detected peaks to known L-function zeros.
    
    Resolution is determined by: Δγ ≈ 2π/(ln N_max - ln N_min)
    """
    matches = []
    
    for name, gamma_known in known_zeros.items():
        # Find closest detected peak
        if len(detected_peaks) > 0:
            distances = np.abs(detected_peaks - gamma_known)
            min_idx = np.argmin(distances)
            gamma_detected = detected_peaks[min_idx]
            offset = abs(gamma_detected - gamma_known)
            rel_diff = offset / gamma_known * 100
            
            status = "✓ Match" if offset < resolution else "? Marginal"
            
            matches.append({
                'L-function': name,
                'γ_known': gamma_known,
                'γ_detected': gamma_detected,
                'offset': offset,
                'rel_diff_%': rel_diff,
                'status': status
            })
    
    return pd.DataFrame(matches)

def plot_riemann_radar(t, signal, detrended, gamma, power, known_zeros, output_path):
    """Create the Riemann Radar visualization."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    # Panel A: Original signal with trend
    ax = axes[0]
    ax.plot(t, signal, 'b-', alpha=0.7, linewidth=0.5)
    trend = C2 / (2 * t)
    ax.plot(t, trend, 'r--', linewidth=2, label=f'Trend: $C_2/(2t)$')
    ax.set_xlabel('$t = \\ln N$')
    ax.set_ylabel('Bias $(G-HL)/HL$')
    ax.set_title('(A) Raw Signal with Decay Trend')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel B: Detrended signal
    ax = axes[1]
    ax.plot(t, detrended, 'g-', alpha=0.7, linewidth=0.5)
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('$t = \\ln N$')
    ax.set_ylabel('Detrended Bias')
    ax.set_title('(B) Detrended Signal')
    ax.grid(True, alpha=0.3)
    
    # Panel C: Power spectrum with zero markers
    ax = axes[2]
    ax.semilogy(gamma, power, 'b-', linewidth=0.8)
    
    # Mark known zeros
    colors = ['red', 'orange', 'green', 'purple']
    for i, (name, gamma_known) in enumerate(known_zeros.items()):
        ax.axvline(gamma_known, color=colors[i % len(colors)], 
                   linestyle='--', alpha=0.7, label=f'{name}: γ={gamma_known:.2f}')
    
    ax.set_xlabel('$\\gamma$ (Imaginary part of zero)')
    ax.set_ylabel('Power')
    ax.set_title('(C) Riemann Radar: Power Spectrum')
    ax.set_xlim(0, 50)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved to {output_path}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Riemann Radar FFT Analysis')
    parser.add_argument('--input', type=str, default='data/goldbach_counts.csv', 
                        help='Input CSV file')
    parser.add_argument('--output', type=str, default='data/fft_results.csv',
                        help='Output results file')
    parser.add_argument('--figure', type=str, default='figures/riemann_radar.pdf',
                        help='Output figure file')
    parser.add_argument('--n_points', type=int, default=2048,
                        help='Number of points for FFT')
    args = parser.parse_args()
    
    print("Riemann Radar: FFT Analysis in Log-Domain")
    print("=" * 50)
    
    # Load data
    print(f"Loading data from {args.input}...")
    df = load_data(args.input)
    print(f"  Loaded {len(df)} data points")
    print(f"  N range: [{df['N'].min()}, {df['N'].max()}]")
    
    # Compute resolution
    ln_range = df['ln_N'].max() - df['ln_N'].min()
    resolution = 2 * np.pi / ln_range
    print(f"  FFT resolution: Δγ ≈ {resolution:.3f}")
    print()
    
    # Uniform sampling in log-domain
    print("Resampling uniformly in ln(N)...")
    t, signal = uniform_log_sampling(df, n_points=args.n_points)
    print(f"  Resampled to {len(t)} points")
    
    # Detrend
    print("Detrending signal...")
    detrended = detrend_signal(t, signal)
    
    # FFT
    print("Computing FFT...")
    gamma, power = compute_fft(t, detrended)
    
    # Find peaks
    print("Finding spectral peaks...")
    peak_gammas, peak_powers = find_spectral_peaks(gamma, power)
    print(f"  Found {len(peak_gammas)} significant peaks")
    
    # Match to known zeros
    print()
    print("Matching to known L-function zeros:")
    print("-" * 60)
    matches = match_zeros(peak_gammas, KNOWN_ZEROS, resolution=resolution)
    print(matches.to_string(index=False))
    
    # Save results
    matches.to_csv(args.output, index=False)
    print(f"\nResults saved to {args.output}")
    
    # Generate figure
    print(f"\nGenerating figure...")
    plot_riemann_radar(t, signal, detrended, gamma, power, KNOWN_ZEROS, args.figure)
    
    print("\nDone!")

if __name__ == '__main__':
    main()
