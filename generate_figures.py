#!/usr/bin/env python3
"""
generate_figures.py - Reproduce all figures for Paper V

This script regenerates all figures used in the paper from the raw data.

Author: Ruqing Chen
Paper: Self-Similar Structure in Goldbach Deviations (Paper V)
Repository: https://github.com/Ruqing1963/goldbach-fractal-analysis
"""

import subprocess
import sys
import os

def run_script(script_name, args=[]):
    """Run a Python script with arguments."""
    cmd = [sys.executable, f'scripts/{script_name}'] + args
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    else:
        print(result.stdout)
    return result.returncode == 0

def main():
    print("=" * 60)
    print("Paper V Figure Generation")
    print("=" * 60)
    print()
    
    # Check if data exists
    if not os.path.exists('data/goldbach_counts.csv'):
        print("Data file not found. Generating...")
        success = run_script('goldbach_count.py', [
            '--nmin', '1000',
            '--nmax', '500000',
            '--output', 'data/goldbach_counts.csv'
        ])
        if not success:
            print("Failed to generate data. Exiting.")
            return
        print()
    
    # Generate Riemann Radar figure
    print("=" * 60)
    print("Generating Riemann Radar (FFT) Figure...")
    print("=" * 60)
    run_script('riemann_radar.py', [
        '--input', 'data/goldbach_counts.csv',
        '--output', 'data/fft_results.csv',
        '--figure', 'figures/fig2_riemann_radar.pdf'
    ])
    print()
    
    # Generate Hurst figure
    print("=" * 60)
    print("Generating Hurst Exponent Figure...")
    print("=" * 60)
    run_script('hurst_analysis.py', [
        '--input', 'data/goldbach_counts.csv',
        '--output', 'data/hurst_results.csv',
        '--figure', 'figures/fig3_hurst.pdf'
    ])
    print()
    
    # Generate Kappa figure
    print("=" * 60)
    print("Generating Kappa Analysis Figure...")
    print("=" * 60)
    run_script('kappa_estimation.py', [
        '--input', 'data/goldbach_counts.csv',
        '--output', 'data/kappa_results.csv',
        '--figure', 'figures/fig4_kappa.pdf'
    ])
    print()
    
    print("=" * 60)
    print("All figures generated successfully!")
    print("=" * 60)
    print()
    print("Figures saved to figures/ directory:")
    for f in os.listdir('figures'):
        if f.endswith('.pdf'):
            print(f"  - {f}")

if __name__ == '__main__':
    main()
