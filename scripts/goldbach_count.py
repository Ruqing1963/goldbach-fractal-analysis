#!/usr/bin/env python3
"""
goldbach_count.py - Compute exact Goldbach representation counts G(N)

This script computes the number of ways to express even integers N as
the sum of two primes: G(N) = |{(p,q) : p + q = N, p ≤ q, p,q prime}|

Author: Ruqing Chen
Paper: Self-Similar Structure in Goldbach Deviations (Paper V)
Repository: https://github.com/Ruqing1963/goldbach-fractal-analysis
"""

import numpy as np
import pandas as pd
from scipy.special import expi
import argparse
from tqdm import tqdm
import time

# Twin prime constant C_2
C2 = 0.6601618158468695739278121100145557784326233602847334133194484233354056

def sieve_of_eratosthenes(limit):
    """Generate all primes up to limit using Sieve of Eratosthenes."""
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[0:2] = False
    for i in range(2, int(np.sqrt(limit)) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    return np.where(is_prime)[0]

def goldbach_count(n, primes, prime_set):
    """
    Count the number of Goldbach representations for even integer n.
    
    Returns G(n) = number of ways to write n = p + q where p ≤ q are primes.
    """
    if n % 2 != 0 or n < 4:
        return 0
    
    count = 0
    half = n // 2
    for p in primes:
        if p > half:
            break
        if (n - p) in prime_set:
            count += 1
    return count

def singular_series(n, primes):
    """
    Compute the singular series S(N) for Hardy-Littlewood formula.
    
    S(N) = ∏_{p|N, p>2} (p-1)/(p-2)
    """
    s = 1.0
    temp = n
    # Remove factor of 2
    while temp % 2 == 0:
        temp //= 2
    
    # Check odd prime factors
    for p in primes:
        if p > temp:
            break
        if p == 2:
            continue
        if temp % p == 0:
            s *= (p - 1) / (p - 2)
            while temp % p == 0:
                temp //= p
    
    # If remaining factor > 1, it's a prime
    if temp > 1:
        s *= (temp - 1) / (temp - 2)
    
    return s

def hardy_littlewood_prediction(n, primes):
    """
    Compute Hardy-Littlewood prediction for G(N).
    
    HL(N) = 2 * C_2 * S(N) * Li_2(N)
    
    where Li_2(N) = ∫_2^N dt/(ln t)^2
    """
    if n < 4:
        return 0
    
    s_n = singular_series(n, primes)
    
    # Li_2(N) approximation using logarithmic integral
    ln_n = np.log(n)
    li2 = n / (ln_n ** 2) * (1 + 2/ln_n + 6/(ln_n**2))
    
    return 2 * C2 * s_n * li2

def compute_batch(n_start, n_end, primes, prime_set, step=2):
    """Compute G(N) and HL(N) for a range of even integers."""
    results = []
    
    for n in tqdm(range(n_start, n_end + 1, step), desc=f"Computing N={n_start}-{n_end}"):
        if n % 2 != 0:
            continue
        
        g_n = goldbach_count(n, primes, prime_set)
        hl_n = hardy_littlewood_prediction(n, primes)
        s_n = singular_series(n, primes)
        
        # Compute bias
        if hl_n > 0:
            bias = (g_n - hl_n) / hl_n
        else:
            bias = 0
        
        # Compute Dirichlet characters
        chi3 = 1 if n % 3 == 1 else (-1 if n % 3 == 2 else 0)
        chi5 = 1 if n % 5 in [1, 4] else (-1 if n % 5 in [2, 3] else 0)
        chi7 = 1 if n % 7 in [1, 2, 4] else (-1 if n % 7 in [3, 5, 6] else 0)
        
        results.append({
            'N': n,
            'G_N': g_n,
            'HL_N': hl_n,
            'S_N': s_n,
            'bias': bias,
            'ln_N': np.log(n),
            'chi3': chi3,
            'chi5': chi5,
            'chi7': chi7
        })
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description='Compute Goldbach counts G(N)')
    parser.add_argument('--nmin', type=int, default=1000, help='Minimum N')
    parser.add_argument('--nmax', type=int, default=100000, help='Maximum N')
    parser.add_argument('--output', type=str, default='data/goldbach_counts.csv', help='Output file')
    args = parser.parse_args()
    
    print(f"Goldbach Count Computation")
    print(f"=" * 50)
    print(f"Range: N ∈ [{args.nmin}, {args.nmax}]")
    print(f"Output: {args.output}")
    print()
    
    # Generate primes
    print("Generating primes...")
    start_time = time.time()
    primes = sieve_of_eratosthenes(args.nmax)
    prime_set = set(primes)
    print(f"  π({args.nmax}) = {len(primes):,} primes")
    print(f"  Time: {time.time() - start_time:.2f}s")
    print()
    
    # Compute Goldbach counts
    print("Computing Goldbach counts...")
    start_time = time.time()
    df = compute_batch(args.nmin, args.nmax, primes, prime_set)
    print(f"  Computed {len(df):,} values")
    print(f"  Time: {time.time() - start_time:.2f}s")
    print()
    
    # Save results
    df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")
    
    # Print summary statistics
    print()
    print("Summary Statistics:")
    print(f"  Mean G(N): {df['G_N'].mean():.2f}")
    print(f"  Mean bias: {df['bias'].mean():.6f}")
    print(f"  Std bias:  {df['bias'].std():.6f}")

if __name__ == '__main__':
    main()
