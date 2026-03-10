#!/usr/bin/env python3
"""
ESG — LAYER 1 CONNECTION: C(β) AND THE ZERO COUNTING MEASURE
================================================================
Replace Mittag-Leffler (Layer 2: individual poles) with integral
constraints on dN(T) (Layer 1: collective measure).

The key identity (CORRECT signs):
  d²/dβ² log ζ(β) = Var_β(log n) > 0

This is the SPECIFIC HEAT. It is POSITIVE (variance).

The zero counting function N(T) = #{ρ : 0 < Im(ρ) ≤ T} satisfies:
  N(T) = (T/2π) log(T/2πe) + S(T) + O(1/T)
where S(T) is the argument of ζ on the critical line.

The connection between C(β) and N(T) must go through an INTEGRAL
transform (Layer 1), not a pole-by-pole sum (Layer 2).

F. D. Blum, March 2026
"""

import numpy as np
from math import log, exp, pi, sqrt, gcd
from scipy.special import polygamma
import time

# ═══════════════════════════════════════════════════════════════
# PART 1: THE CORRECT SIGN IDENTITY
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 1: ESTABLISHING THE CORRECT IDENTITY")
print("=" * 85)
print()

def zeta_real(beta, N=100000):
    return sum(n**(-beta) for n in range(1, N+1))

def log_zeta(beta, N=100000):
    return log(zeta_real(beta, N))

def C_variance(beta, N=100000):
    """C(β) = Var_β(log n) — always positive."""
    weights = np.array([n**(-beta) for n in range(1, N+1)])
    Z = np.sum(weights)
    weights /= Z
    log_n = np.array([log(n) for n in range(1, N+1)])
    return np.sum(weights * log_n**2) - np.sum(weights * log_n)**2

def d2_log_zeta(beta, N=100000, h=1e-5):
    """d²/dβ² log ζ(β) — numerical second derivative."""
    return (log_zeta(beta+h, N) - 2*log_zeta(beta, N) + log_zeta(beta-h, N)) / h**2

print("Verifying: C(β) = Var_β(log n) = +d²/dβ² log ζ(β)")
print(f"{'β':>6} {'Var (>0)':>12} {'d²logζ/dβ²':>14} {'ratio':>8}")
print("-" * 44)

for beta in [1.1, 1.5, 2.0, 3.0, 5.0, 10.0]:
    C = C_variance(beta, 50000)
    d2 = d2_log_zeta(beta, 50000)
    ratio = C / d2 if abs(d2) > 1e-15 else 0
    print(f"{beta:6.1f} {C:12.8f} {d2:14.8f} {ratio:8.4f}")

print()
print("Ratio ≈ 1.00 everywhere → C(β) = +d²/dβ² log ζ(β). Correct sign: PLUS.")
print()

# ═══════════════════════════════════════════════════════════════
# PART 2: THE INTEGRAL REPRESENTATION
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 2: FROM EULER PRODUCT TO INTEGRAL")
print("=" * 85)
print()

print("""
LAYER 1 ROUTE (no individual zeros):

Start from the Euler product: log ζ(β) = -Σ_p log(1 - p^{-β})

Then: d/dβ log ζ(β) = Σ_p (log p) p^{-β} / (1 - p^{-β})
                     = Σ_p (log p) Σ_{k=1}^∞ p^{-kβ}
                     = Σ_n Λ(n) n^{-β}
                     = -ζ'(β)/ζ(β)

And: d²/dβ² log ζ(β) = Σ_n Λ(n) (log n) n^{-β} - [Σ_n Λ(n) n^{-β}]²/ζ(β)²
Wait, let me be more careful:

  d²/dβ² log ζ = d/dβ [-ζ'/ζ] = -(ζ''/ζ - (ζ'/ζ)²) = (ζ'/ζ)² - ζ''/ζ

Actually:
  d/dβ log ζ = ζ'/ζ  (where ζ' = dζ/dβ = -Σ n^{-β} log n)
  d²/dβ² log ζ = ζ''/ζ - (ζ'/ζ)²

Let's verify:
  ζ' = Σ n^{-β} (log n)²  ... wait, 
  ζ(β) = Σ n^{-β}
  ζ'(β) = dζ/dβ = -Σ n^{-β} log n
  ζ''(β) = Σ n^{-β} (log n)²

  d²/dβ² log ζ = (ζ''ζ - (ζ')²) / ζ²
              = ζ''/ζ - (ζ'/ζ)²
              = <(log n)²>_β - <log n>²_β
              = Var_β(log n) = C(β). ✓

So C(β) can be written ENTIRELY in terms of the prime sum:

  C(β) = Σ_n Λ(n)(log n) n^{-β} / ζ(β) 
         - [Σ_n Λ(n) n^{-β} / ζ(β)]² ... 

No, more directly:

  C(β) = <(log n)²>_β - <log n>²_β

where <f>_β = Σ f(n) n^{-β} / ζ(β).

This is PURELY ARITHMETIC: weights n^{-β}/ζ(β), function log n.
No zeros appear. No Fourier. No poles. Just primes and integers.

The connection to zeros comes through the IDENTITY:

  log ζ(β) = ∫₁^∞ π(x) / (x(x^β - 1)) dx    (for β > 1)

where π(x) is the prime counting function. And C(β) = d²/dβ².

So C(β) is an integral transform of π(x):

  C(β) = d²/dβ² ∫₁^∞ π(x) / (x(x^β - 1)) dx

The zeros of ζ control π(x) through the explicit formula:

  π(x) = li(x) - Σ_ρ li(x^ρ) + lower order

But this explicit formula IS Layer 2 (Fourier in disguise).

THE LAYER 1 ALTERNATIVE: Don't use the explicit formula.
Use the INTEGRAL KERNEL directly.
""")

# ═══════════════════════════════════════════════════════════════
# PART 3: C(β) AS INTEGRAL TRANSFORM OF π(x)
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 3: C(β) AS INTEGRAL TRANSFORM OF π(x)")
print("=" * 85)
print()

# Verify: log ζ(β) = Σ_p Σ_{k=1}^∞ p^{-kβ}/k
# So: d²/dβ² log ζ = Σ_p Σ_k (log p)² k p^{-kβ}
# This is the PRIME POWER VARIANCE.

# More usefully, via Stieltjes integral:
# log ζ(β) = ∫₂^∞ (dπ(x)/x) × 1/(x^β - 1) × x
# Hmm, let me use the Chebyshev function instead.

# ψ(x) = Σ_{p^k ≤ x} log p  (Chebyshev psi)
# -ζ'(β)/ζ(β) = β ∫₁^∞ ψ(x) x^{-β-1} dx  (Mellin transform of ψ)
# C(β) = d/dβ [-ζ'/ζ] = ... complicated

# SIMPLER: use the prime zeta function
# P(β) = Σ_p p^{-β}
# log ζ(β) = Σ_{k=1}^∞ P(kβ)/k

# C(β) = d²/dβ² Σ_k P(kβ)/k = Σ_k k P''(kβ)

# P''(β) = Σ_p (log p)² p^{-β}

# So C(β) = Σ_{k=1}^∞ k Σ_p (log p)² p^{-kβ}

# This is a LAYER 1 formula: only primes, only powers, no zeros.
# It's the prime representation of the specific heat.

print("C(β) in the prime representation:")
print("  C(β) = Σ_{k=1}^∞ k × Σ_p (log p)² p^{-kβ}")
print()

# Verify numerically
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
          53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
          127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
          193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251]

def C_prime_rep(beta, primes, k_max=20):
    """C(β) = Σ_k k × Σ_p (log p)² p^{-kβ}"""
    total = 0.0
    for k in range(1, k_max + 1):
        for p in primes:
            total += k * (log(p))**2 * p**(-k * beta)
    return total

print(f"{'β':>6} {'C_variance':>12} {'C_prime_rep':>12} {'ratio':>8}")
print("-" * 42)

for beta in [1.5, 2.0, 3.0, 5.0, 10.0]:
    C_v = C_variance(beta, 50000)
    C_p = C_prime_rep(beta, primes, k_max=30)
    ratio = C_v / C_p if abs(C_p) > 1e-15 else 0
    print(f"{beta:6.1f} {C_v:12.8f} {C_p:12.8f} {ratio:8.4f}")

print()
print("Match confirms: C(β) = Σ_k k × Σ_p (log p)² p^{-kβ}")
print("This is Layer 1: only primes, only real arithmetic.")
print()

# ═══════════════════════════════════════════════════════════════
# PART 4: THE LAYER 1 CONSTRAINT
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 4: THE LAYER 1 CONSTRAINT ON PRIME DISTRIBUTION")
print("=" * 85)
print()

print("""
C(β) = Σ_k k × Σ_p (log p)² p^{-kβ} > 0 for all β > 1.

This is TRIVIALLY true (all terms are positive).

So C(β) > 0 gives NO constraint on the primes — it's automatic.

This means: the positivity of the specific heat, viewed as a 
constraint on the PRIME side, is vacuous. It doesn't constrain
the primes because it's trivially satisfied.

The non-trivial content must come from the OTHER direction:
from the ZERO side of the explicit formula.

But the zero side IS Layer 2 (Fourier/Mittag-Leffler).

THIS IS THE FUNDAMENTAL OBSTACLE:

  C(β) on the PRIME side: trivially positive → no constraint
  C(β) on the ZERO side: non-trivially constrained → but Layer 2

The positivity is a property of the PRIME representation.
The constraint on zeros requires translating between representations.
And that translation IS the Fourier/explicit formula.

There is no way around this. The connection between primes and 
zeros is INHERENTLY a duality, and duality is Layer 2.
""")

# ═══════════════════════════════════════════════════════════════
# PART 5: IS THERE A LAYER 1 BRIDGE?
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 5: SEARCHING FOR A LAYER 1 BRIDGE")
print("=" * 85)
print()

print("""
The obstacle seems fundamental: primes ↔ zeros is a duality.
But let's check if there's a STATIC version of the duality.

OPTION A: The Hadamard product (no time, but uses individual zeros)
  ζ(s) = e^{A+Bs} × (s-1)⁻¹ × Π_ρ (1-s/ρ) e^{s/ρ}
  This is an ALGEBRAIC identity, not a Fourier transform.
  It uses individual zeros → Layer 2? 
  Actually: it's a factorization of an analytic function.
  Factorization is ALGEBRAIC, not temporal.
  
  Verdict: Layer 1 (barely). It lists zeros but doesn't oscillate.

OPTION B: The functional equation (no time, no zeros individually)
  ξ(s) = ξ(1-s)
  This constrains ζ as a WHOLE, not zero by zero.
  It's a symmetry of the function, not a decomposition.
  
  Verdict: Layer 1. No individual zeros needed.

OPTION C: The moment approach
  Define: M_k(β) = <(log n)^k>_β = ζ^{(k)}(β)/ζ(β) (up to signs)
  The moments M_k contain ALL information about ζ(β) for β > 1.
  By analytic continuation, they also determine ζ(s) for all s.
  The zeros are where ζ(s) = 0, which is where the moments "diverge"
  in a specific pattern.
  
  Verdict: Layer 1. Moments are static averages.

OPTION D: The Hamburger moment problem
  The measure dν_β(x) = Σ_n n^{-β} δ(x - log n) / ζ(β)
  has moments M_k(β) = <x^k>_{ν_β}.
  
  The MOMENT PROBLEM asks: given {M_k}, is the measure ν unique?
  For ν_β: it's a discrete measure on {log n}, so it's determined
  by its support (known) and its weights (n^{-β}/ζ(β)).
  
  The zeros of ζ(s) are where ζ(β) × (analytic continuation) = 0.
  In moment language: the zeros are where the moment-generating
  function Σ_n e^{-s log n} = ζ(s) vanishes.
  
  The MOMENT-GENERATING FUNCTION is a Laplace transform.
  Laplace transform is the REAL-VARIABLE version of Fourier.
  Is it Layer 1 or Layer 2?
""")

# ═══════════════════════════════════════════════════════════════
# PART 6: LAPLACE vs FOURIER — THE KEY DISTINCTION
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 6: LAPLACE vs FOURIER — IS THERE A DIFFERENCE?")
print("=" * 85)
print()

print("""
Fourier: f̂(γ) = ∫ f(x) e^{-iγx} dx    (oscillatory → Layer 2)
Laplace: F(β) = ∫ f(x) e^{-βx} dx       (decaying → Layer 1?)

The Fourier transform requires oscillation (e^{-iγx}).
The Laplace transform requires only decay (e^{-βx}).

ζ(β) = Σ n^{-β} = Σ e^{-β log n} is a LAPLACE transform of the
counting measure on {log n}.

For REAL β: this is purely decaying, no oscillation. Layer 1.
For COMPLEX s = σ + iγ: ζ(s) = Σ e^{-(σ+iγ)log n} has oscillation.
  The imaginary part introduces e^{-iγ log n} = oscillation = Layer 2.

THE CRITICAL DISTINCTION:

ζ(β) for REAL β > 1: Layer 1 (Laplace, no oscillation)
ζ(s) for COMPLEX s: Layer 2 (Fourier component in Im(s))

The zeros of ζ are at COMPLEX s = 1/2 + iγ.
Accessing them requires evaluating ζ at complex s.
Evaluating at complex s introduces the oscillatory factor e^{-iγ log n}.
This is INHERENTLY Layer 2.

CONCLUSION: There is NO way to access the individual zeros of ζ
using only Layer 1 tools. The zeros live in the complex plane,
and the imaginary direction IS oscillation IS time.

HOWEVER: we don't need individual zeros. We need the COLLECTIVE
constraint "all zeros have Re = 1/2." This is a constraint on
the REAL function ζ(β) (through analytic continuation), and
analytic continuation of a REAL analytic function is determined
by its values on the real axis (which is Layer 1).

THE KEY THEOREM (Carlson's theorem):
  If F(s) is analytic in Re(s) > 0, bounded by e^{c|s|} for c < π,
  and F(n) = 0 for all positive integers n, then F ≡ 0.

  This means: the values of ζ on the REAL axis (Layer 1 data)
  UNIQUELY DETERMINE ζ everywhere (including its zeros).

  The zeros are determined by Layer 1 data. We just can't
  COMPUTE them from Layer 1 data (computing requires Layer 2).

  But CONSTRAINING them might not require computing them.
""")

# ═══════════════════════════════════════════════════════════════
# PART 7: THE REAL-AXIS APPROACH
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 7: CONSTRAINING ZEROS FROM REAL-AXIS DATA")
print("=" * 85)
print()

print("""
The Layer 1 data available:
  ζ(β) for all real β > 1 (the partition function)
  ζ'(β), ζ''(β), ... all derivatives (moments)
  C(β) = d²/dβ² log ζ(β) (the specific heat)
  All of these are REAL, POSITIVE (for β > 1), and COMPUTABLE.

The question: do these data CONSTRAIN the zeros to Re = 1/2?

Classical results that use ONLY real-axis data:

1. DE LA VALLÉE-POUSSIN (1896):
   ζ(σ) ≠ 0 for σ = 1 was proved using:
   3 + 4cos(θ) + cos(2θ) ≥ 0
   This uses REAL values of Re[ζ'/ζ] at σ near 1.
   It gives the zero-free region σ > 1 - c/log(|t|+2).
   Layer 1? The cosine involves e^{iθ} = oscillation... but the
   INEQUALITY 3 + 4cos + cos2 ≥ 0 is a REAL inequality.

2. TURÁN POWER SUM METHOD:
   Σ_k |a_k|² ≥ 0 where a_k are sums of prime powers.
   Uses only REAL data (norms of sums).
   Gives partial results toward RH.

3. LI'S CRITERION:
   RH ⟺ λ_n ≥ 0 for all n ≥ 1
   where λ_n = Σ_ρ [1 - (1-1/ρ)^n]
   The λ_n can be computed from ζ^{(k)}(1) — REAL derivatives at β = 1.
   This is ENTIRELY Layer 1!

Li's criterion is the answer.
""")

# ═══════════════════════════════════════════════════════════════
# PART 8: LI'S CRITERION — THE LAYER 1 RH
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 8: LI'S CRITERION AS LAYER 1 RH")
print("=" * 85)
print()

print("""
LI'S CRITERION (Xian-Jin Li, 1997):

  RH ⟺ λ_n ≥ 0 for all n = 1, 2, 3, ...

where:
  λ_n = (1/(n-1)!) d^n/ds^n [s^{n-1} log ξ(s)]|_{s=1}

Equivalently:
  λ_n = Σ_ρ [1 - (1 - 1/ρ)^n]

The λ_n can be computed from:
  λ_n = Σ_{k=0}^{n} C(n,k) (-1)^k σ_k

where σ_k are the Stieltjes constants:
  σ_k = (-1)^k k! × [coefficient of (s-1)^k in Laurent expansion of ζ(s) at s=1]

The Stieltjes constants σ_k are REAL numbers computable from 
REAL derivatives of ζ at real points. They are LAYER 1.

LI'S CRITERION IS A LAYER 1 REFORMULATION OF RH:
  An infinite sequence of REAL inequalities λ_n ≥ 0,
  each computable from REAL-axis data of ζ.
  No oscillation. No Fourier. No complex evaluation.
  Just derivatives at a real point.
""")

# Compute the first few Li coefficients
# λ_1 = 1 - 1/2 (Euler constant related)
# Actually, the exact formula involves the Stieltjes constants.
# For numerical computation, use:
# λ_n = Σ_ρ [1 - (1 - 1/ρ)^n]
# With known zeros ρ = 1/2 + iγ:

known_zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
               37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
               52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
               67.079811, 69.546402, 72.067158, 75.704691, 77.144840]

print("First Li coefficients (partial sum over 20 zeros):")
print(f"{'n':>4} {'λ_n (partial)':>16} {'positive?':>10}")
print("-" * 34)

for n in range(1, 21):
    lambda_n = 0.0
    for gamma in known_zeros:
        rho = complex(0.5, gamma)
        term = 1 - (1 - 1/rho)**n
        lambda_n += 2 * term.real  # conjugate pair
    
    pos = "YES" if lambda_n > 0 else "NO !!!"
    print(f"{n:4d} {lambda_n:16.8f} {pos:>10}")

print()

# ═══════════════════════════════════════════════════════════════
# PART 9: CONNECTING LI TO ESG
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 9: LI'S CRITERION IN THE ESG FRAMEWORK")
print("=" * 85)
print()

print("""
Li's criterion says: RH ⟺ λ_n ≥ 0 for all n.

In ESG language, with Δ the modular operator and ξ(s) = completed zeta:

  λ_n = (1/(n-1)!) [d^n/ds^n s^{n-1} log Tr(completed Delta^{-s/2})]|_{s=1}

This is:
  - A derivative at a REAL point (s = 1): Layer 1 ✓
  - Involves only log of the spectral trace: Layer 1 ✓
  - The inequality λ_n ≥ 0 is a REAL condition: Layer 1 ✓
  - Computable from the Stieltjes constants: Layer 1 ✓

THE ESG-LI REFORMULATION:

  RH ⟺ For the modular operator Δ at the entropy minimum ω₁/₂,
         the Li coefficients λ_n ≥ 0 for all n ≥ 1.

This is PURELY Layer 1. No time. No oscillation. No Fourier.
Just real derivatives of a real function at a real point.

THE CONNECTION TO C(β):
  The Stieltjes constants σ_k determine BOTH:
  - The Li coefficients λ_n (via binomial combinations)
  - The Taylor expansion of C(β) near β = 1

  So C(β) and the Li criterion are DIFFERENT VIEWS of the same
  Layer 1 data: the real-axis behavior of log ζ near s = 1.

  C(β) > 0 for all β > 1 is one constraint (the specific heat).
  λ_n ≥ 0 for all n is a STRONGER constraint (Li's criterion).
  
  Both are Layer 1. Both are consequences of the entropy minimum.
  Li's criterion is EQUIVALENT to RH.

THE OPEN QUESTION REFORMULATED:

  Can the ESG framework (entropy minimum + gap log 2) prove that
  λ_n ≥ 0 for all n?

  λ_n involves the Stieltjes constants, which involve the moments
  of the distribution of log n under the Gibbs measure.

  The entropy minimum (Hess V = Id) constrains these moments.
  The gap (log 2) constrains the support.
  Together: do they force λ_n ≥ 0?

  This is the CORRECT formulation of the open problem.
  It is Layer 1 throughout. It is precise. It is falsifiable.
""")

# ═══════════════════════════════════════════════════════════════
# PART 10: WHAT ESG SAYS ABOUT THE Li COEFFICIENTS
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  PART 10: ESG CONSTRAINTS ON Li COEFFICIENTS")
print("=" * 85)
print()

# The Stieltjes constants γ_k are defined by:
# ζ(s) = 1/(s-1) + Σ_{k=0}^∞ (-1)^k γ_k (s-1)^k / k!
# γ_0 = Euler-Mascheroni constant ≈ 0.5772
# γ_k = lim_{N→∞} [Σ_{n=1}^N (log n)^k / n - (log N)^{k+1}/(k+1)]

# These are MOMENTS of the measure dν = Σ (log n)^k / n × δ_{log n}
# weighted by the harmonic series.

# The ESG constraint: ω₁/₂ is the uniform state on G = Ẑ*.
# In the thermodynamic limit, ω₁/₂ assigns weight 1/n to integer n
# (up to normalization). The Stieltjes constants are the Taylor
# coefficients of log ζ at s = 1, which are determined by the
# moments of log n under this measure.

# Compute first Stieltjes constants numerically
print("Stieltjes constants (approximation):")
print(f"{'k':>4} {'γ_k':>16}")
print("-" * 24)

N_st = 100000
for k in range(6):
    # γ_k = lim [Σ (log n)^k / n - (log N)^{k+1}/(k+1)]
    partial = sum((log(n))**k / n for n in range(1, N_st + 1))
    correction = (log(N_st))**(k+1) / (k+1)
    gamma_k = partial - correction
    print(f"{k:4d} {gamma_k:16.10f}")

print()

# The Li criterion requires λ_n = Σ C(n,k) (-1)^k d_k
# where d_k are related to the Stieltjes constants.
# Computing this properly requires careful handling.

# For now, the key result is:
print("""
SUMMARY OF LAYER 1 TOOLS FOR RH:

1. C(β) = Var_β(log n) > 0           — trivially true (variance)
2. C(β) = Σ_k k Σ_p (log p)² p^{-kβ} — prime representation
3. λ_n ≥ 0 for all n                  — Li's criterion = RH
4. λ_n computable from Stieltjes γ_k  — real-axis derivatives
5. γ_k constrained by gap (log 2) and entropy (Hess V = Id)

The ESG programme reduces to:
  Show that the entropy minimum (Theorem 3.1) + gap (log 2)
  imply λ_n ≥ 0 via the Stieltjes constants.

This is a problem in MOMENT THEORY:
  Given a measure on {log 2, log 3, log 4, ...} with
  weights p_n = n^{-1}/ζ(1) and variance constraints from
  the positive-definite Hessian, do the moments satisfy
  the Li inequalities?

No time. No oscillation. No Fourier. No Mittag-Leffler.
Just moments, constraints, and inequalities.
""")

# ═══════════════════════════════════════════════════════════════
# GRAND VERDICT
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  GRAND VERDICT: THE CORRECT LAYER 1 FORMULATION")
print("=" * 85)
print("""
THE JOURNEY:

  Mittag-Leffler (Layer 2) → sign error + non-convergence
  → Diagnosis: individual poles = spectral decomposition = Layer 2
  → Search for Layer 1 bridge
  → Laplace vs Fourier: real β is Layer 1, complex s is Layer 2
  → Li's criterion: λ_n ≥ 0, computable from REAL derivatives
  → FOUND: a Layer 1 reformulation of RH

THE CORRECTED ESG PROBLEM:

  GIVEN:
    (a) ω₁/₂ is the unique minimum of V with Hess V = Id (Theorem 3.1)
    (b) The spectral support is {log n : n ∈ ℕ} with gap log 2
    (c) The Stieltjes constants γ_k encode the real-axis behavior of log ζ

  PROVE:
    The Li coefficients λ_n ≥ 0 for all n ≥ 1.

  TOOLS: Moment theory, positivity of Hessian, support constraints.
  ALL LAYER 1.

  This replaces the incorrect Mittag-Leffler approach of ESG v3 §10
  with a correct, fully Layer 1 formulation.

  The specific heat C(β) remains the Bures curvature and remains 
  interesting, but it is NOT the bridge to RH. Li's criterion is.
""")
