#!/usr/bin/env python3
"""
ESG I - Numerical Verification of Symmetry Relations
=====================================================
Truncated Bost-Connes system at beta = 1

We work on the GNS Hilbert space, which for a truncation to n = 1,...,N
is the space of Hilbert-Schmidt operators: dim = N^2.

Basis: |m><n| for m,n = 1,...,N  (vectorized as e_{(m-1)*N + (n-1)})

Key operators:
  K  = modular Hamiltonian (log Delta), eigenvalue log(n/m) on |m><n|
  J  = Tomita-Takesaki conjugation: J(|m><n|) = |n><m| (antilinear = transpose)
  gamma = arithmetic parity grading (candidate)
  T  = functional equation candidate (antilinear)
  C  = Galois conjugation candidate (antilinear)

We verify:
  (1) J^2 = Id, JKJ = -K                     [must hold - TT theorem]
  (2) gamma^2 = Id, {gamma, K} = 0 ?         [test - determines chiral structure]
  (3) T^2 = +1, C^2 = +1                     [test]
  (4) J = TC ?                                [test - central decomposition question]
  (5) KO-dimension signs                      [derive from results]

F. D. Blum, March 2026
"""

import numpy as np
from sympy import factorint
from collections import defaultdict
import sys

# ═══════════════════════════════════════════════════════════════
# SETUP
# ═══════════════════════════════════════════════════════════════

N = 30  # Truncation: n = 1, ..., N
DIM = N * N  # GNS Hilbert space dimension

print("=" * 70)
print(f"ESG I - NUMERICAL VERIFICATION")
print(f"Truncation: N = {N}, GNS dimension = {DIM}")
print("=" * 70)

# ─── Arithmetic functions ───────────────────────────────────────
def big_omega(n):
    """Number of prime factors with multiplicity: Omega(n)"""
    if n <= 1:
        return 0
    return sum(factorint(n).values())

def liouville_lambda(n):
    """Liouville function: (-1)^Omega(n)"""
    return (-1) ** big_omega(n)

# Precompute
omega_vals = [big_omega(n) for n in range(N + 1)]  # omega_vals[n] = Omega(n)
liouville = [liouville_lambda(n) for n in range(N + 1)]

print(f"\nArithmetic parity Omega(n) mod 2 for n = 1..{N}:")
for n in range(1, N + 1):
    parity = "+" if liouville[n] == 1 else "-"
    print(f"  n={n:2d}: Omega={omega_vals[n]}, parity={parity}", end="")
    if n % 5 == 0:
        print()
print()

# Count even/odd
n_even = sum(1 for n in range(1, N+1) if omega_vals[n] % 2 == 0)
n_odd = N - n_even
print(f"Parity distribution: {n_even} even, {n_odd} odd (ratio = {n_even/n_odd:.3f})")
print()

# ─── Index mapping ──────────────────────────────────────────────
def idx(m, n):
    """Map (m,n) pair (1-indexed) to GNS vector index (0-indexed)"""
    return (m - 1) * N + (n - 1)

def inv_idx(k):
    """Map GNS index back to (m,n) pair (1-indexed)"""
    m = k // N + 1
    n = k % N + 1
    return m, n

# ─── State weights ──────────────────────────────────────────────
# KMS_1 state: w_n = n^{-1} / Z_N  where Z_N = sum_{k=1}^N 1/k
Z_N = sum(1.0 / k for k in range(1, N + 1))
weights = np.array([1.0 / n / Z_N for n in range(1, N + 1)])

print(f"Partition function Z_N = H_{N} = {Z_N:.6f}")
print(f"Weights: w_1 = {weights[0]:.6f}, w_N = {weights[-1]:.8f}")
print()

# ═══════════════════════════════════════════════════════════════
# OPERATOR 1: MODULAR HAMILTONIAN K
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("OPERATOR 1: Modular Hamiltonian K = log(Delta)")
print("=" * 70)

# K is diagonal in the |m><n| basis with eigenvalue log(n/m)
K_diag = np.zeros(DIM)
for m in range(1, N + 1):
    for n in range(1, N + 1):
        K_diag[idx(m, n)] = np.log(n) - np.log(m)

K = np.diag(K_diag)

# Verify spectrum is symmetric around 0
K_eigs = np.sort(K_diag)
print(f"K spectrum: min = {K_eigs[0]:.4f}, max = {K_eigs[-1]:.4f}")
print(f"K spectrum symmetric? max + min = {K_eigs[-1] + K_eigs[0]:.10f}")

# Eigenvalue 0 has multiplicity N (the diagonal |n><n|)
n_zero = np.sum(np.abs(K_diag) < 1e-12)
print(f"Zero eigenvalue multiplicity: {n_zero} (expected: {N})")
print()

# ═══════════════════════════════════════════════════════════════
# OPERATOR 2: TOMITA-TAKESAKI CONJUGATION J
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("OPERATOR 2: Tomita-Takesaki conjugation J")
print("=" * 70)

# J acts as transpose: J(|m><n|) = |n><m|
# As a REAL matrix on the vectorized basis (since J is antilinear,
# we represent it as a real matrix acting on coefficients):
# J maps coefficient of |m><n| to coefficient of |n><m|
J_mat = np.zeros((DIM, DIM))
for m in range(1, N + 1):
    for n in range(1, N + 1):
        i = idx(m, n)
        j = idx(n, m)  # transposed
        J_mat[j, i] = 1.0

# Verify J^2 = Id
J_squared = J_mat @ J_mat
J2_error = np.max(np.abs(J_squared - np.eye(DIM)))
print(f"J^2 = Id? Max error: {J2_error:.2e}")

# Verify JKJ = -K
# Since J is antilinear, JKJ means: J * K * J (as real matrices, since
# K is diagonal real and J is a real permutation matrix)
JKJ = J_mat @ K @ J_mat
JKJ_plus_K = JKJ + K
JKJ_error = np.max(np.abs(JKJ_plus_K))
print(f"JKJ = -K? Max error: {JKJ_error:.2e}")

# Verify J preserves cyclic vector
# Omega = rho^{1/2} = sum_n w_n^{1/2} |n><n|
Omega = np.zeros(DIM)
for n in range(1, N + 1):
    Omega[idx(n, n)] = np.sqrt(weights[n - 1])

J_Omega = J_mat @ Omega
J_Omega_error = np.max(np.abs(J_Omega - Omega))
print(f"J Omega = Omega? Max error: {J_Omega_error:.2e}")
print()

# ═══════════════════════════════════════════════════════════════
# OPERATOR 3: ARITHMETIC GRADING gamma
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("OPERATOR 3: Arithmetic parity grading gamma")
print("=" * 70)
print("Testing THREE candidate gradings on the GNS space:")
print()

# --- Candidate A: Left action only ---
# gamma_L(|m><n|) = (-1)^Omega(m) |m><n|
print("--- Candidate A: gamma_L (left Liouville) ---")
gamma_L_diag = np.zeros(DIM)
for m in range(1, N + 1):
    for n in range(1, N + 1):
        gamma_L_diag[idx(m, n)] = liouville[m]
gamma_L = np.diag(gamma_L_diag)

# Check gamma_L^2 = Id
gL2_error = np.max(np.abs(gamma_L @ gamma_L - np.eye(DIM)))
print(f"  gamma_L^2 = Id? Max error: {gL2_error:.2e}")

# Check {gamma_L, K} = 0 (anticommutation)
anticomm_L = gamma_L @ K + K @ gamma_L
anticomm_L_max = np.max(np.abs(anticomm_L))
print(f"  {{gamma_L, K}} = 0? Max |entry|: {anticomm_L_max:.6f}")

# Check [gamma_L, K] (commutation)
comm_L = gamma_L @ K - K @ gamma_L
comm_L_max = np.max(np.abs(comm_L))
print(f"  [gamma_L, K] = 0? Max |entry|: {comm_L_max:.6f}")

# Check J gamma_L J vs gamma_L
JgLJ = J_mat @ gamma_L @ J_mat
JgLJ_vs_gL = np.max(np.abs(JgLJ - gamma_L))
JgLJ_vs_mgL = np.max(np.abs(JgLJ + gamma_L))
print(f"  J gamma_L J = +gamma_L? Max error: {JgLJ_vs_gL:.6f}")
print(f"  J gamma_L J = -gamma_L? Max error: {JgLJ_vs_mgL:.6f}")
print()

# --- Candidate B: Both sides (tensor product parity) ---
# gamma_LR(|m><n|) = (-1)^{Omega(m) + Omega(n)} |m><n|
print("--- Candidate B: gamma_LR (left x right Liouville) ---")
gamma_LR_diag = np.zeros(DIM)
for m in range(1, N + 1):
    for n in range(1, N + 1):
        gamma_LR_diag[idx(m, n)] = liouville[m] * liouville[n]
gamma_LR = np.diag(gamma_LR_diag)

gLR2_error = np.max(np.abs(gamma_LR @ gamma_LR - np.eye(DIM)))
print(f"  gamma_LR^2 = Id? Max error: {gLR2_error:.2e}")

anticomm_LR = gamma_LR @ K + K @ gamma_LR
anticomm_LR_max = np.max(np.abs(anticomm_LR))
print(f"  {{gamma_LR, K}} = 0? Max |entry|: {anticomm_LR_max:.6f}")

comm_LR = gamma_LR @ K - K @ gamma_LR
comm_LR_max = np.max(np.abs(comm_LR))
print(f"  [gamma_LR, K] = 0? Max |entry|: {comm_LR_max:.6f}")

JgLRJ = J_mat @ gamma_LR @ J_mat
JgLRJ_vs_gLR = np.max(np.abs(JgLRJ - gamma_LR))
print(f"  J gamma_LR J = +gamma_LR? Max error: {JgLRJ_vs_gLR:.6f}")
print()

# --- Candidate C: Right action only ---
# gamma_R(|m><n|) = (-1)^Omega(n) |m><n|
print("--- Candidate C: gamma_R (right Liouville) ---")
gamma_R_diag = np.zeros(DIM)
for m in range(1, N + 1):
    for n in range(1, N + 1):
        gamma_R_diag[idx(m, n)] = liouville[n]
gamma_R = np.diag(gamma_R_diag)

gR2_error = np.max(np.abs(gamma_R @ gamma_R - np.eye(DIM)))
print(f"  gamma_R^2 = Id? Max error: {gR2_error:.2e}")

anticomm_R = gamma_R @ K + K @ gamma_R
anticomm_R_max = np.max(np.abs(anticomm_R))
print(f"  {{gamma_R, K}} = 0? Max |entry|: {anticomm_R_max:.6f}")

comm_R = gamma_R @ K - K @ gamma_R
comm_R_max = np.max(np.abs(comm_R))
print(f"  [gamma_R, K] = 0? Max |entry|: {comm_R_max:.6f}")

JgRJ = J_mat @ gamma_R @ J_mat
JgRJ_vs_gR = np.max(np.abs(JgRJ - gamma_R))
JgRJ_vs_mgR = np.max(np.abs(JgRJ + gamma_R))
print(f"  J gamma_R J = +gamma_R? Max error: {JgRJ_vs_gR:.6f}")
print(f"  J gamma_R J = -gamma_R? Max error: {JgRJ_vs_mgR:.6f}")
print()

# ═══════════════════════════════════════════════════════════════
# DEEPER ANALYSIS: WHY ANTICOMMUTATION FAILS
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("ANALYSIS: Why does {gamma, K} = 0 fail?")
print("=" * 70)

# For gamma_LR (most natural candidate), {gamma_LR, K} = 0 requires:
# For each basis vector |m><n|:
# gamma_LR * K_eigenvalue + K_eigenvalue * gamma_LR = 0
# i.e., liouville[m]*liouville[n] * log(n/m) + log(n/m) * liouville[m]*liouville[n] = 0
# i.e., 2 * liouville[m]*liouville[n] * log(n/m) = 0
# This is zero only when m = n (log = 0) or when we consider off-diagonal
# Wait - gamma_LR is diagonal and K is diagonal in the SAME basis.
# Two diagonal operators ALWAYS commute: [gamma_LR, K] = 0.
# They can NEVER anticommute (unless one is zero).

print("\nCRITICAL OBSERVATION:")
print("K and all gamma candidates are DIAGONAL in the |m><n| basis.")
print("Two diagonal operators always COMMUTE: [gamma, K] = 0.")
print("They can NEVER ANTICOMMUTE: {gamma, K} = 0 is impossible")
print("(unless gamma or K is identically zero).")
print()
print("This means: NO grading defined by arithmetic parity of m, n")
print("can serve as a chiral symmetry that anticommutes with K.")
print()

# ═══════════════════════════════════════════════════════════════
# ALTERNATIVE: OFF-DIAGONAL GRADING
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("ALTERNATIVE: Off-diagonal grading (particle-hole structure)")
print("=" * 70)

# The natural Z_2 grading on the GNS space that DOES anticommute with K
# is the "particle-hole" grading that distinguishes |m><n| from |n><m|.
# Since K has eigenvalue log(n/m) on |m><n| and log(m/n) = -log(n/m) on |n><m|,
# a grading that swaps these would anticommute with K.
# But that's exactly what J does! And J is antilinear, not linear.
#
# For a LINEAR grading, consider:
# sigma(|m><n|) = sign(log(n/m)) * |m><n|  for m != n
# sigma(|n><n|) = 0 (or undefined)
#
# But this is not well-defined on the zero-eigenvalue subspace.
#
# The proper way to think about this: the modular Hamiltonian K
# naturally decomposes the GNS space as H = H_+ + H_0 + H_-
# where H_+ = {K > 0}, H_0 = {K = 0}, H_- = {K < 0}.
# J maps H_+ <-> H_-.
# The zero subspace H_0 = span{|n><n|} is the diagonal (classical) part.

print("Natural decomposition of GNS space by sign of K eigenvalue:")
n_pos = np.sum(K_diag > 1e-12)
n_neg = np.sum(K_diag < -1e-12)
n_zero_K = np.sum(np.abs(K_diag) < 1e-12)
print(f"  H_+ (K > 0): dim = {n_pos}")
print(f"  H_0 (K = 0): dim = {n_zero_K}")
print(f"  H_- (K < 0): dim = {n_neg}")
print(f"  H_+ and H_- related by J (transpose)")
print(f"  Expected: dim H_+ = dim H_- = N(N-1)/2 = {N*(N-1)//2}")
print()

# ═══════════════════════════════════════════════════════════════
# THE SIGN OPERATOR
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("THE SIGN OPERATOR: Sigma = sign(K)")
print("=" * 70)

# Define Sigma = sign(K) on H_+ and H_-, zero on H_0
Sigma_diag = np.zeros(DIM)
for i in range(DIM):
    if K_diag[i] > 1e-12:
        Sigma_diag[i] = 1.0
    elif K_diag[i] < -1e-12:
        Sigma_diag[i] = -1.0
    else:
        Sigma_diag[i] = 0.0
Sigma = np.diag(Sigma_diag)

# On H_+ + H_- (excluding H_0), Sigma^2 = Id
# On full space, Sigma^2 = Projection onto H_+ + H_-
print("Sigma^2 on H_+ + H_-:")
Sigma2 = Sigma @ Sigma
proj_offdiag = np.diag(np.array([0.0 if abs(K_diag[i]) < 1e-12 else 1.0 for i in range(DIM)]))
S2_error = np.max(np.abs(Sigma2 - proj_offdiag))
print(f"  Sigma^2 = Proj(H_+ + H_-)? Max error: {S2_error:.2e}")

# Check {Sigma, K} on the off-diagonal subspace
# Sigma is diagonal with +1/-1, K is diagonal with log(n/m)
# {Sigma, K} is diagonal with entries sign(log(n/m)) * log(n/m) + log(n/m) * sign(log(n/m))
# = 2 |log(n/m)| != 0 in general
# So {Sigma, K} != 0 either! Sigma COMMUTES with K (both diagonal).

print("\nSigma commutes with K (both diagonal): [Sigma, K] = 0")
print("Sigma does NOT anticommute with K: {Sigma, K} != 0")
print()

# ═══════════════════════════════════════════════════════════════
# FUNDAMENTAL OBSTRUCTION
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("FUNDAMENTAL OBSTRUCTION TO CHIRAL GRADING")
print("=" * 70)

print("""
THEOREM (Numerical): No diagonal operator on the GNS space H_1 can
serve as a chiral grading for the modular Hamiltonian K.

Proof: K is diagonal in the |m><n| basis. Any diagonal operator gamma
commutes with K: [gamma, K] = 0. The anticommutation relation
{gamma, K} = 0 would require gamma*K + K*gamma = 2*gamma*K = 0,
which forces either gamma = 0 or K = 0 on every basis vector.
Since K != 0 on the off-diagonal part, gamma must vanish there,
leaving at most a grading on the N-dimensional diagonal subspace
H_0 = ker(K), which is useless for classification purposes. QED

CONSEQUENCE: Any chiral symmetry for the BC system must be
OFF-DIAGONAL in the |m><n| basis -- it must MIX basis vectors,
not merely assign signs to them.

The only natural off-diagonal operator we have is J itself
(the Tomita-Takesaki conjugation), but J is ANTILINEAR.
A linear chiral operator must be found elsewhere.
""")

# ═══════════════════════════════════════════════════════════════
# SEARCH FOR OFF-DIAGONAL CHIRAL OPERATOR
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("SEARCH: Off-diagonal chiral operator using Liouville function")
print("=" * 70)

# Idea: define S(|m><n|) = liouville[m] * liouville[n] * |n><m|
# This is LINEAR, squares to gamma_LR (since applying twice gives
# liouville[m]*liouville[n]*liouville[n]*liouville[m] * |m><n| = |m><n|),
# and maps K eigenvalue log(n/m) to eigenvalue log(m/n) = -log(n/m).

print("Candidate: S_L(|m><n|) = lambda(m)*lambda(n) * |n><m|")
print("where lambda(n) = (-1)^Omega(n) is the Liouville function.")
print()

S_L = np.zeros((DIM, DIM))
for m in range(1, N + 1):
    for n in range(1, N + 1):
        i_in = idx(m, n)
        i_out = idx(n, m)
        S_L[i_out, i_in] = liouville[m] * liouville[n]

# Check S_L^2
S_L2 = S_L @ S_L
S_L2_vs_Id = np.max(np.abs(S_L2 - np.eye(DIM)))
print(f"S_L^2 = Id? Max error: {S_L2_vs_Id:.2e}")

# Check S_L K + K S_L = 0 (anticommutation with K)
anticomm_SL_K = S_L @ K + K @ S_L
anticomm_SL_K_max = np.max(np.abs(anticomm_SL_K))
print(f"{{S_L, K}} = 0? Max |entry|: {anticomm_SL_K_max:.2e}")

# Check S_L K - K S_L (commutation)
comm_SL_K = S_L @ K - K @ S_L
comm_SL_K_max = np.max(np.abs(comm_SL_K))
print(f"[S_L, K] = 0? Max |entry|: {comm_SL_K_max:.6f}")

# Check J S_L J vs S_L
JS_LJ = J_mat @ S_L @ J_mat
JS_LJ_vs_SL = np.max(np.abs(JS_LJ - S_L))
JS_LJ_vs_mSL = np.max(np.abs(JS_LJ + S_L))
print(f"J S_L J = +S_L? Max error: {JS_LJ_vs_SL:.6f}")
print(f"J S_L J = -S_L? Max error: {JS_LJ_vs_mSL:.6f}")
print()

# Also try without Liouville: plain transpose
print("Candidate: S_plain(|m><n|) = |n><m| (= J, but treated as linear)")
S_plain = np.zeros((DIM, DIM))
for m in range(1, N + 1):
    for n in range(1, N + 1):
        S_plain[idx(n, m), idx(m, n)] = 1.0

# Note: S_plain = J_mat as a matrix, but J is antilinear while S_plain is linear
Sp2 = S_plain @ S_plain
Sp2_vs_Id = np.max(np.abs(Sp2 - np.eye(DIM)))
print(f"S_plain^2 = Id? Max error: {Sp2_vs_Id:.2e}")

anticomm_Sp_K = S_plain @ K + K @ S_plain
anticomm_Sp_K_max = np.max(np.abs(anticomm_Sp_K))
print(f"{{S_plain, K}} = 0? Max |entry|: {anticomm_Sp_K_max:.2e}")
print()

# ═══════════════════════════════════════════════════════════════
# ANALYSIS OF RESULTS
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("ANALYSIS OF S_L AND S_plain")
print("=" * 70)

print("""
Both S_L (with Liouville signs) and S_plain (plain transpose) satisfy:
  - S^2 = Id  (involution)
  - {S, K} = 0  (anticommutation with modular Hamiltonian)

The key difference:
  S_plain = J_mat (same matrix), so "linear J" already anticommutes.
  S_L adds arithmetic signs via Liouville function.

For S_plain (= matrix form of J):
  This is the transpose operator, which maps |m><n| -> |n><m|.
  Treated as LINEAR, it anticommutes with K because it swaps
  the eigenvalue log(n/m) to log(m/n) = -log(n/m).
  
  BUT: in quantum mechanics, the physical J is ANTILINEAR.
  The LINEAR version of "transpose" is the SUPERTRANSPOSE, which
  does not have the same physical interpretation.

For S_L (Liouville-weighted transpose):
  This also anticommutes with K and also satisfies S^2 = +1.
  The Liouville signs are "transparent" to the anticommutation
  because they multiply both sides equally.
""")

# ═══════════════════════════════════════════════════════════════
# DECOMPOSITION TEST: J = T * C ?
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("DECOMPOSITION TEST: J = T * C")
print("=" * 70)

# T candidate: complex conjugation composed with spectral reflection
# On the GNS space with REAL coefficients, complex conjugation is trivial.
# We need to work with COMPLEX coefficients to make this meaningful.
#
# T (functional equation): antilinear, maps spectral parameter t -> -t
# On the GNS space: T(|m><n|) = something involving n^{-1}, m^{-1}
# 
# C (Galois sigma_{-1}): antilinear, complex conjugation on cyclotomic values
# On the GNS space with integer labels: C is essentially complex conjugation
# of coefficients, since |n> has no complex phase in this basis.

# Let's define them concretely:
# C = plain complex conjugation (antilinear): in the real basis, C = Id
# T = J composed with C^{-1} = J (if C = Id)

# This is the crux: if C is trivial on the real basis, then T = J.
# There's no non-trivial decomposition J = TC in this representation.

print("In the standard basis |m><n| with REAL coefficients:")
print("  C (complex conjugation) acts trivially: C = Id")
print("  Therefore J = TC forces T = J")
print("  The decomposition J = TC is trivially satisfied but UNINFORMATIVE")
print()
print("To make the decomposition non-trivial, we need a basis where")
print("the cyclotomic structure of KMS_1 states is visible.")
print("This requires the Fourier/Peter-Weyl decomposition on G = Z-hat*.")
print()

# ═══════════════════════════════════════════════════════════════
# PETER-WEYL: Dirichlet characters
# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("FOURIER ANALYSIS: Dirichlet characters on (Z/qZ)*")
print("=" * 70)

# For a finite approximation, we work with G_q = (Z/qZ)* for q = product of first few primes
# G_q approximates Z-hat* = lim (Z/qZ)*
#
# The Dirichlet characters mod q form the dual group of (Z/qZ)*
# In this basis, the Galois action sigma_{-1} sends chi -> chi-bar
# This is non-trivial for complex characters!

from math import gcd

def euler_phi(n):
    return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)

def get_group_elements(q):
    """Elements of (Z/qZ)*"""
    return [k for k in range(1, q+1) if gcd(k, q) == 1]

def compute_characters(q):
    """Compute all Dirichlet characters mod q via DFT on (Z/qZ)*"""
    elems = get_group_elements(q)
    phi = len(elems)
    
    # Build the group table: elem_to_idx
    elem_to_idx = {e: i for i, e in enumerate(elems)}
    
    # The group (Z/qZ)* is abelian, so we find generators and compute chars
    # For simplicity, use the DFT approach: characters are the eigenvectors
    # of the regular representation
    
    # Regular representation matrix for generator g:
    # For each a in G, left multiplication by g sends a -> ga mod q
    # We need a generator. For prime q, any primitive root works.
    
    # General approach: compute the character table directly
    # chi_j(g^k) = omega^{jk} where omega = exp(2pi i / ord(g))
    
    # Even more general: decompose (Z/qZ)* as product of cyclic groups
    # and use tensor product of 1D DFTs
    
    # SIMPLEST: just find characters by eigenvalues of multiplication operators
    # For each a in (Z/qZ)*, define M_a: f(x) -> f(ax mod q)
    # Characters are simultaneous eigenvectors of all M_a
    
    # Build multiplication matrix for the smallest generator
    M = np.zeros((phi, phi), dtype=complex)
    for i, a in enumerate(elems):
        for j, b in enumerate(elems):
            ab = (a * b) % q
            k = elem_to_idx[ab]
            M[k, j] += 1  # but this mixes things...
    
    # Actually, just use the regular representation of a generator
    # Find a generator (or set of generators)
    # For prime q, we can find a primitive root
    
    # Let's just compute characters for small q directly
    chars = np.zeros((phi, phi), dtype=complex)
    
    # Use eigendecomposition of the shift operator
    # Pick the first non-trivial element as "generator"
    if len(elems) < 2:
        return np.ones((1, 1), dtype=complex), elems
    
    g = elems[1]  # first element > 1
    
    # Build shift matrix: S_{ij} = 1 if elems[j]*g = elems[i] mod q
    S = np.zeros((phi, phi))
    for j, b in enumerate(elems):
        gb = (g * b) % q
        i = elem_to_idx[gb]
        S[i, j] = 1.0
    
    # Eigenvalues of S are the character values chi(g)
    evals, evecs = np.linalg.eig(S)
    
    # Each eigenvector gives a character
    # chi_k(elems[j]) proportional to evecs[j, k]
    # Normalize so chi_k(1) = 1
    idx_1 = elem_to_idx[1 % q if 1 % q in elem_to_idx else 1]
    for k in range(phi):
        if abs(evecs[idx_1, k]) > 1e-10:
            evecs[:, k] /= evecs[idx_1, k]
    
    return evecs, elems

# Test with q = 5 (phi = 4, so 4 characters)
q = 5
chars, elems = compute_characters(q)
phi_q = euler_phi(q)
print(f"Working with q = {q}, phi({q}) = {phi_q}")
print(f"Group elements (Z/{q}Z)*: {elems}")
print()

# Display character table
print(f"Character table (Z/{q}Z)*:")
print(f"{'chi\\g':>8}", end="")
for g in elems:
    print(f"  {g:>8}", end="")
print()
for k in range(phi_q):
    print(f"  chi_{k}:", end="")
    for j in range(phi_q):
        v = chars[j, k]
        if abs(v.imag) < 1e-8:
            print(f"  {v.real:>8.4f}", end="")
        else:
            print(f"  {v.real:+.3f}{v.imag:+.3f}i", end="")
    print()
print()

# Now check: sigma_{-1} sends g -> g^{-1} mod q
# In terms of characters: sigma_{-1}(chi)(g) = chi(g^{-1}) = chi-bar(g) for unitary chars
print("Action of sigma_{-1} (g -> g^{-1} mod q):")
for k in range(phi_q):
    # Compute chi_k(g^{-1}) for each g
    chi_k_conj = np.zeros(phi_q, dtype=complex)
    for j, g in enumerate(elems):
        # Find g^{-1} mod q
        g_inv = pow(g, -1, q)
        j_inv = elems.index(g_inv)
        chi_k_conj[j] = chars[j_inv, k]
    
    # Find which chi_l this matches
    matched = -1
    for l in range(phi_q):
        if np.max(np.abs(chi_k_conj - chars[:, l])) < 1e-6:
            matched = l
            break
    
    is_real = (matched == k)
    is_complex = (matched != k and matched >= 0)
    print(f"  sigma_{{-1}}(chi_{k}) = chi_{matched}"
          f"  {'(real character)' if is_real else '(complex pair)' if is_complex else '(??)'}")

print()
print("OBSERVATION: sigma_{-1} acts non-trivially on complex characters,")
print("pairing chi with chi-bar. On the character basis, C is NOT trivial.")
print("The decomposition J = TC becomes non-trivial in this basis.")

# ═══════════════════════════════════════════════════════════════
# LARGER q: q = 12 (more structure)
# ═══════════════════════════════════════════════════════════════
print()
print("=" * 70)
print(f"EXTENDED: q = 12, phi(12) = {euler_phi(12)}")
print("=" * 70)

q2 = 12
chars2, elems2 = compute_characters(q2)
phi_q2 = euler_phi(q2)
print(f"Group elements (Z/{q2}Z)*: {elems2}")

print(f"\nAction of sigma_{{-1}} on characters mod {q2}:")
n_real_chars = 0
n_complex_pairs = 0
for k in range(phi_q2):
    chi_k_conj = np.zeros(phi_q2, dtype=complex)
    for j, g in enumerate(elems2):
        g_inv = pow(g, -1, q2)
        j_inv = elems2.index(g_inv)
        chi_k_conj[j] = chars2[j_inv, k]
    
    matched = -1
    for l in range(phi_q2):
        if np.max(np.abs(chi_k_conj - chars2[:, l])) < 1e-6:
            matched = l
            break
    
    is_real = (matched == k)
    if is_real:
        n_real_chars += 1
    else:
        n_complex_pairs += 1
    
    print(f"  sigma_{{-1}}(chi_{k}) = chi_{matched}"
          f"  {'(REAL)' if is_real else '(COMPLEX pair)'}")

print(f"\nReal characters: {n_real_chars}")
print(f"Characters in complex pairs: {n_complex_pairs}")
print(f"Number of complex pairs: {n_complex_pairs // 2}")

# ═══════════════════════════════════════════════════════════════
# SUMMARY OF FINDINGS
# ═══════════════════════════════════════════════════════════════
print()
print("=" * 70)
print("SUMMARY OF NUMERICAL FINDINGS")
print("=" * 70)

print("""
1. TOMITA-TAKESAKI (CONFIRMED):
   J^2 = +1, JKJ = -K, J*Omega = Omega
   All verified to machine precision. This is theorem, not surprise.

2. DIAGONAL GRADING (RULED OUT):
   No diagonal operator on the GNS space (including arithmetic parity
   via Liouville function) can anticommute with the modular Hamiltonian K.
   This is because K is itself diagonal in the |m><n| basis.
   CONSEQUENCE: The naive BDI grading proposed in the original BDI
   document DOES NOT WORK.

3. OFF-DIAGONAL OPERATORS (FOUND):
   The transpose operator S(|m><n|) = |n><m| (and its Liouville-weighted
   variant) DOES anticommute with K. It satisfies S^2 = +1 and {S,K} = 0.
   However, S is the MATRIX FORM of J (the antilinear Tomita conjugation).
   Treating it as a linear operator is the "supertranspose" and has
   different physical meaning.

4. DECOMPOSITION J = TC (REQUIRES CHARACTER BASIS):
   In the integer basis |n>, the Galois conjugation C is trivial
   (because all weights are real). The decomposition J = TC becomes
   non-trivial ONLY in the Fourier basis (Dirichlet characters),
   where sigma_{-1} pairs complex characters chi <-> chi-bar.
   
5. KO-DIMENSION (PARTIALLY DETERMINED):
   eps = +1 (from J^2 = +1)  -- PROVED
   eps' = -1 (from JKJ = -K) -- PROVED
   eps'' (from J*gamma*J vs gamma) -- UNDETERMINED because no valid
   linear chiral grading gamma has been found.
   
   WITHOUT eps'', the KO-dimension is ambiguous between:
   - KO-dim 1 (BDI) if eps'' exists
   - KO-dim 6 (C) if no grading exists
   
6. CRITICAL OPEN QUESTION:
   The existence of a chiral grading requires a LINEAR operator that
   anticommutes with K and is not simply the matrix form of J.
   This operator, if it exists, must come from a structure BEYOND
   the Hecke algebra -- possibly from the interaction between the
   Galois action (which becomes non-trivial in the character basis)
   and the spectral flow.
   
   The Peter-Weyl decomposition of L^2(Z-hat*) into Dirichlet
   characters is the natural framework for this search.
""")

print("=" * 70)
print("VERDICT")
print("=" * 70)
print("""
The original BDI classification document was WRONG in its identification
of the chiral symmetry (Proposition 3.6). The numerical verification
confirms: no arithmetic parity grading anticommutes with the modular
Hamiltonian in the integer basis.

However, the partial KO-dimension result (eps = +1, eps' = -1) stands
on solid ground (Tomita-Takesaki). The unresolved question is whether
the Bost-Connes system admits a chiral structure (BDI, KO-dim 1) or
not (class C, KO-dim 6, or possibly no Z invariant at all).

NEXT STEP: Reformulate the entire calculation in the Dirichlet
character basis, where the Galois action sigma_{-1} is non-trivial
and may provide the missing chiral structure.
""")
