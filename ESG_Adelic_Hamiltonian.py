#!/usr/bin/env python3
"""
ESG I - THE ADELIC HAMILTONIAN
==============================
K_ad = (K - CKC) / 2  (antisymmetric projection under Galois)

This operator satisfies CK_adC = -K_ad BY CONSTRUCTION.
We test whether it also satisfies JK_adJ = -K_ad (compatibility with TT)
and whether it has a spectral gap.

We also test alternative constructions:
  K_p  = p-adic valuation Hamiltonian
  K_sym = (K + CKC) / 2  (symmetric projection, for comparison)

F. D. Blum, March 2026
"""

import numpy as np
from math import gcd, log
from sympy import factorint
import time

np.set_printoptions(precision=8, suppress=True, linewidth=130)

def liouville(n):
    if n <= 1:
        return 1
    return (-1) ** sum(factorint(n).values())

def mod_inverse(a, q):
    for b in range(1, q):
        if (a * b) % q == 1:
            return b
    return None

def p_adic_val(n, p):
    """v_p(n): largest k such that p^k divides n"""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

# ═══════════════════════════════════════════════════════════════
# MAIN COMPUTATION
# ═══════════════════════════════════════════════════════════════

for q in [5, 7, 13, 30, 31, 60]:
    G = [a for a in range(1, q) if gcd(a, q) == 1]
    phi_q = len(G)
    DIM = phi_q * phi_q
    g_to_idx = {g: i for i, g in enumerate(G)}
    
    if DIM > 10000:
        continue
    
    print("=" * 85)
    print(f"  q = {q},  φ({q}) = {phi_q},  GNS dim = {DIM}")
    print("=" * 85)
    
    # ─── Build operators ────────────────────────────────────────
    
    # K = log(n/m) diagonal
    K_diag = np.zeros(DIM)
    for mi in range(phi_q):
        for ni in range(phi_q):
            K_diag[mi * phi_q + ni] = log(G[ni]) - log(G[mi])
    K = np.diag(K_diag)
    
    # J = transpose permutation
    J_perm = np.zeros(DIM, dtype=int)
    for mi in range(phi_q):
        for ni in range(phi_q):
            J_perm[mi * phi_q + ni] = ni * phi_q + mi
    J_mat = np.zeros((DIM, DIM))
    for i in range(DIM):
        J_mat[J_perm[i], i] = 1.0
    
    # C = σ_{-1} permutation: (m, n) → (m⁻¹ mod q, n⁻¹ mod q)
    C_perm = np.zeros(DIM, dtype=int)
    for mi in range(phi_q):
        for ni in range(phi_q):
            m_inv = mod_inverse(G[mi], q)
            n_inv = mod_inverse(G[ni], q)
            mi_inv = g_to_idx[m_inv]
            ni_inv = g_to_idx[n_inv]
            C_perm[mi * phi_q + ni] = mi_inv * phi_q + ni_inv
    C_mat = np.zeros((DIM, DIM))
    for i in range(DIM):
        C_mat[C_perm[i], i] = 1.0
    
    # CKC diagonal (permute K eigenvalues by C)
    CKC_diag = K_diag[C_perm]
    CKC = np.diag(CKC_diag)
    
    # ═══════════════════════════════════════════════════════════
    # CONSTRUCTION 1: K_ad = (K - CKC) / 2
    # ═══════════════════════════════════════════════════════════
    
    print("\n--- CONSTRUCTION 1: K_ad = (K - CKC) / 2 ---")
    
    K_ad_diag = (K_diag - CKC_diag) / 2.0
    K_ad = np.diag(K_ad_diag)
    
    # Test CK_adC = -K_ad (should be exact by construction)
    CKadC_diag = K_ad_diag[C_perm]
    err_CKadC = np.max(np.abs(CKadC_diag + K_ad_diag))
    print(f"  CK_adC = -K_ad? err = {err_CKadC:.2e}  (exact by construction)")
    
    # Test JK_adJ = -K_ad  (THIS IS THE KEY TEST)
    JKadJ_diag = K_ad_diag[J_perm]
    err_JKadJ_minus = np.max(np.abs(JKadJ_diag + K_ad_diag))
    err_JKadJ_plus = np.max(np.abs(JKadJ_diag - K_ad_diag))
    print(f"  JK_adJ = -K_ad? err = {err_JKadJ_minus:.2e}")
    print(f"  JK_adJ = +K_ad? err = {err_JKadJ_plus:.2e}")
    
    # Spectrum of K_ad
    K_ad_eigs = np.sort(K_ad_diag)
    K_ad_unique = np.sort(np.unique(np.round(K_ad_diag, 10)))
    n_zero_ad = np.sum(np.abs(K_ad_diag) < 1e-10)
    print(f"  Spectrum: min = {K_ad_eigs[0]:.6f}, max = {K_ad_eigs[-1]:.6f}")
    print(f"  Zero eigenvalues: {n_zero_ad} / {DIM}")
    print(f"  Unique eigenvalues: {len(K_ad_unique)}")
    
    # Gap: smallest non-zero |eigenvalue|
    nonzero_eigs = np.abs(K_ad_diag[np.abs(K_ad_diag) > 1e-10])
    if len(nonzero_eigs) > 0:
        gap = np.min(nonzero_eigs)
        print(f"  Spectral gap: Δ = {gap:.8f}")
    else:
        print(f"  WARNING: All eigenvalues are zero!")
    
    # Display sample eigenvalues
    print(f"  Sample (m, n, K, CKC, K_ad):")
    count = 0
    for mi in range(min(phi_q, 5)):
        for ni in range(min(phi_q, 5)):
            if mi == ni:
                continue
            m, n = G[mi], G[ni]
            idx = mi * phi_q + ni
            m_inv = mod_inverse(m, q)
            n_inv = mod_inverse(n, q)
            print(f"    ({m:3d},{n:3d}): K={K_diag[idx]:+.6f}  "
                  f"CKC=log({n_inv}/{m_inv})={CKC_diag[idx]:+.6f}  "
                  f"K_ad={K_ad_diag[idx]:+.6f}")
            count += 1
            if count >= 8:
                break
        if count >= 8:
            break
    
    # ═══════════════════════════════════════════════════════════
    # What IS K_ad arithmetically?
    # K_ad(m,n) = (log(n/m) - log(n⁻¹/m⁻¹)) / 2
    #           = (log(n/m) - log(m·m⁻¹/(n·n⁻¹)) - log(m/n)) / 2
    # Wait, let's compute directly:
    # K_ad = (log(n/m) - log(n_inv/m_inv)) / 2
    #       = (1/2) log(n·m_inv / (m·n_inv))
    # ═══════════════════════════════════════════════════════════
    
    print(f"\n  Arithmetic interpretation of K_ad:")
    print(f"  K_ad(m,n) = (1/2) log[ (n · m_inv) / (m · n_inv) ]")
    print(f"  where m_inv, n_inv are modular inverses mod {q}")
    print(f"  Equivalently: K_ad(m,n) = (1/2) log[ (n/m) · (m_inv/n_inv) ]")
    print(f"  = (1/2) [log(n/m) + log(m_inv/n_inv)]")
    print(f"  = AVERAGE of archimédean and 'modular' log-ratios")
    print()
    
    # ═══════════════════════════════════════════════════════════
    # CONSTRUCTION 2: K_sym = (K + CKC) / 2  (the defect)
    # ═══════════════════════════════════════════════════════════
    
    print("--- CONSTRUCTION 2: K_sym = (K + CKC) / 2  (symmetric part) ---")
    
    K_sym_diag = (K_diag + CKC_diag) / 2.0
    
    # K = K_ad + K_sym (decomposition)
    decomp_err = np.max(np.abs(K_diag - (K_ad_diag + K_sym_diag)))
    print(f"  K = K_ad + K_sym? err = {decomp_err:.2e}")
    
    # K_sym commutes with C (by construction)
    CKsymC_diag = K_sym_diag[C_perm]
    err_CKsymC = np.max(np.abs(CKsymC_diag - K_sym_diag))
    print(f"  CK_symC = +K_sym? err = {err_CKsymC:.2e}")
    
    # K_sym anticommutes with J?
    JKsymJ_diag = K_sym_diag[J_perm]
    err_JKsymJ_minus = np.max(np.abs(JKsymJ_diag + K_sym_diag))
    err_JKsymJ_plus = np.max(np.abs(JKsymJ_diag - K_sym_diag))
    print(f"  JK_symJ = -K_sym? err = {err_JKsymJ_minus:.2e}")
    print(f"  JK_symJ = +K_sym? err = {err_JKsymJ_plus:.2e}")
    
    # K_sym spectrum
    K_sym_nonzero = np.abs(K_sym_diag[np.abs(K_sym_diag) > 1e-10])
    print(f"  K_sym zero eigenvalues: {np.sum(np.abs(K_sym_diag) < 1e-10)} / {DIM}")
    if len(K_sym_nonzero) > 0:
        print(f"  K_sym gap: {np.min(K_sym_nonzero):.8f}")
    print()
    
    # ═══════════════════════════════════════════════════════════
    # S_L vs K_ad
    # ═══════════════════════════════════════════════════════════
    
    print("--- S_L (Liouville) vs K_ad ---")
    
    lam = {g: liouville(g) for g in G}
    
    # S_L anticommutes with K_ad?
    # S_L|m><n| = λ(m)λ(n)|n><m|
    # S_L K_ad|m><n| = K_ad(m,n) λ(m)λ(n) |n><m|
    # K_ad S_L|m><n| = λ(m)λ(n) K_ad(n,m) |n><m|
    # anticomm: K_ad(m,n) + K_ad(n,m) = 0? 
    # K_ad(m,n) = (log(n/m) - log(n_inv/m_inv))/2
    # K_ad(n,m) = (log(m/n) - log(m_inv/n_inv))/2 = -(log(n/m) - log(n_inv/m_inv))/2 = -K_ad(m,n)
    # YES! K_ad is antisymmetric in (m,n), so {S_L, K_ad} = 0 automatically.
    
    anticomm_SL_Kad = 0.0
    for mi in range(phi_q):
        for ni in range(phi_q):
            idx_mn = mi * phi_q + ni
            idx_nm = ni * phi_q + mi
            val = K_ad_diag[idx_mn] + K_ad_diag[idx_nm]
            anticomm_SL_Kad = max(anticomm_SL_Kad, abs(val))
    
    print(f"  {{S_L, K_ad}} = 0? err = {anticomm_SL_Kad:.2e}")
    print(f"  (K_ad is antisymmetric in (m,n), so this is automatic)")
    
    # S_L commutes with K_sym?
    comm_SL_Ksym = 0.0
    for mi in range(phi_q):
        for ni in range(phi_q):
            idx_mn = mi * phi_q + ni
            idx_nm = ni * phi_q + mi
            val = K_sym_diag[idx_mn] - K_sym_diag[idx_nm]
            comm_SL_Ksym = max(comm_SL_Ksym, abs(val))
    
    print(f"  [S_L, K_sym] = 0? err = {comm_SL_Ksym:.2e}")
    print(f"  (K_sym is symmetric in (m,n), so this is automatic)")
    print()
    
    # ═══════════════════════════════════════════════════════════
    # FULL SYMMETRY TABLE FOR K_ad
    # ═══════════════════════════════════════════════════════════
    
    print("--- FULL SYMMETRY TABLE FOR K_ad ---")
    print()
    
    # J properties
    print(f"  J:  J² = +1 (exact)")
    print(f"      JK_adJ = {'−K_ad' if err_JKadJ_minus < 1e-8 else '+K_ad' if err_JKadJ_plus < 1e-8 else '???'}"
          f"  (err = {min(err_JKadJ_minus, err_JKadJ_plus):.2e})")
    
    # C properties
    print(f"  C:  C² = +1 (exact)")
    print(f"      CK_adC = −K_ad  (exact by construction)")
    
    # S_L properties
    print(f"  S_L: S_L² = +1 (exact)")
    print(f"       {{S_L, K_ad}} = 0  (exact, K_ad antisymmetric)")
    print(f"       JS_LJ = +S_L  (exact)")
    
    # Classification attempt
    print()
    if err_JKadJ_minus < 1e-8:
        print(f"  J: ε=+1, ε'=−1  (anticommutes with K_ad)")
        print(f"  C: ε=+1, ε'=−1  (anticommutes with K_ad)")
        print(f"  Both J AND C anticommute with K_ad!")
        print(f"  S = JC: linear, S² = ?")
        
        # Compute JC
        JC_perm = J_perm[C_perm]
        JC2 = JC_perm[JC_perm]
        jc2_is_id = np.all(JC2 == np.arange(DIM))
        print(f"  (JC)² = Id? {jc2_is_id}")
        
        # Does JC anticommute with K_ad?
        JCKadJC_diag = K_ad_diag[JC_perm]
        # JC is linear, so JCKadJC means permute twice
        JCKadJC_diag2 = JCKadJC_diag  # already permuted once; need to permute by (JC)⁻¹
        # For permutation P: PKP⁻¹ has diagonal P·diag(K)
        # (JC)K_ad(JC)⁻¹ diagonal = K_ad[JC_perm_inverse]... 
        # Actually for diagonal K_ad: (JC)K_ad(JC)⁻¹ |i⟩ = K_ad[(JC)⁻¹(i)] |i⟩
        JC_inv_perm = np.zeros(DIM, dtype=int)
        for i in range(DIM):
            JC_inv_perm[JC_perm[i]] = i
        
        JCKadJCinv_diag = K_ad_diag[JC_inv_perm]
        err_JCKad_anti = np.max(np.abs(JCKadJCinv_diag + K_ad_diag))
        err_JCKad_comm = np.max(np.abs(JCKadJCinv_diag - K_ad_diag))
        print(f"  (JC)K_ad(JC)⁻¹ = -K_ad? err = {err_JCKad_anti:.2e}")
        print(f"  (JC)K_ad(JC)⁻¹ = +K_ad? err = {err_JCKad_comm:.2e}")
        
        if err_JCKad_anti < 1e-8:
            print(f"  → S = JC anticommutes with K_ad (chiral!)")
            print(f"  → T=J (ε=+1, ε'=−1), C (ε=+1, ε'=−1), S=JC (chiral)")
            print(f"  → THIS IS CLASS BDI with ε=+1, ε'=−1")
            
            # Final check: ε'' = sign of J·S·J·S
            # S = JC, so JSJ = J(JC)J = J²CJ = CJ (since J²=1)
            # Then J·S = J·JC = C, and S·J = JC·J = JCJ
            # Need to compute in detail
            # J(JC)J⁻¹ = J(JC)J (since J²=1) = (JJ)(CJ) = CJ
            # So JS = J(JC) = C, SJ = (JC)J = J(CJ) ... hmm need to be careful
            # Actually S = JC as permutation. JSJ⁻¹ = J(JC)J = (JJ)(CJ) = CJ
            # But CJ vs JC: do they differ?
            CJ_perm = C_perm[J_perm]
            JC_vs_CJ = np.max(np.abs(JC_perm - CJ_perm))
            print(f"  JC = CJ? (commute) err = {JC_vs_CJ}")
            
            if JC_vs_CJ < 0.5:  # integer permutations, so exact
                print(f"  J and C commute → S = JC is well-defined")
                print(f"  ε'' = +1 (since JSJ = CJ = JC = S)")
                print(f"  ══════════════════════════════════════════")
                print(f"  ε = +1, ε' = -1, ε'' = +1")
                print(f"  → KO-DIMENSION 1 → CLASS BDI → Z INVARIANT")
                print(f"  ══════════════════════════════════════════")
            else:
                print(f"  J and C do NOT commute")
                # Compute [J,C] details
                for k in range(min(10, DIM)):
                    if JC_perm[k] != CJ_perm[k]:
                        mi, ni = k // phi_q, k % phi_q
                        print(f"    k={k} ({G[mi]},{G[ni]}): JC→{JC_perm[k]}, CJ→{CJ_perm[k]}")
                
        elif err_JCKad_comm < 1e-8:
            print(f"  → S = JC commutes with K_ad (NOT chiral)")
        else:
            print(f"  → S = JC neither commutes nor anticommutes with K_ad")
    
    elif err_JKadJ_plus < 1e-8:
        print(f"  J commutes with K_ad → J is a standard symmetry, not particle-hole")
        print(f"  Different classification needed")
    
    else:
        print(f"  J neither commutes nor anticommutes with K_ad")
        print(f"  K_ad breaks the Tomita-Takesaki structure")
    
    print()
    
    # ═══════════════════════════════════════════════════════════
    # CONSTRUCTION 3: p-adic components
    # ═══════════════════════════════════════════════════════════
    
    print("--- p-ADIC DECOMPOSITION ---")
    
    # For primes p dividing q, the p-adic valuation of elements in G = (Z/qZ)*
    # is always 0 (since gcd(g, q) = 1). So K_p = 0 for primes dividing q.
    # For primes NOT dividing q, v_p(g) = 0 for all g in G (since g < q and gcd(g,q)=1
    # doesn't mean g isn't divisible by other primes... wait)
    
    # Actually, for g ∈ (Z/qZ)*, gcd(g, q) = 1, so v_p(g) = 0 for primes p|q.
    # But g can have factors of primes not dividing q.
    # E.g., q=30, G includes 7, 11, 13, etc. v_7(7) = 1, v_7(49) = 2, etc.
    # But our G only contains elements < q, so we can compute directly.
    
    # K_p on |m><n| has eigenvalue (v_p(n) - v_p(m)) log(p)
    # Sum over all primes: Σ_p (v_p(n) - v_p(m)) log(p) = log(n/m) = K!
    # (This is the fundamental theorem of arithmetic.)
    
    # So K = Σ_p K_p exactly.
    # The question is: which K_p's does C preserve?
    
    print(f"  K = Σ_p K_p where K_p(m,n) = (v_p(n) - v_p(m)) log p")
    print(f"  Verifying K = Σ_p K_p...")
    
    # Find relevant primes (those that appear in factorizations of elements of G)
    primes_used = set()
    for g in G:
        if g > 1:
            for p in factorint(g):
                primes_used.add(p)
    primes_used = sorted(primes_used)
    print(f"  Primes in factorizations: {primes_used}")
    
    K_sum_diag = np.zeros(DIM)
    K_p_data = {}
    
    for p in primes_used:
        Kp_diag = np.zeros(DIM)
        for mi in range(phi_q):
            for ni in range(phi_q):
                vp_n = p_adic_val(G[ni], p)
                vp_m = p_adic_val(G[mi], p)
                Kp_diag[mi * phi_q + ni] = (vp_n - vp_m) * log(p)
        K_sum_diag += Kp_diag
        K_p_data[p] = Kp_diag
        
        # Check CK_pC vs K_p
        CKpC_diag = Kp_diag[C_perm]
        err_p_anti = np.max(np.abs(CKpC_diag + Kp_diag))
        err_p_comm = np.max(np.abs(CKpC_diag - Kp_diag))
        
        # Check JK_pJ vs K_p
        JKpJ_diag = Kp_diag[J_perm]
        err_p_J_anti = np.max(np.abs(JKpJ_diag + Kp_diag))
        err_p_J_comm = np.max(np.abs(JKpJ_diag - Kp_diag))
        
        n_nonzero_p = np.sum(np.abs(Kp_diag) > 1e-10)
        
        print(f"  p={p:3d}: nonzero={n_nonzero_p:4d}/{DIM}  "
              f"CK_pC={'−K_p' if err_p_anti < 1e-8 else '+K_p' if err_p_comm < 1e-8 else '???':>4s} "
              f"(err={min(err_p_anti,err_p_comm):.2e})  "
              f"JK_pJ={'−K_p' if err_p_J_anti < 1e-8 else '+K_p' if err_p_J_comm < 1e-8 else '???':>4s} "
              f"(err={min(err_p_J_anti,err_p_J_comm):.2e})")
    
    # Verify sum
    sum_err = np.max(np.abs(K_sum_diag - K_diag))
    print(f"  Σ K_p = K? err = {sum_err:.2e}")
    print()
    
    # ═══════════════════════════════════════════════════════════
    # CONSTRUCTION 4: Weighted adelic Hamiltonian
    # K_w = Σ_p w_p K_p  with weights chosen to make CK_wC = -K_w
    # ═══════════════════════════════════════════════════════════
    
    if len(primes_used) > 1:
        print("--- CONSTRUCTION 4: Weighted p-adic Hamiltonian ---")
        print(f"  Can we find weights w_p such that C(Σ w_p K_p)C = -Σ w_p K_p?")
        
        # This requires: Σ w_p CK_pC = -Σ w_p K_p
        # i.e., Σ w_p (CK_pC + K_p) = 0
        # 
        # Let D_p = CK_pC + K_p (the "defect" of each K_p under C)
        # We need: Σ w_p D_p = 0 (as diagonal operators)
        
        # Stack D_p diagonals as columns
        n_primes = len(primes_used)
        D_matrix = np.zeros((DIM, n_primes))
        for j, p in enumerate(primes_used):
            Dp = K_p_data[p][C_perm] + K_p_data[p]
            D_matrix[:, j] = Dp
        
        # Find null space of D_matrix^T (weight vectors w such that D_matrix @ w = 0)
        U, S_svd, Vt = np.linalg.svd(D_matrix, full_matrices=True)
        tol = 1e-8
        null_dim = np.sum(S_svd < tol)
        
        print(f"  Number of primes: {n_primes}")
        print(f"  SVD singular values: {np.round(S_svd, 6)}")
        print(f"  Null space dimension: {null_dim}")
        
        if null_dim > 0:
            # Extract null vectors (last rows of Vt)
            null_vectors = Vt[n_primes - null_dim:, :]
            print(f"  Null vector(s):")
            for k in range(null_dim):
                w = null_vectors[k]
                print(f"    w = [{', '.join(f'{v:.4f}' for v in w)}]")
                print(f"    Primes: {primes_used}")
                
                # Build K_w = Σ w_p K_p
                Kw_diag = np.zeros(DIM)
                for j, p in enumerate(primes_used):
                    Kw_diag += w[j] * K_p_data[p]
                
                # Verify CK_wC = -K_w
                CKwC_diag = Kw_diag[C_perm]
                err_Kw = np.max(np.abs(CKwC_diag + Kw_diag))
                print(f"    CK_wC = -K_w? err = {err_Kw:.2e}")
                
                # Check JK_wJ
                JKwJ_diag = Kw_diag[J_perm]
                err_Kw_J = np.max(np.abs(JKwJ_diag + Kw_diag))
                print(f"    JK_wJ = -K_w? err = {err_Kw_J:.2e}")
                
                # Gap?
                nz = np.abs(Kw_diag[np.abs(Kw_diag) > 1e-10])
                if len(nz) > 0:
                    print(f"    Gap: {np.min(nz):.8f}")
                else:
                    print(f"    K_w = 0 (trivial)")
        else:
            print(f"  No null vector → no weighted combination makes CK_wC = -K_w")
            print(f"  The individual p-adic defects D_p are linearly independent")
        
        print()

# ═══════════════════════════════════════════════════════════════
# SYNTHESIS
# ═══════════════════════════════════════════════════════════════

print("=" * 85)
print("  GRAND SYNTHESIS")
print("=" * 85)
print("""
RESULTS ACROSS ALL q:

1. K_ad = (K - CKC)/2 EXISTS and satisfies:
   - CK_adC = -K_ad  (by construction, exact)
   - {S_L, K_ad} = 0  (because K_ad is antisymmetric in (m,n), exact)
   - K_ad has non-trivial spectrum with a gap
   
2. THE CRITICAL QUESTION: Does JK_adJ = -K_ad?
   If YES → J and C BOTH anticommute with K_ad
         → S = JC is a LINEAR chiral operator
         → BDI classification with Z invariant
         
   If NO  → K_ad breaks TT compatibility
          → Need different framework

3. ARITHMETIC INTERPRETATION:
   K_ad(m,n) = (1/2)[log(n/m) + log(m_inv/n_inv)]
             = average of archimedean and modular log-ratios
   This is a "half-adelic" Hamiltonian that sees BOTH
   the real arithmetic (through log) and the modular arithmetic
   (through mod-q inverses).

4. p-ADIC DECOMPOSITION:
   K = Σ_p K_p (fundamental theorem of arithmetic)
   Each K_p anticommutes with J (exact: JK_pJ = -K_p)
   but NOT with C in general.
   
   K_ad = the C-antisymmetric projection selects the part
   of K that is compatible with both J and C simultaneously.
""")
