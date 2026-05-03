#!/usr/bin/env python3
"""
Validation de la série K = S + Σ b_n I_{n-1}
avec I_m évalué de manière stable par mp.quad.
Accélération d'Aitken pour estimer la limite.
"""

import mpmath as mp

mp.mp.dps = 50

# ----------------------------------------------------------------------
# 1. S = ∫_0^1 sin(x)/x dx (série rapide)
# ----------------------------------------------------------------------
def compute_S():
    return mp.nsum(lambda k: (-1)**k / ((2*k+1) * mp.factorial(2*k+1)),
                   [0, mp.inf])

# ----------------------------------------------------------------------
# 2. Coefficients b_n (récurrence stable)
# ----------------------------------------------------------------------
def coefficients_b(N):
    b = [mp.mpf(1)]  # b_0
    for n in range(1, N+1):
        s = mp.mpf(0)
        for j in range(1, n+1):
            a_j = (-1)**j / (j+1)
            s += a_j * b[n-j]
        b.append(-s)
    return b

# ----------------------------------------------------------------------
# 3. I_m par intégration directe (stable)
# ----------------------------------------------------------------------
def I_m_quad(m):
    """Retourne ∫_0^1 x^m sin x dx via mp.quad (évaluation robuste)."""
    if m == 0:
        return 1 - mp.cos(1)
    if m == 1:
        return mp.sin(1) - mp.cos(1)
    return mp.quad(lambda x: x**m * mp.sin(x), [0, 1])

# ----------------------------------------------------------------------
# 4. Somme partielle K_N
# ----------------------------------------------------------------------
def K_partial(N):
    S = compute_S()
    b = coefficients_b(N)
    total = S
    for n in range(1, N+1):
        total += b[n] * I_m_quad(n-1)
    return total

# ----------------------------------------------------------------------
# 5. Accélération d'Aitken
# ----------------------------------------------------------------------
def aitken(seq):
    s0, s1, s2 = seq[0], seq[1], seq[2]
    return s0 - (s1 - s0)**2 / (s2 - 2*s1 + s0)

# ----------------------------------------------------------------------
# Référence
# ----------------------------------------------------------------------
ref = mp.quad(lambda x: mp.sin(x) / mp.log(1+x), [0, 1])
print("Référence (quad) :", ref)
print()

# Calcul de quelques sommes partielles (N petit pour limiter le temps)
N_values = [500,750, 1000, 1250, 1500, 1750, 2000]
partials = []
for N in N_values:
    val = K_partial(N)
    partials.append(val)
    print(f"N = {N:3d}   K_N = {val}")

# On prend les trois dernières sommes pour Aitken
acc = aitken(partials[-3:])
print(f"\nAccélération d'Aitken avec N={N_values[-3:]}:")
print(f"K ≈ {acc}")
print(f"Écart / référence : {abs(acc - ref)}")