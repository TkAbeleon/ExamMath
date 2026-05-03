#!/usr/bin/env python3
"""
Calcul de I = ∫_0^1 cos(x² + eˣ) dx
via la série double exacte, avec quadrature multi‑précision pour éviter
les instabilités numériques.
"""

import mpmath as mp

mp.mp.dps = 50   # 50 chiffres décimaux

def I_quad(m, a):
    """∫_0^1 x^m e^{a x} dx évaluée par quadrature robuste."""
    if a == 0:
        return mp.mpf(1) / (m + 1)
    f = lambda x: x**m * mp.e**(a * x)
    return mp.quad(f, [0, 1])

def serie_double(N_max):
    """Somme la série double jusqu'à n = N_max."""
    total = mp.mpf('0')
    for n in range(N_max + 1):
        signe = mp.mpf(-1) ** n
        fact_2n = mp.factorial(2 * n)
        somme_k = mp.mpf('0')
        for k in range(0, 2 * n + 1):
            coeff = mp.binomial(2 * n, k)
            m_val = 2 * k
            a_val = 2 * n - k
            I_val = I_quad(m_val, a_val)
            somme_k += coeff * I_val
        total += signe * somme_k / fact_2n
    return total

def main():
    # Valeur de référence via scipy (optionnelle)
    try:
        import scipy.integrate as spi
        import numpy as np
        ref, err = spi.quad(lambda x: np.cos(x**2 + np.exp(x)), 0, 1)
        print(f"Référence (scipy.quad) : {ref:.15f} ± {err:.2e}\n")
    except ImportError:
        ref = None
        print("scipy non trouvé, pas de référence.\n")

    print("  N   |   Somme partielle (série double)")
    print("----------------------------------------")
    for N in [5, 10, 15, 20, 25, 30]:
        val = serie_double(N)
        print(f" {N:3d}  |   {val}")

    final = serie_double(30)
    print(f"\nValeur finale (N=30) : {final}")
    if ref is not None:
        diff = abs(mp.mpf(ref) - final)
        print(f"Écart / référence    : {diff}")

if __name__ == '__main__':
    main()