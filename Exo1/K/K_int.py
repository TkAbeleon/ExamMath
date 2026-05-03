#!/usr/bin/env python3
"""
Calcul direct (numérique) de K = ∫_0^1 sin(x) / ln(1+x) dx
avec scipy.quad et mpmath.quad.
"""

import numpy as np
from scipy.integrate import quad
import mpmath as mp

# ----------------------------------------------------------------------
# 1. Intégration avec scipy (double précision)
# ----------------------------------------------------------------------
# Définition robuste de la fonction : sin(0)/ln(1) → 1 (limite)
def f_scipy(x):
    # On traite x=0 explicitement pour éviter 0/0
    if np.isscalar(x):
        if x == 0.0:
            return 1.0
        return np.sin(x) / np.log1p(x)  # log1p(x) = ln(1+x)
    else:
        # Cas vectoriel pour les appels internes de quad
        res = np.empty_like(x)
        mask = (x == 0.0)
        res[mask] = 1.0
        res[~mask] = np.sin(x[~mask]) / np.log1p(x[~mask])
        return res

val_quad, err_quad = quad(f_scipy, 0, 1, limit=200, epsabs=1e-14, epsrel=1e-14)
print("=== Intégration directe de K = ∫_0^1 sin(x)/ln(1+x) dx ===\n")
print(f"[scipy.quad]  K ≈ {val_quad:.15f}   (erreur estimée = {err_quad:.2e})")

# ----------------------------------------------------------------------
# 2. Vérification avec mpmath (précision arbitraire)
# ----------------------------------------------------------------------
mp.mp.dps = 50   # 50 chiffres décimaux

def f_mp(x):
    # mpmath gère naturellement la limite en 0 car mp.sin(0)=0 et mp.log(1)=0
    # mais l'évaluation directe 0/0 produit nan. On utilise une petite astuce :
    if x == 0:
        return mp.mpf(1)
    return mp.sin(x) / mp.log(1 + x)

val_mp = mp.quad(f_mp, [0, 1])
print(f"[mpmath]      K ≈ {val_mp}")

# ----------------------------------------------------------------------
# 3. Comparaison et valeur retenue
# ----------------------------------------------------------------------
diff = abs(float(val_mp) - val_quad)
print(f"\nÉcart scipy vs mpmath : {diff:.2e}")

print("\n=== Valeur approchée retenue (15 décimales) ===")
print(f"K = {float(val_mp):.15f}")