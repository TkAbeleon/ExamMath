#!/usr/bin/env python3
"""
Calcul de J = ∫_0^∞ exp(-(x + sin(x))) dx
Intégration directe sur [0, ∞) avec scipy et mpmath.
"""

import numpy as np
from scipy.integrate import quad
import mpmath as mp

# ----------------------------------------------------------------------
# 1. Méthode scipy (quad) – supporte les bornes infinies
# ----------------------------------------------------------------------
def f(x):
    return np.exp(-(x + np.sin(x)))

val_quad, err_quad = quad(f, 0, np.inf, limit=500, epsabs=1e-12, epsrel=1e-12)

print("=== Intégration directe de J = ∫_0^∞ exp(-(x+sin x)) dx ===\n")
print(f"[scipy.quad]  J ≈ {val_quad:.15f}   (erreur estimée = {err_quad:.2e})")

# ----------------------------------------------------------------------
# 2. Vérification avec mpmath (précision arbitraire)
# ----------------------------------------------------------------------
mp.mp.dps = 40   # 40 chiffres

def f_mp(t):
    return mp.e**(-(t + mp.sin(t)))

val_mp = mp.quad(f_mp, [0, mp.inf])
print(f"[mpmath]      J ≈ {val_mp}")

# Comparaison (conversion explicite)
diff = abs(float(val_mp) - val_quad)
print(f"\nÉcart scipy vs mpmath : {diff:.2e}")

# ----------------------------------------------------------------------
# 3. Valeur finale retenue
# ----------------------------------------------------------------------
print("\n=== Valeur approchée retenue ===")
print(f"J = {float(val_mp):.15f}")