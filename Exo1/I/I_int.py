#!/usr/bin/env python3
"""
Calcul de I = ∫_0^1 cos(x² + e^x) dx
par plusieurs méthodes numériques (scipy, numpy, Monte Carlo).
Version corrigée (utilise trapezoid au lieu de trapz déprécié).
"""

import numpy as np
from scipy import integrate

# ----------------------------------------------------------------------
# Fonction à intégrer
# ----------------------------------------------------------------------
def f(x):
    return np.cos(x**2 + np.exp(x))

# ----------------------------------------------------------------------
# Méthodes de calcul
# ----------------------------------------------------------------------
def method_quad(a, b):
    """Méthode 1 : scipy.integrate.quad (adaptative, référence)"""
    val, err = integrate.quad(f, a, b, limit=200, epsabs=1e-12, epsrel=1e-12)
    return val, err

def method_fixed_quad(a, b, n=200):
    """Méthode 2 : quadrature de Gauss-Legendre à n points fixes"""
    val, _ = integrate.fixed_quad(f, a, b, n=n)
    return val

def method_simpson(a, b, N=10_000):
    """Méthode 3 : règle de Simpson sur N intervalles"""
    x = np.linspace(a, b, N + 1)
    y = f(x)
    val = integrate.simpson(y, x)
    return val



def method_trapezoidal(a, b, N=10_000):
    """Méthode 5 : règle des trapèzes (utilise trapezoid si dispo, sinon trapz)"""
    x = np.linspace(a, b, N + 1)
    y = f(x)
    # Choix compatible avec les anciennes et nouvelles versions de NumPy
    if hasattr(np, 'trapezoid'):
        val = np.trapezoid(y, x)
    else:
        val = np.trapz(y, x)  # pragma: no cover
    return val

def method_montecarlo(a, b, N=100_000, seed=42):
    """Méthode 6 : Monte Carlo simple (échantillonnage uniforme)"""
    rng = np.random.default_rng(seed)
    x_sample = rng.uniform(a, b, N)
    y_sample = f(x_sample)
    val = (b - a) * np.mean(y_sample)
    err = (b - a) * np.std(y_sample, ddof=1) / np.sqrt(N)
    return val, err

# ----------------------------------------------------------------------
# Programme principal
# ----------------------------------------------------------------------
def main():
    a, b = 0.0, 1.0
    print("Calcul de I = ∫_0^1 cos(x² + e^x) dx\n")

    ref, err_ref = method_quad(a, b)
    print(f"{'Méthode':<25} {'Résultat':<20} {'Erreur estimée':<16} {'Différence / ref'}")
    print("-" * 75)
    print(f"{'quad (référence)':<25} {ref:<20.15f} ±{err_ref:<14.2e} {'---':<10}")

    methods = [
        ("Gauss-Legendre (200pts)", lambda: method_fixed_quad(a, b, n=200)),
        ("Simpson (10000 interv.)", lambda: method_simpson(a, b, N=10_000)),
        ("Trapèzes (10000 interv.)", lambda: method_trapezoidal(a, b, N=10_000)),
        ("Monte Carlo (100k éch.)", lambda: method_montecarlo(a, b, N=100_000)),
    ]

    for name, func in methods:
        if "Monte Carlo" in name:
            val, err = func()
            diff = abs(val - ref)
            print(f"{name:<25} {val:<20.15f} ±{err:<14.2e} {diff:<10.2e}")
        else:
            val = func()
            diff = abs(val - ref)
            print(f"{name:<25} {val:<20.15f} {'---':<16} {diff:<10.2e}")

if __name__ == '__main__':
    main()
