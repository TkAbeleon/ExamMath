#!/usr/bin/env python3
"""
Calcul de J = ∫_0^∞ exp(-(x+sin x)) dx
avec une très haute précision (100 ou 200 décimales) en utilisant
la série exacte issue de la linéarisation de sin^n x.
"""

import mpmath as mp

# ---------------------------- PARAMÈTRES ---------------------------------
DECIMALS = 100  # Nombre de décimales souhaité (modifiable)
# ------------------------------------------------------------------------

mp.mp.dps = DECIMALS + 5  # quelques chiffres de garde

def I_2m(m):
    """I_{2m} exact (pair) en précision arbitraire."""
    m_val = mp.mpf(m)
    terme_central = mp.binomial(2*m, m) / (2**(2*m))
    somme = mp.mpf('0')
    for r in range(m):
        signe = mp.mpf(-1)**(m - r)
        coeff = mp.binomial(2*m, r)
        denom = 1 + 4 * (m - r)**2
        somme += signe * coeff / denom
    somme /= 2**(2*m - 1)
    return terme_central + somme

def I_2m_plus_1(m):
    """I_{2m+1} exact (impair) en précision arbitraire."""
    somme = mp.mpf('0')
    for r in range(m + 1):
        signe = mp.mpf(-1)**(m - r)
        coeff = mp.binomial(2*m + 1, r)
        a = mp.mpf(2*m + 1 - 2*r)
        somme += signe * coeff * a / (1 + a**2)
    return somme / (2**(2*m))

def J_haute_precision(termes_max=200, tolerance=None):
    """
    Somme la série J = Σ (I_{2m}/(2m)! - I_{2m+1}/(2m+1)!).
    S'arrête si le terme courant est plus petit que `tolerance` ou si
    `termes_max` est atteint.
    """
    total = mp.mpf('0')
    for m in range(termes_max + 1):
        terme = (I_2m(m) / mp.factorial(2*m)
                 - I_2m_plus_1(m) / mp.factorial(2*m + 1))
        total += terme
        if tolerance is not None and abs(terme) < tolerance:
            break
    return total

if __name__ == '__main__':
    print(f"Calcul de J avec {DECIMALS} décimales...\n")
    # On fixe la tolérance à 10^{-(DECIMALS+2)} pour avoir DECIMALS chiffres exacts
    tolerance = mp.mpf('1e-{}'.format(DECIMALS + 2))
    J_val = J_haute_precision(termes_max=200, tolerance=tolerance)

    # Affichage avec le nombre de décimales voulu
    mp.nprint(J_val, DECIMALS)
    print(f"\nNombre de chiffres exacts affichés : {DECIMALS}")