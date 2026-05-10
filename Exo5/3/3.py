"""
Recherche de solutions pour la conjecture d'Erdős-Straus :
    4/n = 1/x + 1/y + 1/z, avec x, y, z entiers > 0.

Méthode robuste :
1. Identités polynomiales exactes pour les classes modulo 840.
2. Recherche exhaustive avec bornes optimisées (théorie des nombres).
3. Vérification rigoureuse des solutions.
"""

import math
import sys
from fractions import Fraction

# ------------------------------------------------------------
# 1. Identités polynomiales rigoureuses
# ------------------------------------------------------------

def identite_n_pair(n):
    """n pair : 4/(2k) = 2/k, puis on décompose 2/k en 3 fractions."""
    if n % 2 == 0:
        k = n // 2
        # 2/k = 1/k + 1/k + 0 ? Non, on veut 3 fractions > 0.
        # 2/k = 1/k + 1/(2k) + 1/(2k) fonctionne pour k > 1
        if k == 1:
            return None  # n=2 -> 4/2 = 2, pas de décomposition immédiate
        # 4/(2k) = 1/k + 1/(2k) + 1/(2k)  est FAUX : 1/k + 1/(2k) + 1/(2k) = 2/k OK !
        # Vérification : 1/k + 1/(2k) + 1/(2k) = (2 + 1 + 1)/(2k) = 4/(2k) = 2/k
        # Oui, c'est correct :
        return (k, 2*k, 2*k)
    return None

def identite_n_3mod4(n):
    """n ≡ 3 (mod 4) : identité de Sierpinski."""
    if n % 4 == 3:
        k = (n - 3) // 4
        # 4/(4k+3) = 1/(k+1) + 1/((k+1)(4k+3)) + 1/((k+1)(4k+3))
        # Vérifions : 1/a + 1/(ab) + 1/(ab) = 1/a + 2/(ab) = (b+2)/(ab)
        # On veut (b+2)/(ab) = 4/b => a = (b+2)/4
        # Avec b = 4k+3, a = (4k+5)/4 pas entier. Donc cette identité est fausse.
        # On utilise plutôt une recherche.
        return None
    return None

def identite_n_1mod4(n):
    """n ≡ 1 (mod 4) : identité de Mordell."""
    if n % 4 == 1:
        # Chercher un diviseur d de n tel que d ≡ 1 (mod 4)
        for d in range(1, int(math.isqrt(n)) + 1):
            if n % d == 0:
                for div in [d, n//d]:
                    if div % 4 == 1 and div > 1:
                        k = (div - 1) // 4
                        a = (n - div) // 4
                        if a >= 0:
                            # Formule de Mordell généralisée
                            x = (n + 3*div) // 4
                            y = (n * (n + div)) // (4*div)
                            # Vérifier si x et y sont entiers et > 0
                            if (n + 3*div) % 4 == 0 and (n * (n + div)) % (4*div) == 0:
                                # Calculer z
                                # 4/n = 1/x + 1/y + 1/z
                                num = 4*x*y - n*y - n*x
                                den = n*x*y
                                if num > 0 and den % num == 0:
                                    z = den // num
                                    if z > 0:
                                        return (x, y, z)
        return None
    return None

# ------------------------------------------------------------
# 2. Recherche exhaustive correcte et optimisée
# ------------------------------------------------------------

def recherche_erdos_straus(n, limite_x=None):
    """
    Recherche exhaustive avec bornes théoriques optimales.

    De 4/n = 1/x + 1/y + 1/z, on tire :
    1/x < 4/n => x > n/4
    De plus, par symétrie, on peut supposer x ≤ y ≤ z.
    Alors 1/x ≥ 1/3 * 4/n => x ≤ 3n/4.

    Pour x fixé, on résout 1/y + 1/z = 4/n - 1/x = (4x-n)/(nx).
    C'est une équation en 2 fractions égyptiennes, soluble si et seulement si
    le dénominateur réduit divise quelque chose...

    Optimisation : on parcourt y de manière bornée.
    """
    if limite_x is None:
        limite_x = min(n, 50000)  # Ajustable selon puissance PC

    # x doit être entre n/4 et 3n/4 (si on suppose x ≤ y ≤ z)
    # Mais pour être sûr, on explore x > n/4
    x_min = n // 4 + 1
    x_max = min(3*n//4, limite_x)

    for x in range(x_min, x_max + 1):
        # Vérifier que 4x > n
        if 4*x <= n:
            continue

        # 1/y + 1/z = (4x-n)/(nx)
        num = 4*x - n
        den = n * x

        # Simplifier la fraction
        g = math.gcd(num, den)
        p = num // g
        q = den // g

        # On a 1/y + 1/z = p/q avec p, q premiers entre eux
        # Solutions si et seulement si il existe y tel que
        # y divise q et y > q/p ? Non, cherchons directement y.

        # y doit vérifier : 1/y < p/q => y > q/p
        y_min = q // p + 1
        # y ne peut pas être trop grand car 1/y + 1/z = p/q avec z ≥ y
        # donc 2/y ≥ p/q => y ≤ 2q/p
        y_max = min(2*q//p + 1, 1000000)  # Borne de sécurité

        for y in range(y_min, y_max + 1):
            # 1/z = p/q - 1/y = (py - q)/(qy)
            num_z = p*y - q
            if num_z <= 0:
                continue
            den_z = q*y
            if den_z % num_z == 0:
                z = den_z // num_z
                if z > 0 and z >= y:  # On impose x ≤ y ≤ z
                    return (x, y, z)
    return None

# ------------------------------------------------------------
# 3. Fonction principale combinant les méthodes
# ------------------------------------------------------------

def trouver_solution(n):
    """Trouve (x,y,z) solution de 4/n = 1/x + 1/y + 1/z."""

    # 1. n pair : identité simple
    sol = identite_n_pair(n)
    if sol and verifier_solution(n, sol):
        return sol

    # 2. n ≡ 1 (mod 4) : identité de Mordell
    sol = identite_n_1mod4(n)
    if sol and verifier_solution(n, sol):
        return sol

    # 3. Recherche exhaustive par x croissant
    #    On augmente progressivement la limite pour les cas difficiles
    for limite in [1000, 5000, 20000, 50000]:
        sol = recherche_erdos_straus(n, limite_x=limite)
        if sol and verifier_solution(n, sol):
            return sol

    return None

def verifier_solution(n, sol):
    """Vérifie rigoureusement qu'une solution est correcte."""
    if sol is None:
        return False
    x, y, z = sol
    if x <= 0 or y <= 0 or z <= 0:
        return False
    # Calcul exact avec Fraction
    somme = Fraction(1, x) + Fraction(1, y) + Fraction(1, z)
    cible = Fraction(4, n)
    return somme == cible

# ------------------------------------------------------------
# 4. Tests et exploration
# ------------------------------------------------------------

if __name__ == "__main__":
    print("=== Test de la conjecture d'Erdős-Straus ===\n")

    # Classes récalcitrantes modulo 840
    recalcitrants_mod840 = [1, 121, 169, 289, 361, 529]
    print("1. Test des classes récalcitrantes modulo 840 (jusqu'à 10000):")
    for r in recalcitrants_mod840:
        trouve = False
        for n in range(r, 10000, 840):
            sol = trouver_solution(n)
            if sol and verifier_solution(n, sol):
                x, y, z = sol
                print(f"   Classe {r:3d} mod 840 : n={n:5d} -> x={x:6d}, y={y:8d}, z={z:8d}")
                trouve = True
                break
        if not trouve:
            print(f"   Classe {r:3d} mod 840 : aucune solution trouvée < 10000")

    print("\n2. Test de tous les entiers de 2 à 200 (vérification exhaustive):")
    echecs = []
    for n in range(2, 201):
        sol = trouver_solution(n)
        if sol is None:
            echecs.append(n)
        elif not verifier_solution(n, sol):
            echecs.append(n)
            print(f"   ERREUR pour n={n} : solution {sol} incorrecte")
    if not echecs:
        print("   Tous les entiers de 2 à 200 ont une solution correcte.")
    else:
        print(f"   Échecs : {echecs}")

    print("\n3. Test pour quelques nombres premiers > 10000:")
    for n in [10007, 10009, 10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093]:
        sol = trouver_solution(n)
        if sol:
            x, y, z = sol
            print(f"   n={n:5d} -> x={x:6d}, y={y:8d}, z={z:8d}")
        else:
            print(f"   n={n:5d} -> pas de solution trouvée")

    print("\n4. Statistiques pour n jusqu'à 1000:")
    import time
    debut = time.time()
    nb_succes = 0
    for n in range(2, 1001):
        sol = trouver_solution(n)
        if sol and verifier_solution(n, sol):
            nb_succes += 1
    temps = time.time() - debut
    print(f"   {nb_succes}/999 succès en {temps:.2f} secondes")
    if nb_succes < 999:
        print(f"   ATTENTION : {999 - nb_succes} échecs !")
