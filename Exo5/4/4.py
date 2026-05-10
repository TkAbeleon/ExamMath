"""
Exploration optimisée de |a^x - b^y| = k, k>=2.
Recherche avec bases jusqu'à 1000, exposants >=2, valeur max paramétrable.
Utilise un dictionnaire pour éviter O(N²).
"""

import math
from collections import defaultdict
import time

# Paramètres ajustables
MAX_BASE = 1000          # bases a, b de 2 à MAX_BASE
MAX_EXP = 50             # exposant maximal (peu sera atteint pour les grandes bases)
MAX_VAL = 10**15         # générer a^x jusqu'à cette valeur
MAX_K = 50               # rechercher k de 2 à MAX_K

def generate_power_dict(max_base, max_exp, max_val):
    """
    Génère un dictionnaire { valeur : liste de (base, exposant) }
    pour toutes les puissances base^exposant avec 2 ≤ base ≤ max_base,
    2 ≤ exposant ≤ max_exp et valeur ≤ max_val.
    """
    power_dict = defaultdict(list)
    for a in range(2, max_base + 1):
        val = a * a          # exposant 2
        exp = 2
        while val <= max_val and exp <= max_exp:
            power_dict[val].append((a, exp))
            val *= a
            exp += 1
    return power_dict

def find_solutions(max_k):
    print(f"Génération des puissances (base max {MAX_BASE}, valeur max {MAX_VAL})...")
    start = time.time()
    pwr = generate_power_dict(MAX_BASE, MAX_EXP, MAX_VAL)
    elapsed = time.time() - start
    print(f"Nombre de valeurs distinctes : {len(pwr)} (temps : {elapsed:.2f} s)")

    # Récupération et tri des valeurs
    sorted_values = sorted(pwr.keys())
    sol_by_k = defaultdict(list)

    # Parcours des valeurs triées, en regardant pour chaque v si v+k est présent
    for i, v in enumerate(sorted_values):
        # Pour chaque k, on vérifie si v+k existe
        for k in range(2, max_k + 1):
            target = v + k
            if target > sorted_values[-1]:
                break  # on dépasse la plus grande valeur
            # Recherche dichotomique ou accès dict
            if target in pwr:
                # On a trouvé une différence = k
                # On vérifie tous les couples (a1,x1) pour v et (a2,x2) pour target
                for a1, x1 in pwr[v]:
                    for a2, x2 in pwr[target]:
                        if math.gcd(a1, a2) == 1:
                            # Éviter les doublons en stockant le couple ordonné
                            sol_by_k[k].append((a1, x1, a2, x2, v, target))
                            # On peut arrêter pour ne pas surcharger l'affichage
                            break
                    if sol_by_k[k]:
                        break
    return sol_by_k

def main():
    solutions = find_solutions(MAX_K)
    total_sols = 0
    for k in sorted(solutions):
        lst = solutions[k]
        total_sols += len(lst)
        print(f"\n--- k = {k} ({len(lst)} solution(s)) ---")
        for a1, x1, a2, x2, v1, v2 in lst[:5]:  # affiche au plus 5 exemples
            print(f"  {a1}^{x1} = {v1},  {a2}^{x2} = {v2}   (pgcd={math.gcd(a1,a2)})")
        if len(lst) > 5:
            print(f"  ... et {len(lst)-5} autre(s)")
    print(f"\nNombre total de solutions trouvées (k=2..{MAX_K}) : {total_sols}")

if __name__ == "__main__":
    main()
