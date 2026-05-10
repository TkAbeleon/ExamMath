from sympy.ntheory.modular import crt

# Modules et restes
moduli = [7, 11, 13, 17, 23]
remainders = [3, 1, 5, 15, 12]

# Résolution directe
remainder, modulus = crt(moduli, remainders)

print(f"Solution : x ≡ {remainder} (mod {modulus})")
# Vérification
for r, m in zip(remainders, moduli):
    assert remainder % m == r, f"Erreur pour le modulo {m}"
print("Vérification : toutes les congruences sont satisfaites.")
