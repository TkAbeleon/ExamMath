import math

def extended_gcd(a, b):
    """
    Algorithme d'Euclide étendu.
    Retourne (d, u, v) tel que d = pgcd(a, b) et au + bv = d.
    """
    if b == 0:
        return a, 1, 0
    d, u1, v1 = extended_gcd(b, a % b)
    u = v1
    v = u1 - (a // b) * v1
    return d, u, v

def modinv(a, m):
    """
    Inverse modulaire de a modulo m (si pgcd(a,m)=1).
    """
    d, u, _ = extended_gcd(a, m)
    if d != 1:
        raise ValueError(f"Pas d'inverse modulaire car pgcd({a},{m}) = {d} ≠ 1")
    return u % m

def solve_congruence(a, b, n):
    """
    Résout l'équation a*x ≡ b (mod n).
    Retourne une liste de toutes les solutions distinctes modulo n.
    """
    # 1. Calcul du pgcd(a, n)
    d = math.gcd(a, n)
    if b % d != 0:
        return []  # Aucune solution

    # 2. Simplification
    a_prime = a // d
    b_prime = b // d
    n_prime = n // d

    # 3. Solution particulière de a'*x ≡ b' (mod n')
    # a' et n' sont premiers entre eux, on cherche l'inverse de a'
    inv_a = modinv(a_prime, n_prime)
    x0 = (b_prime * inv_a) % n_prime

    # 4. Toutes les solutions modulo n
    solutions = []
    for k in range(d):
        x = x0 + k * n_prime
        solutions.append(x)

    return solutions

# ---------------------------
# Application à l'exercice
# ---------------------------
a = 750005
b = 1500010
n = 9876545

solutions = solve_congruence(a, b, n)

print(f"Résolution de {a} x ≡ {b} (mod {n})")
print(f"Nombre de solutions : {len(solutions)}")
print("Solutions :")
for i, x in enumerate(solutions):
    # Vérification
    assert (a * x - b) % n == 0
    print(f"  x_{i+1} = {x} (mod {n})")
