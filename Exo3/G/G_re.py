"""
Calcul de g(t) = L^{-1}{ ln(s)/sqrt(1+s^2) }
Formule exacte : g(t) = -J0(t)*ln(t) + Σ (-1)^n/(4^n*(n!)^2) * t^{2n} * ψ(2n+1)
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import j0, psi  # Bessel J0 et digamma

def g_series(t, N=30):
    """Calcule g(t) pour un scalaire t > 0 avec N termes."""
    if t <= 0:
        raise ValueError("t doit être > 0")
    val = -j0(t) * np.log(t)
    for n in range(N + 1):
        coeff = (-1)**n / ( (4**n) * (math.factorial(n)**2) )
        val += coeff * (t**(2*n)) * psi(2*n + 1)
    return val

# Vectorisation pour array
g_t_vec = np.vectorize(g_series, excluded=['N'])

# Analyse et description
print("=== Analyse de g(t) ===")
print("g(t) = L^{-1}{ ln(s)/sqrt(1+s^2) }")
print("Expression exacte (série convergente) :")
print("g(t) = -J0(t) ln(t) + Σ_{n=0}^∞ (-1)^n/(4^n (n!)^2) t^{2n} ψ(2n+1)")
print("Propriétés :")
print("  - t -> 0+ : g(t) ~ -ln(t) - γ  (J0(0)=1)")
print("    → g(t) → +∞")
print("  - t modéré : oscillations amorties.")
print("  - t -> +∞ : décroissance lente avec oscillations.")
print("N.B. : La série converge bien pour t < ~20 avec N=30")

# Exemple numérique
t_vals = np.linspace(0.01, 15, 200)
g_vals = g_t_vec(t_vals, N=30)

plt.figure(figsize=(6,4))
plt.plot(t_vals, g_vals)
plt.xlabel('t')
plt.ylabel('g(t)')
plt.title('Transformée inverse de $G(s)=\\ln s/\\sqrt{1+s^2}$')
plt.grid(True)
plt.tight_layout()
plt.savefig('g_t.png')
plt.show()
