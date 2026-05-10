"""
G(s) = ln|s| / sqrt(1+s^2)   (pour s > 0)
"""
import numpy as np
import matplotlib.pyplot as plt

def G(s):
    """Calcule G(s) pour s > 0 (réel)"""
    return np.log(s) / np.sqrt(1 + s**2)

# Analyse et description
print("=== Analyse de G(s) ===")
print("Forme : G(s) = ln(s) / sqrt(1+s^2)")
print("Domaine : s > 0 (réel) ; pas analytique sur C à cause de |s|")
print("Comportement :")
print("  - s -> 0+ : ln(s) -> -∞, donc G(s) -> -∞ lentement.")
print("  - s -> +∞ : G(s) ~ ln(s)/s -> 0+.")
print("Dérivée : G'(s) = (1 - s*ln(s)/(1+s^2)) / (s*sqrt(1+s^2))")
print("Singularité logarithmique en s = 0.")

# Tracé
s_vals = np.logspace(-2, 2, 200)
G_vals = G(s_vals)

plt.figure(figsize=(6,4))
plt.plot(s_vals, G_vals)
plt.xscale('log')
plt.xlabel('s')
plt.ylabel('G(s)')
plt.title(r'$G(s) = \frac{\ln s}{\sqrt{1+s^2}}$')
plt.grid(True)
plt.tight_layout()
plt.savefig('G_s.png')
plt.show()
