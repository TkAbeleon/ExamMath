import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as plt

# Domaine s > s_min où s_min = -1 + exp(-1/e)
s_min = -1 + np.exp(-1/np.e) + 1e-6
s_vals = np.linspace(s_min, 5, 1000)

# Calcul de F(s) = W(ln(s+1)), branche principale réelle
ln_s1 = np.log(s_vals + 1)
F_vals = lambertw(ln_s1).real

# Dérivée exacte (formule close) : F'(s) = F(s) / ((s+1) ln(s+1) (1+F(s)))
# En s=0, la limite est 1.
with np.errstate(divide='ignore', invalid='ignore'):
    F_prime_exact = np.where(
        np.abs(s_vals) < 1e-12,
        1.0,
        F_vals / ((s_vals + 1) * ln_s1 * (1 + F_vals))
    )

# ---------- Vérifications ----------
err_orig = np.max(np.abs(F_vals * np.exp(F_vals) - ln_s1))
print(f"Erreur max sur F(s) exp(F(s)) = ln(s+1) : {err_orig:.2e}")

# Équation différentielle sous sa forme originale : e^{F} (1+F) F' = 1/(s+1)
lhs_diff = np.exp(F_vals) * (1 + F_vals) * F_prime_exact
rhs_diff = 1 / (s_vals + 1)
err_diff = np.max(np.abs(lhs_diff - rhs_diff))
print(f"Erreur max sur l'équation différentielle (forme originale) : {err_diff:.2e}")

# ---------- Affichage et export ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(s_vals, F_vals, 'b-', linewidth=2, label=r'$F(s) = W(\ln(s+1))$')
ax1.set_xlabel('s')
ax1.set_ylabel('F(s)')
ax1.set_title('Fonction F(s)')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend()

ax2.plot(s_vals, F_prime_exact, 'r-', linewidth=2, label=r"$F'(s)$ (exacte)")
ax2.set_xlabel('s')
ax2.set_ylabel("F'(s)")
ax2.set_title("Dérivée de F(s)")
ax2.grid(True, linestyle='--', alpha=0.6)
ax2.legend()

plt.tight_layout()
plt.savefig('Exo2_plot.png', dpi=300, bbox_inches='tight')
plt.show()