import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as plt

# Domaine (borne inférieure où ln(s+1) >= -1/e)
s0 = -1 + np.exp(-1/np.e)
s_vals = np.linspace(s0 + 1e-6, 5, 1000)  # plus de points pour précision

ln_s1 = np.log(s_vals + 1)
F_vals = lambertw(ln_s1).real

# 1. Vérification de l'équation originale
lhs = F_vals * np.exp(F_vals)
err1 = np.max(np.abs(lhs - ln_s1))
print(f"Erreur max sur F(s) exp(F(s)) = ln(s+1) : {err1:.2e}")

# 2. Dérivée numérique (différences centrées, plus précis que gradient sur bords uniformes)
h = s_vals[1] - s_vals[0]
F_prime = np.gradient(F_vals, h)

# Vérification de l'équation différentielle sous forme produit
residu_produit = (s_vals + 1) * ln_s1 * (1 + F_vals) * F_prime - F_vals
err2 = np.max(np.abs(residu_produit))
print(f"Erreur max forme produit : {err2:.2e}")

# Focalisons sur un sous-domaine excluant le voisinage de s=0 (où F(0)=ln(1)=0)
mask = np.abs(s_vals) > 0.1   # on garde |s| > 0.1
err2_masked = np.max(np.abs(residu_produit[mask]))
print(f"Erreur max forme produit pour |s|>0.1 : {err2_masked:.2e}")

# Comparaison avec la forme division (attendue) sur le même sous-domaine
rhs_diff = F_vals[mask] / ((s_vals[mask] + 1) * ln_s1[mask] * (1 + F_vals[mask]))
# Erreur entre dérivée numérique et formule exacte
err_derivative = np.max(np.abs(F_prime[mask] - rhs_diff))
print(f"Erreur max sur la dérivée pour |s|>0.1 : {err_derivative:.2e}")

# Tracé
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(s_vals, F_vals, label='F(s)')
plt.xlabel('s'); plt.ylabel('F(s)'); plt.legend()

plt.subplot(1,2,2)
plt.plot(s_vals, F_prime, label='F\'(s) numérique')
plt.plot(s_vals[mask], rhs_diff, '--', label='Formule exacte (|s|>0.1)')
plt.xlabel('s'); plt.legend()
plt.tight_layout()
plt.show()