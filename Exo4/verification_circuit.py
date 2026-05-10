import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import TransferFunction, step
from scipy.integrate import odeint

# =============================================================================
# 1. Vérification symbolique (méthode des impédances)
# =============================================================================
s = sp.symbols('s', complex=True)
R1, R2, R3 = sp.symbols('R1 R2 R3', positive=True)
L1, L2, C  = sp.symbols('L1 L2 C', positive=True)

# Impédances
Z1 = 1 / (1/R1 + 1/(L1*s))          # R1 // L1
Z2 = 1 / (1/R2 + 1/(L2*s) + C*s)    # R2 // L2 // C
Z_R3 = R3

# Fonction de transfert
H_s = (Z2 + Z_R3) / (Z1 + Z2 + Z_R3)
H_s_simpl = sp.together(H_s)
Num_H, Den_H = H_s_simpl.as_numer_denom()

# Coefficients "manuels"
a3 = L1*R2*R3*C*L2
a2 = R1*R2*R3*C*L2 + L1*L2*(R2+R3)
a1 = R1*L2*(R2+R3) + L1*R2*R3
a0 = R1*R2*R3
b3 = R2*C*L2*L1*(R1+R3)
b2 = L1*L2*(R1+R2+R3) + R1*R2*R3*C*L2
b1 = R1*R2*(L1+L2) + R1*L2*R3 + L1*R2*R3
b0 = R1*R2*R3
Num_manuel = a3*s**3 + a2*s**2 + a1*s + a0
Den_manuel = b3*s**3 + b2*s**2 + b1*s + b0

print("=== Vérification symbolique ===")
if sp.simplify(Num_H*Den_manuel - Den_H*Num_manuel) == 0:
    print("✅ Transfert identique à l'expression manuelle.")
else:
    print("❌ Différence détectée !")

# =============================================================================
# 2. Calcul des conditions initiales exactes pour l'échelon unitaire
# =============================================================================
# S(s) = H(s)/s
S_s = H_s / s
k0 = sp.limit(s * S_s, s, sp.oo)                         # s(0+)
k1 = sp.limit(s * (s*S_s - k0), s, sp.oo)                # s'(0+)
k2 = sp.limit(s * (s**2*S_s - k0*s - k1), s, sp.oo)      # s''(0+)

print(f"\nConditions initiales (t=0+) :")
print(f" s(0+)  = {k0}")
print(f" s'(0+) = {k1}")
print(f" s''(0+)= {k2}")

# =============================================================================
# 3. Vérification numérique avec les bonnes CI
# =============================================================================
subs_dict = {R1: 10e3, R2: 5e3, R3: 2e3,
             L1: 10e-3, L2: 5e-3, C: 0.1e-6}

a3_n = float(a3.subs(subs_dict))
a2_n = float(a2.subs(subs_dict))
a1_n = float(a1.subs(subs_dict))
a0_n = float(a0.subs(subs_dict))
b3_n = float(b3.subs(subs_dict))
b2_n = float(b2.subs(subs_dict))
b1_n = float(b1.subs(subs_dict))
b0_n = float(b0.subs(subs_dict))

# Conditions initiales numériques
s0  = float(k0.subs(subs_dict))
s0p = float(k1.subs(subs_dict))
s0pp= float(k2.subs(subs_dict))

# Système via scipy
num = [a3_n, a2_n, a1_n, a0_n]
den = [b3_n, b2_n, b1_n, b0_n]
sys = TransferFunction(num, den)

t_end = 0.5e-3
t = np.linspace(0, t_end, 5000)
t, y_scipy = step(sys, T=t)

# Équation différentielle avec odeint
def equa_diff(y, t, a3,a2,a1,a0, b3,b2,b1,b0):
    s_val, ds, d2s = y
    # e(t) = 1 pour t>0, dérivées nulles
    d3s = (a0 - b0*s_val - b1*ds - b2*d2s) / b3
    return [ds, d2s, d3s]

# ICI la correction : y0 = [s(0+), s'(0+), s''(0+)]
y0 = [s0, s0p, s0pp]
sol = odeint(equa_diff, y0, t, args=(a3_n,a2_n,a1_n,a0_n,b3_n,b2_n,b1_n,b0_n))
s_ode = sol[:, 0]

# Comparaison
erreur_max = np.max(np.abs(y_scipy - s_ode))
print(f"\nErreur max entre scipy et odeint : {erreur_max:.2e}")
if erreur_max < 1e-6:
    print("✅ L'équation différentielle est parfaitement vérifiée.")
else:
    print("⚠️ Écart non négligeable – vérifier les calculs.")

# Graphique
plt.figure(figsize=(10,6))
plt.plot(t, y_scipy, 'b-', lw=2, label='scipy (TransferFunction)')
plt.plot(t, s_ode, 'r--', lw=2, label='odeint (équa. diff.)')
plt.xlabel('Temps (s)')
plt.ylabel('Sortie s(t)')
plt.title('Réponse à un échelon unitaire (CI corrigées)')
plt.legend()
plt.grid(True)
plt.savefig('reponse_echelon.png')
plt.show()
