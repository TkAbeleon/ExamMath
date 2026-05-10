import numpy as np
from scipy.special import sici   # contient si et ci
import matplotlib.pyplot as plt

a = np.log(2)

def f_analytique(t):
    """
    Calcule f(t) = (1/π) * [ -cos(a t) Ci(a |t|) + sin(a |t|) (π/2 - Si(a |t|)) ]
    """
    x = a * np.abs(t)            # a |t|
    if x == 0:
        return np.inf            # singularité logarithmique en 0
    si_val, ci_val = sici(x)     # renvoie (Si(x), Ci(x))
    # si_val = Si(x), ci_val = Ci(x)
    terme1 = -np.cos(a * t) * ci_val
    terme2 =  np.sin(a * np.abs(t)) * (np.pi / 2 - si_val)
    return (terme1 + terme2) / np.pi

# Vectorisation pour travailler sur des tableaux
f_vec = np.vectorize(f_analytique)

# Intervalle de tracé (on évite la singularité exacte en 0)
t_vals = np.linspace(-5, 5, 1001)
t_vals = t_vals[t_vals != 0]    # retire le point exact 0 pour éviter inf dans le graphe
f_vals = f_vec(t_vals)

# Graphique
plt.figure(figsize=(10, 5))
plt.plot(t_vals, f_vals, label=r'$f(t)$ (formule analytique)')
plt.xlabel('t')
plt.ylabel('f(t)')
plt.title(r'Transformée inverse de $F(s)=\frac{1}{\ln 2 + |s|}$')
plt.grid(True)
plt.legend()
plt.xlim(t_vals[0], t_vals[-1])
plt.show()

# Affichage de quelques valeurs numériques
test_t = np.array([0.1, 0.5, 1.0, 2.0, 3.0])
print("t\t f(t) analytique")
for t in test_t:
    print(f"{t:.2f}\t {f_analytique(t):.8f}")
