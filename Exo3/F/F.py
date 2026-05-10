import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

a = np.log(2)

def f_quad(t):
    """
    Calcule f(t) = (1/π) ∫_0^∞ cos(ωt) / (a + ω) dω
    en utilisant la quadrature de Fourier de QUAD (weight='cos').
    """
    if np.abs(t) < 1e-12:
        return np.inf          # singularité logarithmique en 0

    # QUAD avec weight='cos' intègre automatiquement
    #   ∫_0^∞ (1/(a+ω)) * cos(ωt) dω
    res, err = quad(lambda w: 1.0 / (a + w), 0, np.inf,
                    weight='cos', wvar=t,
                    limit=500, epsabs=1e-10, epsrel=1e-10)
    return res / np.pi

# Vectorisation pour tracer
f_vec = np.vectorize(f_quad)

# Points de calcul (on évite 0)
t_vals = np.linspace(-5, 5, 1001)
t_vals = t_vals[t_vals != 0]
f_vals = f_vec(t_vals)

# Tracé
plt.figure(figsize=(10, 5))
plt.plot(t_vals, f_vals, label=r'$f(t)$ (quadrature directe)')
plt.xlabel('t')
plt.ylabel('f(t)')
plt.title(r'Transformée inverse de $F(s)=\frac{1}{\ln 2 + |s|}$')
plt.grid(True)
plt.legend()
plt.xlim(t_vals[0], t_vals[-1])
plt.show()

# Quelques valeurs
test_values = np.array([0.1, 0.5, 1.0, 2.0, 3.0])
print("t\t f(t)")
for t in test_values:
    print(f"{t:.2f}\t {f_quad(t):.8f}")
