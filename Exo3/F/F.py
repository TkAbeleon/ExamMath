import numpy as np
import matplotlib.pyplot as plt

# Domaine : demi-plan Re(s) > -2
x_min, x_max = -1.9, 5.0
y_min, y_max = -4.0, 4.0
nx, ny = 201, 201

x = np.linspace(x_min, x_max, nx)
y = np.linspace(y_min, y_max, ny)
X, Y = np.meshgrid(x, y)
S = X + 1j*Y

# Calcul de F(s) = 1 / ln(2 + |s|)
# |s| = sqrt(x^2 + y^2)
modS = np.abs(S)
F = 1.0 / np.log(2.0 + modS)

# Partie réelle et imaginaire
u = np.real(F)
v = np.imag(F)

# Vérification que v est partout nulle (aux erreurs numériques près)
print(f"Max |Im(F)| = {np.max(np.abs(v)):.3e}")

# Calcul des dérivées partielles par différences finies
dx = x[1] - x[0]
dy = y[1] - y[0]

# du/dx, du/dy
dudx = (u[2:,1:-1] - u[:-2,1:-1]) / (2*dx)
dudy = (u[1:-1,2:] - u[1:-1,:-2]) / (2*dy)

# dv/dx, dv/dy (devraient être nuls mais calculés quand même)
dvdx = (v[2:,1:-1] - v[:-2,1:-1]) / (2*dx)
dvdy = (v[1:-1,2:] - v[1:-1,:-2]) / (2*dy)

# Cauchy-Riemann : du/dx == dv/dy  et  du/dy == -dv/dx
erreur_CR1 = np.abs(dudx - dvdy)
erreur_CR2 = np.abs(dudy + dvdx)
print(f"Erreur max Cauchy-Riemann (1) : {np.max(erreur_CR1):.3e}")
print(f"Erreur max Cauchy-Riemann (2) : {np.max(erreur_CR2):.3e}")
print("Si ces erreurs ne sont pas ~0, la fonction n'est pas holomorphe.")
print("Ici elles sont grandes car la fonction est réelle et dépend de |s|.")

# Visualisation : partie réelle u(x,y) et sa dépendance en |s|
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
cp1 = ax1.contourf(X, Y, u, levels=30, cmap='viridis')
plt.colorbar(cp1, ax=ax1)
ax1.set_title(r"$u(x,y) = \mathrm{Re}\,F(s)$")
ax1.set_xlabel('Re(s)')
ax1.set_ylabel('Im(s)')

# Coupe pour y=0
ax2.plot(x, u[:, (ny-1)//2], 'r-', label='coupe Im(s)=0')
ax2.set_title('F(s) sur l\'axe réel')
ax2.set_xlabel('s (réel)')
ax2.set_ylabel('F(s)')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig('preuve_non_holomorphe.png', dpi=150)
plt.show()
