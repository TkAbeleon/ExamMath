import sympy as sp

# Symboles
R1, R2, R3, L1, L2, C, s = sp.symbols('R1 R2 R3 L1 L2 C s')
V1 = sp.symbols('V1')
V2, V3 = sp.symbols('V2 V3')

# Admittances / impédances
Y1 = 1/R1 + 1/(L1*s)            # entre nœud 1 et nœud 2
Y_par = 1/R2 + 1/(L2*s) + C*s   # entre nœud 2 et nœud 3

# Lois des nœuds
eq1 = (V2 - V1) * Y1 + (V2 - V3) * Y_par
eq2 = (V3 - V2) * Y_par + V3 / R3

sol = sp.solve([sp.Eq(eq1, 0), sp.Eq(eq2, 0)], [V2, V3], dict=True)[0]
H_num = sp.simplify(sol[V2] / V1)

# Fonction de transfert théorique (d'après les formules complètes ci-dessus)
b3 = C*L1*L2*R2*R3
b2 = C*L2*R1*R2*R3 + L1*L2*R2 + L1*L2*R3
b1 = L1*R2*R3 + L2*R1*R2 + L2*R1*R3
b0 = R1*R2*R3

a3 = C*L1*L2*R1*R2 + C*L1*L2*R2*R3
a2 = C*L2*R1*R2*R3 + L1*L2*R1 + L1*L2*R2 + L1*L2*R3
a1 = L1*R1*R2 + L1*R2*R3 + L2*R1*R2 + L2*R1*R3
a0 = R1*R2*R3

num_th = b3*s**3 + b2*s**2 + b1*s + b0
den_th = a3*s**3 + a2*s**2 + a1*s + a0
H_th = num_th / den_th

# Vérification : la différence doit être nulle
print("Fonction de transfert nodale :")
sp.pprint(sp.simplify(H_num))
print("\nFonction de transfert théorique :")
sp.pprint(sp.simplify(H_th))
print("\nÉgalité ? ", sp.simplify(H_num - H_th) == 0)
