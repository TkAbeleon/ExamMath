import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as plt

def verify_exo2():
    # Définition du domaine pour s, on commence un peu après la singularité
    s_min = np.exp(-np.exp(-1)) - 1 + 1e-4
    s_vals = np.linspace(s_min, 5.0, 1000)
    
    # 1. Calcul analytique de F(s)
    ln_s1 = np.log(s_vals + 1)
    F_vals = np.real(lambertw(ln_s1))
    
    # 2. Vérification de l'équation différentielle (F'(s))
    # Formule théorique 2 : F'(s) = F(s) / ((s+1) * ln(s+1) * (1+F(s)))
    # On évite s=0 où ln(1)=0 et F(0)=0 (la limite vaut 1)
    mask = np.abs(s_vals) > 1e-4
    F_prime_theo = np.ones_like(s_vals) # Initialise à 1 (pour s=0)
    F_prime_theo[mask] = F_vals[mask] / ((s_vals[mask] + 1) * ln_s1[mask] * (1 + F_vals[mask]))
    
    # Dérivée numérique
    ds = s_vals[1] - s_vals[0]
    F_prime_num = np.gradient(F_vals, ds)
    
    # Pour comparer la dérivée numérique, on ignore les bords et on s'éloigne de la singularité (-0.37)
    valid_range = (s_vals > -0.3) & (s_vals < 4.9)
    err_diff = np.max(np.abs(F_prime_num[valid_range] - F_prime_theo[valid_range]))
    print(f"Erreur max entre dérivée numérique et dérivée théorique (pour s > -0.3) : {err_diff:.2e}")

    # ==========================
    # GÉnÉRATION DE LA COURBE
    # ==========================
    plt.figure(figsize=(12, 5))
    
    # Courbe 1 : F(s)
    plt.subplot(1, 2, 1)
    plt.plot(s_vals, F_vals, 'b-', lw=2, label=r'$F(s) = W(\ln(s+1))$')
    plt.xlabel('s')
    plt.ylabel('F(s)')
    plt.title("Solution F(s)")
    plt.grid(True, linestyle='--')
    plt.legend()
    
    # Courbe 2 : F'(s) théorique vs numérique
    plt.subplot(1, 2, 2)
    plt.plot(s_vals, F_prime_theo, 'g-', lw=3, label=r"$F'(s)$ analytique")
    plt.plot(s_vals, F_prime_num, 'r--', lw=2, label=r"$F'(s)$ numérique")
    
    # On restreint l'axe y car la dérivée tend vers l'infini près de s = -0.37
    plt.ylim(-1, 5)
    plt.xlabel('s')
    plt.ylabel("F'(s)")
    plt.title("Vérification de l'équation différentielle")
    plt.grid(True, linestyle='--')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('Exo2_verif_plot.png', dpi=300)
    print("Graphique de vérification sauvegardé sous : Exo2_verif_plot.png")

if __name__ == '__main__':
    verify_exo2()
