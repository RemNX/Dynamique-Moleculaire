import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Définition de L pour la boîte de délimitation
L = 11  # Définissez la longueur du côté souhaitée pour la boîte

# Chargement des données de température
time_data = []
# Chargement des données de positions
positions = []

with open('positions_data.txt', 'r') as pos_file:
    for line in pos_file:
        data = list(map(float, line.split()))
        time_data.append(data[0])
        positions.append(data[1:])  # Les autres éléments sont les positions

# Création de la figure et des axes
fig, ax1 = plt.subplots()
nbx_particules = len(positions[1]) // 2

# Configuration des limites de l'axe pour correspondre à la taille de la boîte
ax1.set_xlim(-L/2 - 5, L/2 + 5)
ax1.set_ylim(-L/2 - 5, L/2 + 5)

# Dessin de la boîte de délimitation avec un contour plus épais
boite = plt.Rectangle((-L/2, -L/2), L, L, fill=False, color="red", linestyle="--", linewidth=2)
ax1.add_patch(boite)

# Initialisation des particules
lines = [ax1.plot([], [], "ob")[0] for _ in range(nbx_particules)]

# Ajout du texte pour le compteur de temps
time_text = ax1.text(0.95, 0.95, '', transform=ax1.transAxes, ha='right', va='top', fontsize=12, color="black")

for s in range(nbx_particules):
    pos_x = positions[0][2 * s]     # Position x initiale de la particule
    pos_y = positions[0][2 * s + 1] # Position y initiale de la particule
    lines[s].set_data([pos_x], [pos_y])

# Sauvegarde de l'image initiale
fig.set_size_inches(8, 6)  # Définit la taille en pouces pour obtenir 800x600 pixels à 100 DPI
plt.savefig("etat_initial.png", dpi=300)  # Enregistrer avec une résolution de 300 DPI

# Fonction d'initialisation
def init():
    for line in lines:
        line.set_data([], [])
    time_text.set_text('')  # Texte initial vide pour le compteur de temps
    return lines + [time_text]

# Fonction d'animation
def update(i):
    # Mise à jour des positions des particules
    for s in range(nbx_particules):
        pos_t = positions[i]
        pos_x = pos_t[2 * s]     # Indice pair pour x
        pos_y = pos_t[2 * s + 1] # Indice impair pour y
        lines[s].set_data([pos_x], [pos_y])

    # Mise à jour du compteur de temps
    time_text.set_text(f'Temps: {time_data[i]:.2f}s')
    return lines + [time_text]

# Création de l'animation
ani = FuncAnimation(fig, update, frames=range(0,len(time_data),10), interval=1, init_func=init, blit=True)

# Affichage de l'animation
plt.tight_layout()
plt.show()
