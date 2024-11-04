import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
nbx_particules = len(positions[1])//2

# Configuration de l'axe de la particule
ax1.set_xlim(-20, 35)
ax1.set_ylim(-20, 35)

lines = [ax1.plot([], [], "ob")[0] for _ in range(nbx_particules)]
# Fonction d'initialisation
def init():
    for line in lines:
        line.set_data([], [])
    return lines 

# Fonction d'animation
def update(i):
    for s in range(nbx_particules):
        pos_t=positions[i]
        pos_x = pos_t[2*s]     # Indice pair pour x
        pos_y = pos_t[2*s + 1] # Indice impair pour y
        # Mise à jour des positions des particules
        lines[s].set_data([pos_x], [pos_y])  # Supposons que y = 0 pour toutes les particules
    return  lines 

# Création de l'animation
ani = FuncAnimation(fig, update, frames=len(time_data), interval=1,init_func=init, blit=True)

# Affichage de l'animation
plt.tight_layout()
plt.show()
