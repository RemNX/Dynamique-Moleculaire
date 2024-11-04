import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random as rd


delta_t = 10e-5
t = np.linspace(0,1,100)
particules=[]
nbx_particules = 2

#Initialisation d'un vecteur avec les positions+vitesse initiale au hasard pour chaque particule :

for k in range(0,nbx_particules):
    if( k==0):
        x_ini = 1.4
        y_ini = 2
        v_x_ini = 6
        v_y_ini = 9
        u_ini = np.array([x_ini, y_ini, v_x_ini, v_y_ini])
        particules.append(u_ini)
    if( k==1):
        x_ini = 1.2
        y_ini = 1.1
        v_x_ini = 2
        v_y_ini = 1
        u_ini = np.array([x_ini, y_ini, v_x_ini, v_y_ini])
        particules.append(u_ini)

#Définition des forces selon x et y :

def force_selon_x(xi,xj,yi,yj):
    rij = np.sqrt(((xi-xj)**2) + ((yi-yj)**2))
    terme_commun = (24*1/(rij**2))*((1/rij)**6)*(2*((1/rij)**6)-1)
    return terme_commun*(xj-xi)

def force_selon_y(xi,xj,yi,yj):
    rij = np.sqrt(((xi-xj)**2) + ((yi-yj)**2))
    terme_commun = (24/(rij**2))*((1/rij)**6)*(2*((1/rij)**6)-1)
    return terme_commun*(yj-yi)
def EM(u1,u2):
    xi,yi,vxi,vyi=u1
    xj,yj,vxj,vyj=u2
    rij = np.sqrt(((xi-xj)**2) + ((yi-yj)**2))
    E=(vxi**2+vyi**2)/2 + (vxj**2+vyj**2)/2
    U=4*((1/rij)**12-(1/rij)**6)
    return U+E

def distance(u1,u2):
    xi,yi,vxi,vyi=u1
    xj,yj,vxj,vyj=u2
    return np.sqrt(((xi-xj)**2) + ((yi-yj)**2))

#résolution :

def resolution(u, indice):
    #je commence par déterminer f_indice_x :
    f_x = 0
    f_y =0
    for l in range(0,nbx_particules):
        if l != indice :
            f_x += force_selon_x(u[indice][0], u[l][0], u[indice][1], u[l][1])
            f_y += force_selon_y(u[indice][0], u[l][0], u[indice][1], u[l][1])
    v_x = u[indice][2] + f_x*delta_t
    v_y = u[indice][3] + f_y*delta_t
    x = u[indice][0] + v_x*delta_t
    y = u[indice][1] + v_y*delta_t
    return np.array([x,y,v_x,v_y])


#fig, ax = plt.subplots()
#lines = [ax.plot([], [], "o")[0] for _ in range(nbx_particules)]
plt.ylim(-1.3, 2)
y_plot=[]
r_plot=[]#np.linspace(1,2.5,100)
#def animate(i): 
for i in range(2000):
    for s in range(nbx_particules):
        u_nouveau = resolution(particules,s)
        particules[s] = u_nouveau
        #lines[s].set_data(particules[s][0], particules[s][1])
    print(EM(particules[0],particules[1]))
    #print(particules[0][2])
    y_plot.append(EM(particules[0],particules[1]))
    r_plot.append(distance(particules[0],particules[1]))
    #print(r_plot)

plt.plot(r_plot,y_plot)


    #return lines
 
#ani = animation.FuncAnimation(fig, animate, frames=1000,interval=100, blit=True, repeat=False)
plt.show()
