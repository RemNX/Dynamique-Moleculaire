import matplotlib.pyplot as plt
import numpy as np


#Je prends le minimum comme temps de référence t* et ce t* est à 0.00047

#Enregistrer les fichiers au même endroit que le programme py : 
nom_fichiers = ['Etot17.txt', 'Etot16.txt', 'Etot15.txt', 'Etot14.txt', 'Etot13.txt', 'Etot12.txt', 'Etot11.txt', 'Etot10.txt', 'Etot.txt', 'Etot2.txt', 'Etot3.txt', 'Etot4.txt', 'Etot5.txt', 'Etot6.txt', 'Etot7.txt', 'Etot8.txt', 'Etot9.txt', 'Etot18.txt', 
                'Etot19.txt', 'Etot20.txt', 'Etot21.txt', 'Etot22.txt', 'Etot23.txt', 'Etot24.txt', 'Etot25.txt',  'Etot26.txt',  'Etot27.txt']
liste_temps=[]
pas_de_temps=[10e-14,10e-13,10e-12,10e-11,10e-10,10e-9,10e-8,10e-7,10e-6, 10e-5, 10e-4, 10e-3, 10e-2, 10e-1, 10e0, 10e1, 10e2, 10e3, 10e4, 10e5, 10e6, 10e7, 10e8, 10e9, 10e10,10e20, 10e30]
data = []
liste_E_t_ref=[]

#récupéré avec un min pour un pas de 10e-5 :
temp_ref = 0.00047 
E0=64559.103

for i in nom_fichiers :
    with open(i, 'r') as file :

        for line in file :
            left,right = line.split(':')
            data.append(float(right.strip()))
            liste_temps.append(float(left.strip()))

        closest_index = min(range(len(liste_temps)), key=lambda i: abs(liste_temps[i] - temp_ref))
        closest_value = data[closest_index]
        liste_E_t_ref.append(closest_value)

        data=[]
        liste_temps=[]

liste_diff = [np.log(abs(i - E0)) for i in liste_E_t_ref]
liste_diff_deltat =[np.log(j) for j in pas_de_temps]
droite=[k for k in liste_diff_deltat]


plt.plot( liste_diff_deltat, liste_diff, "o")
plt.plot(liste_diff_deltat, droite)
plt.xlabel("ln(δt)")
plt.ylabel("ln(|ΔE(t*)|)")
plt.show()