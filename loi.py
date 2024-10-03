import matplotlib.pyplot as plt
import numpy as np


#Je prends le minimum comme temps de référence t* et ce t* est à 0.00047

#Enregistrer les fichiers au même endroit que le programme py : 
nom_fichiers = []
pas_de_temps=[]

for i in range(-5,-1):
    for j in np.arange(1,9.1,0.4):
        j_rounded = round(j, 1)
        nom_fichiers.append(f'Etot_{j_rounded}e{i}.txt')
        pas_de_temps.append(float(f'{j_rounded}e{i}'))

data_ref=[]
temps_ref=[]

with open('Etot_5.0e-4.txt','r') as file :
    for line in file :
        left,right = line.split(':')
        data_ref.append(float(right.strip()))
        temps_ref.append(float(left.strip()))
    maxi = max(data_ref)
    for h in range(0,len(data_ref)):
        if data_ref[h] == maxi :
            indice=h

E0 = data_ref[0]
temps_ref_final=temps_ref[indice]

liste_temps=[]
data = []
liste_E_t_ref=[]

for i in nom_fichiers :
    with open(i, 'r') as file :

        for line in file :
            left,right = line.split(':')
            data.append(float(right.strip()))
            liste_temps.append(float(left.strip()))

        init = abs(liste_temps[0]-temps_ref_final)
        for l in range(0,len(liste_temps)):
            if abs(liste_temps[l]-temps_ref_final) <= init :
                closest_index=l

        closest_value = data[closest_index]
        print(closest_value)
        liste_E_t_ref.append(closest_value)

        data=[]
        liste_temps=[]

liste_diff = [np.log(abs(i - E0)) for i in liste_E_t_ref]
liste_diff_deltat =[np.log(j) for j in pas_de_temps]
droite = [k-1.6 for k in liste_diff_deltat]

plt.plot( liste_diff_deltat, liste_diff, "o")
plt.plot(liste_diff_deltat, droite)
plt.xlabel("ln(δt)")
plt.ylabel("ln(|ΔE(t*)|)")
#plt.xlim(-9.5,-7)
plt.show()
