import matplotlib.pyplot as plt
import numpy as np


#Je prends le minimum comme temps de référence t* et ce t* est à 0.00047

#Enregistrer les fichiers au même endroit que le programme py : 
nom_fichiers = []
pas_de_temps=[0.007625  ,0.003812  ,0.002542  ,0.001906  ,0.001525  ,0.001271  ,0.001089  ,0.000953  ,0.000847  ,0.000762  ,0.000693  ,0.000635  ,0.000587  ,0.000545  ,0.000508  ,0.000477  ,0.000449  ,0.000424  ,0.000401  ,0.000381  ,0.000363  ,0.000347  ,0.000332  ,0.000318  ,0.000305  ,0.000293  ,0.000282  ,0.000272  ,0.000263  ,0.000254]

for i in range(30):
        nom_fichiers.append(f'Etot_{i}.txt')
       

E0 = 293.67526 
temps_ref_final=0.762500

liste_temps=[]
data = []
liste_E_t_ref=[]

for i in nom_fichiers :
    with open(i, 'r') as file :

        for line in file :
            left,right = line.split(':')
            data.append(float(right.strip()))
            liste_temps.append(float(left.strip()))

        for l in range(0,len(liste_temps)):
            if abs(liste_temps[l]-temps_ref_final) <= 0 :
                closest_index=l

        closest_value = data[closest_index]
        print(closest_value)
        liste_E_t_ref.append(closest_value)

        data=[]
        liste_temps=[]

liste_diff = [np.log(abs(i - E0)) for i in liste_E_t_ref]
liste_diff_deltat =[np.log(j) for j in pas_de_temps]
droite = [k+4.6 for k in liste_diff_deltat]

plt.plot( liste_diff_deltat, liste_diff, "o")
plt.plot(liste_diff_deltat, droite)
plt.xlabel("ln(δt)")
plt.ylabel("ln(|ΔE(t*)|)")
#plt.xlim(-9.5,-7)
plt.show()
