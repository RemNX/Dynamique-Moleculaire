import matplotlib.pyplot as plt
import numpy as np


#Je prends le minimum comme temps de référence t* et ce t* est à 0.00047

#Enregistrer les fichiers au même endroit que le programme py : 
nom_fichiers = []
pas_de_temps=[0.7625000000 
 ,0.3812500000  ,0.2541666667  ,0.1906250000  ,0.1525000000  ,0.1270833333 
 ,0.1089285714  ,0.0953125000  ,0.0847222222  ,0.0762500000  ,0.0381250000 
 ,0.0254166667  ,0.0190625000  ,0.0152500000  ,0.0127083333  ,0.0108928571 
 ,0.0095312500  ,0.0084722222  ,0.0076250000  ,0.0038125000  ,0.0025416667 
 ,0.0019062500  ,0.0015250000  ,0.0012708333  ,0.0010892857  ,0.0009531250 
 ,0.0008472222  ,0.0007625000  ,0.0003812500  ,0.0002541667  ,0.0001906250 
 ,0.0001525000  ,0.0001270833  ,0.0001089286  ,0.0000953125  ,0.0000847222 
 ,0.0000762500  ,0.0000381250  ,0.0000254167  ,0.0000190625  ,0.0000152500 
 ,0.0000127083  ,0.0000108929  ,0.0000095312  ,0.0000084722  ,0.0000076250 
 ,0.0000038125  ,0.0000025417  ,0.0000019062  ,0.0000015250  ,0.0000012708 
 ,0.0000010893  ,0.0000009531  ,0.0000008472]

for i in range(len(pas_de_temps)):
        nom_fichiers.append(f'Energies/Etot_{i}.txt')
       

E0 = 293.67526 
temps_ref_final=0.76250

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
