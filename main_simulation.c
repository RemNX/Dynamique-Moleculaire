#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//-----------------------------{Définition des constantes}-----------------------------
#define nbx_particules 100   //Nombre de particules du système
#define nbx_particules_actives 100 //Nombre de particules actives du système (sert juste pour définir des défauts)
#define Rc 2.5 * sigma              //Rayon de coupure 
#define pi M_PI                     //Constante pi
#define L  11               // Taille de la boîte
const int sigma = 1;                // sigma du lennard jones adimentionné
const double Lmoitie = L * 0.5;    // calcul de la valeur de la moitié de L
const double rho = nbx_particules_actives / (L * L); // rho la densité volumique
const double epsilon = 1;          //le E0 du lennard jones adimentionné

/*Calcul pour Rc pour eviter de les refaire, vu qu'ils sont constants*/
const double Rc2 = Rc * Rc;        
const double inv_Rc_2 = 1 / Rc2;
const double inv_Rc_4 = inv_Rc_2 * inv_Rc_2;
const double inv_Rc_6 = inv_Rc_2 * inv_Rc_2 * inv_Rc_2;
const double inv_Rc_10 = inv_Rc_6 * inv_Rc_4;
const double inv_Rc_12 = inv_Rc_6 * inv_Rc_6;

/*Calcul des terme à ajouter pour les grandeurs calculées, qui sont constants*/
const double terme_commun_frac = 2/5;
const double U_tail = pi * epsilon * nbx_particules_actives * rho * (terme_commun_frac * inv_Rc_10 - inv_Rc_4);
const double U_decal  = 4 * (inv_Rc_12 - inv_Rc_6);
const double P_tail = 6*pi*epsilon*rho*rho*(terme_commun_frac*inv_Rc_10 - 0.5*inv_Rc_4);

//-----------------------------{Définition de la structure pour la particule}-----------------------------
/**
 * Définition de la classe Particule avec les attribus.
 *
 * @param x Coordonnée x de la particule.
 * @param y Coordonnée y de la particule
 * @param vx vitesse selon x de la particule.
 * @param vy vitesse selon y de la particule.
 * @param fx force appliquée selon x sur la particule.
 * @param fy force appliquée selon y sur la particule.
 * @param actif  indique si la particule existe(sert pour les défauts)
 */
typedef struct {
    double x; 
    double y;
    double vx;
    double vy;
    double fx;
    double fy;
    int actif; //1=active et 0= inactive, permet d'inclure des défauts
} Particule;

//-----------------------------{Constructeur pour la particule}---------------------------------------
/**
 * @brief Initialise les attributs d'une instance de Particule.
 *
 * Cette fonction configure les coordonnées et les vitesses d'une particule 
 * en initialisant ses attributs avec les valeurs fournies.
 *
 * @param p Pointeur vers l'instance de Particule à initialiser.
 * @param x Coordonnée x de la particule.
 * @param y Coordonnée y de la particule.
 * @param vx Vitesse en x de la particule.
 * @param vy Vitesse en y de la particule.
 */
void init_Particule(Particule *p, double x, double y, double vx, double vy) {
    // Initialisation des coordonnées de la particule
    p->x = x;
    p->y = y;
    // Initialisation des vitesses de la particule
    p->vx = vx;
    p->vy = vy;
}

//-----------------------------{Initialise la configuration}------------------------------------------
/**
 * @brief Initialise une configuration cristalline pour un ensemble de particules.
 *
 * Cette fonction dispose les particules dans une structure cristalline,
 * typiquement en utilisant une grille régulière dans l'espace. 
 * Cela signifie que chaque particule est placée à une position fixe et ordonnée.
 *
 * @param tab_par Pointeur vers un tableau de structures `Particule`, 
 *                représentant les particules du système à initialiser.
 *                Le tableau doit être préalablement alloué avec le nombre
 *                de particules souhaité.
 *
 * @note Cette fonction suppose que le tableau `tab_par` est suffisamment 
 *       grand pour contenir toutes les particules nécessaires.
 * 
 */
void initialiserConfigurationCristalline(Particule *tab_par) {
    // Calculer le nombre de particules par ligne et par colonne
    int particules_par_ligne = (int)sqrt(nbx_particules); 
    int particules_par_colonne = (nbx_particules + particules_par_ligne - 1) / particules_par_ligne; 
    
    // Espacement entre les particules
    double espacement_x = 1.1; // Ajustez pour obtenir l'espacement souhaité
    double espacement_y = 1.1;

    // Calcul de la largeur et hauteur du réseau
    double largeur_reseau = espacement_x * (particules_par_ligne - 1);
    double hauteur_reseau = espacement_y * (particules_par_colonne - 1);

    // Décalage pour centrer le réseau autour de (0, 0)
    double decalage_x = -largeur_reseau / 2.0;
    double decalage_y = -hauteur_reseau / 2.0;

    for (int i = 0; i <= particules_par_ligne; i++) 
    {
        for (int j = 0; j < particules_par_colonne; j++) 
        {
            int index = i * particules_par_colonne + j;
            if (index < nbx_particules) 
            {
                // Position initiale centrée autour de (0, 0)
                tab_par[index].x = i * espacement_x + decalage_x;
                tab_par[index].y = j * espacement_y + decalage_y;
                //tab_par[index].vx = 0; 
                //tab_par[index].vy = 0;
                tab_par[index].fx=0;
                tab_par[index].fy=0;

                // Générer une vitesse initiale entre -0.5 et 0.5, puis la multiplier par le facteur
                tab_par[index].vx = ((double)rand() / RAND_MAX - 0.5) * 10;
                tab_par[index].vy = ((double)rand() / RAND_MAX - 0.5) * 10;

                tab_par[index].actif = 1; 
                //si jamais on veut ajouter des defauts 
                /*if(index==12 || index==85)
                {
                    tab_par[index].actif=0; 
                }*/
            }
        }
    }
}



//-----------------------------{Calcul de la température}-----------------------------------------
/**
 * @brief Calcule les grandeurs thermodynamiques du système.
 *
 * Cette fonction effectue les calculs pour déterminer la température, 
 * la pression, l'énergie potentielle et l'énergie totale du système.
 *
 * @param tab_par Tableau des particules du système.
 * @param temperature Pointeur vers la variable où sera stockée la température calculée.
 * @param pression Pointeur vers la variable où sera stockée la pression calculée.
 * @param energie_potentielle Pointeur vers la variable où sera stockée l'énergie potentielle calculée.
 * @param energie_totale Pointeur vers la variable où sera stockée l'énergie totale calculée.
 * @param viriel Valeur du viriel du système, qui représente la somme des forces interparticulaires pondérées.
 */
void calcul_grandeurs(Particule tab_par[], double *temperature , double *pression , double *enegrie_potentielle, double *enegrie_totale,double viriel) 
{
    /*définition des variables locales à utiliser dans les calculs*/
    double energie_cinetique = 0.0, U = 0.0; // energie cinétique et potentielle
    double rij2, inv_r2, inv_r2x3, inv_r2x6; //puissance de rij
    double V = L * L; //volume
    int dimension = 2;

    for (int i = 0; i < nbx_particules; i++) 
    {
        /*condition qui verifie que le calcul se fait que par rapports aux partiucules et ignore les défauts*/
        if (tab_par[i].actif) 
        {
            energie_cinetique += (tab_par[i].vx * tab_par[i].vx + tab_par[i].vy * tab_par[i].vy) / 2.0; //calcul de l'energie cinétique
            for (int j = i + 1; j < nbx_particules; j++) 
            {
                if (!tab_par[j].actif) continue; 

                double x_sous = tab_par[j].x - tab_par[i].x; //Dx
                double y_sous = tab_par[j].y - tab_par[i].y; //
                /*Calcul de la distance entre deux particules*/
                rij2 = x_sous * x_sous + y_sous * y_sous;   

                /*Calcul sur les puissances de rij pour eviter l'appel de pow*/
                inv_r2 = 1.0 / rij2;                  //   1/rij^2
                inv_r2x3 = inv_r2 * inv_r2 * inv_r2; //    1/rji^6
                inv_r2x6 = inv_r2x3 * inv_r2x3;     //    1/rij^12

                /*Conditions qui permet de tronquer le potentiel à Rc*/
                if (rij2 < Rc2) 
                {
                    U += 4 * (inv_r2x6 - inv_r2x3) - U_decal; //calcul du potentiel entre deux particules
                }
            }
        }
        
    }
    *temperature = (energie_cinetique / nbx_particules_actives); //calcul et mise à jour de la temperature
    *pression = (nbx_particules_actives * (*temperature) / V) + (viriel / (dimension * V*nbx_particules_actives))+P_tail; //calcul et mise à jour de la pression
    *enegrie_potentielle = U+U_tail;    //calcul et mise à jour l'energie potentielle
    *enegrie_totale=(*enegrie_potentielle)+energie_cinetique; //calcul et mise à jour de l'energie totale
}
//-----------------------------{Enregistrer les positions}-------------------------------------------
void enregistrer_Positions(Particule *tab_par,double t, FILE *file) {
    fprintf (file, "%22.8g",t);
        for ( int i=0;i<nbx_particules;i++)
        {

        if (!tab_par[i].actif) continue;

         fprintf (file, "%22.8g %22.8g",tab_par[i].x,tab_par[i].y);
        }
        fprintf (file, "\n");
}
//-----------------------------{Calcul_forces}--------------------------------------------------------
/**
 * @brief Calcule les forces entre toutes les paires de particules dans le système.
 *
 * Cette fonction parcourt le tableau des particules et calcule les forces d'interaction
 * entre toutes les paires, mettant ainsi à jour les composantes de force de chaque particule.
 * Le calcul utilise des conditions périodiques aux bords pour gérer les interactions
 * et contribue également au calcul du viriel, utilisé dans le calcul de la pression.
 *
 * @param tab_par Tableau contenant les particules du système, chacune avec ses
 *                attributs de position, de vitesse, et de force.
 * @param viriel Pointeur vers une variable où le viriel total du système sera stocké.
 *               Le viriel est une somme pondérée des forces entre particules et sert
 *               dans les calculs thermodynamiques (notamment la pression).
 *
 * @note La fonction suppose que `tab_par` est initialisé avec les particules en positions
 *       valides dans l'espace de simulation.
 */
void calcul_forces(Particule tab_par[], double *viriel) {
    // Initialiser le viriel (sert pour le calcul de la pression)
    *viriel = 0.0;

    // Réinitialiser les forces et l'énergie potentielle
    for (int i = 0; i < nbx_particules; i++) {
        tab_par[i].fx = 0.0;
        tab_par[i].fy = 0.0;
    }

    // Calcul des forces d'interaction entre les particules
    for (int i = 0; i < nbx_particules; i++) {
        if (!tab_par[i].actif) continue;

        for (int j = i + 1; j < nbx_particules; j++) {
            if (!tab_par[j].actif) continue;

            double x_sous = tab_par[i].x - tab_par[j].x;
            double y_sous = tab_par[i].y - tab_par[j].y;

            double abs_distx = fabs(x_sous);
            double abs_disty = fabs(y_sous);

            if (abs_distx > Lmoitie) {
                x_sous = - (x_sous / abs_distx) * (L - abs_distx);
            }

            if (abs_disty > Lmoitie) {
                y_sous = - (y_sous / abs_disty) * (L - abs_disty);
            }

            double rij2 = x_sous * x_sous + y_sous * y_sous;

            if (rij2 < Rc2) {
                double r6 = rij2 * rij2 * rij2;
                double r12 = r6 * r6;
                double inv_r2 = 1.0 / rij2;
                double inv_r6 = 1.0 / r6;
                double inv_r12 = 1.0 / r12;
                double terme_commun = 24.0 * (2 * inv_r12 - inv_r6) * inv_r2;

                double force_x = terme_commun * x_sous;
                double force_y = terme_commun * y_sous;

                tab_par[i].fx += force_x;
                tab_par[j].fx -= force_x;
                tab_par[i].fy += force_y;
                tab_par[j].fy -= force_y;

                *viriel += terme_commun * rij2;
            }
        }
    }
}


//-----------------------------{resolution}---------------------------------------------------
/**
 * @brief Effectue l'intégration des équations de mouvement des particules.
 *
 * Cette fonction met à jour les positions et vitesses des particules en appliquant un 
 * schéma d'intégration de Vertlet pour simuler leur évolution temporelle
 * en fonction des forces appliquées. Le viriel est passé comme argument dans calcul_forces()
 * pour le mettre à jour.
 *
 * @param tab_par Tableau de structures `Particule`, contenant les particules du système 
 *                avec leurs positions, vitesses, et forces initiales.
 * @param delta_t Pas de temps pour l'intégration. Un petit pas améliore la précision mais
 *                augmente le temps de calcul.
 * @param viriel Pointeur vers le viriel calculé, qui sera mis à jour pendant l'intégration.
 *
 * @note Avant d'appeler cette fonction, les forces sur chaque particule doivent avoir été
 *       calculées (par exemple, en utilisant `calcul_forces`), et les valeurs de `fx` et `fy`
 *       doivent être prêtes à l’emploi.
*/
void resolution(Particule tab_par[], double delta_t, double *viriel) 
{

    //calcul des force à t
    calcul_forces(tab_par, viriel);

    for (int i = 0; i < nbx_particules; i++) 
    {
        if (!tab_par[i].actif) continue; //pour verifer qu'on fait pas la resolution pour un defaut si y'en a

        // on calcul la nouvelle position en x et y en utilisant la méthode de Verlet
        tab_par[i].x += tab_par[i].vx * delta_t + 0.5 * tab_par[i].fx * (delta_t * delta_t);
        tab_par[i].y += tab_par[i].vy * delta_t + 0.5 * tab_par[i].fy * (delta_t * delta_t);

        /*puis on pose les conditions periodiques qui font que si la particule sort de la boite, elle 
        rentre de l'autre coté*/
        //selon x
        if (tab_par[i].x > Lmoitie) {tab_par[i].x -= L;}
        else if (tab_par[i].x < -Lmoitie) {tab_par[i].x += L;}
        //selon y
        if (tab_par[i].y > Lmoitie) {tab_par[i].y -= L;}
        else if (tab_par[i].y < -Lmoitie) {tab_par[i].y += L; }

        tab_par[i].vx += 0.5*tab_par[i].fx*delta_t;
        tab_par[i].vy += 0.5*tab_par[i].fy*delta_t;
    }
    //calcul des force à t+dt
    calcul_forces(tab_par,viriel);

    //on calcul ensuite les vitesses en utilisant toujours la méthode de Verlet
    for (int i = 0; i < nbx_particules; i++) 
    {    
        if (!tab_par[i].actif) continue ; //pour verifer qu'on fait pas la resolution pour un defaut si y'en a

        tab_par[i].vx += 0.5 *tab_par[i].fx * delta_t; // vitesse selon x
        tab_par[i].vy += 0.5 * tab_par[i].fy * delta_t; // vitesse selon y
    }
}

//-----------------------------{main}----------------------------------------------------------------
int main() 
{
    //srand(42);

    const double delta_t = 5.0e-4; //definition du pas de temps 
    double viriel = 0.0,pression,energie_potentielle,temperature,energie_totale; //definition des variables pour le calcul grandeurs thermodynamiques

    /*Pour le calcul des valeurs moyennes*/
    double moyenne_temp=0;
    double moyenne_press=0;
    double moyenne_ep=0; 
    double count=0;

    double t = 0.0; //definition de la variable dont on va stocker le temps incrementé
    // Creation d'un tableau de particules
    Particule tab_particules[nbx_particules];
    // Initialisation des particules avec une configuration cristalline
    initialiserConfigurationCristalline(tab_particules);   
    
    // Ouvre le fichier "positions_data.txt" en mode écriture ("w")
    // Ce fichier stockera les positions des particules au fil du temps dans le format suivant:
    // temps   x1 y1    x2 y2    x3 y3    x4 y4.....
    FILE *pos_file = fopen("positions_data.txt", "w");

    // Ouvre le fichier "simulation_data.txt" en mode écriture ("w")
    // Ce fichier stockera  : temps, energie_totale, energie_potentielle, temperature, pression.
    FILE *data_file = fopen("simulation_data.txt", "w");

    // Appelle la fonction enregistrer_Positions pour écrire les positions initiales des particules
    enregistrer_Positions(tab_particules, t, pos_file);

    for (int j = 0; j < 3*1/delta_t ; j++) 
    {
        t += delta_t; // Incrémenter le temps de simulation de delta_t

        /* Appel de la fonction resolution pour resoudre les equations de mouvement et mettre 
        à jour les positions et les vitesses des particules ,et '&viriel' passe l'adresse de 
        la variable 'viriel' (par référence), permettant à la fonction resolution de modifier sa valeur directement*/
        resolution(tab_particules, delta_t, &viriel); 

        // Calculer les grandeurs thermodynamiques
        calcul_grandeurs(tab_particules,&temperature,&pression,&energie_potentielle,&energie_totale, viriel);
        if (t>0.5){
            moyenne_ep += energie_potentielle;
            moyenne_press += pression;
            moyenne_temp += temperature;
            count +=1;
        }
        // Enregistrer les positions de chaque pas de temps
        enregistrer_Positions(tab_particules,t,pos_file);         
        // Écrire l'énergie totale, la température et la pression dans le fichier
        fprintf(data_file, "%22.8g  %22.8g  %22.8g  %22.8g  %22.8g\n", t, energie_totale, energie_potentielle, temperature, pression);
    }
    fclose(data_file); // Ferme le fichier de simulation_data.txt pour s'assurer que toutes les données sont écrites et libérer les ressources associées
    fclose(pos_file);  // De même pour le fichier positions_data.txt

    moyenne_ep = moyenne_ep/count;
    moyenne_press = moyenne_press/count;
    moyenne_temp = moyenne_temp/count;

    printf("Moyenne de l'énergie potentielle: %.6f\n", moyenne_ep);
    printf("Moyenne de la température: %.6f\n", moyenne_temp);
    printf("Moyenne de la pression: %.6f\n", moyenne_press);


    // Ouvre un processus Gnuplot en mode écriture ("w") pour envoyer des commandes de tracé
    // "gnuplot -persistent" permet de garder la fenêtre du graphique ouverte après l'exécution des commandes
    FILE *gnuplot = popen("gnuplot -persistent", "w");

    if (gnuplot) {
    // Script Gnuplot pour afficher trois sous-graphiques séparés horizontalement avec des couleurs différentes
    fprintf(gnuplot, "set terminal pngcairo size 1800,600\n");  // Largeur augmentée pour chaque graphique
    fprintf(gnuplot, "set output 'simulation_graphs.png'\n");
    fprintf(gnuplot, "set xlabel '{/Times-Italic t}' font ',18'\n");

    // Configuration de multiplot pour trois sous-graphiques alignés horizontalement
    fprintf(gnuplot, "set multiplot layout 1, 3 \n");
    fprintf(gnuplot, "unset key \n");

    // Sous-graphe pour l'énergie potentielle (bleu)
    fprintf(gnuplot, "set ylabel '{/Times-Italic U}' font ',18'\n");
    fprintf(gnuplot, "set title '{/Times-Italic Énergie Potentielle}' font ',20' \n");
    fprintf(gnuplot, "plot 'simulation_data.txt' using 1:3 with lines linecolor rgb '#1976D2' title 'Énergie Potentielle', \\\n");
    fprintf(gnuplot, "     (%.6f) title 'Moyenne' with lines linestyle 2 linecolor rgb '#FF9800'\n", moyenne_ep);

    // Sous-graphe pour la température (rouge)
    fprintf(gnuplot, "set ylabel '{/Times-Italic T}' font ',18'\n");
    fprintf(gnuplot, "set title '{/Times-Italic Température}' font ',20'\n"); 
    fprintf(gnuplot, "plot 'simulation_data.txt' using 1:4 with lines linecolor rgb '#D32F2F' title 'Température', \\\n");
    fprintf(gnuplot, "     (%.6f) title 'Moyenne' with lines linestyle 2 linecolor rgb '#FF9800'\n", moyenne_temp);

    // Sous-graphe pour la pression (vert)
    fprintf(gnuplot, "set ylabel '{/Times-Italic P}' font ',18'\n");
    fprintf(gnuplot, "set title '{/Times-Italic Pression}' font ',20' \n");
    fprintf(gnuplot, "plot 'simulation_data.txt' using 1:5 with lines linecolor rgb '#4CAF50' title 'Pression', \\\n");
    fprintf(gnuplot, "     (%.6f) title 'Moyenne' with lines linestyle 2 linecolor rgb '#FF9800'\n", moyenne_press);

    // Fin de multiplot
    fprintf(gnuplot, "unset multiplot\n");

    // Fermeture du processus Gnuplot
    pclose(gnuplot);
}

    return 0;
}
