#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//-----------------------------{Définition des constantes}-----------------------------
#define nbx_particules 10   //Nombre de particules du système
#define nbx_particules_actives 10 //Nombre de particules actives du système (sert juste pour définir des défauts)
#define Rc 2.5 * sigma              //Rayon de coupure 
#define pi M_PI                     //Constante pi
#define L  50               // Taille de la boîte
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
 * @param actife  indique si la particule existe(sert pour les défauts)
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
    p->actif=1;
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



//-----------------------------{Calcul Energie_totale}-----------------------------------------
double Energie_totale(Particule tab_par[]) 
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

    return U+U_tail+energie_cinetique; //retourne l'energie totale
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
 * Le calcul utilise des conditions périodiques aux bords pour gérer les interactions.
 *
 * @param tab_par Tableau contenant les particules du système, chacune avec ses
 *                attributs de position, de vitesse, et de force.
 *
 * @note La fonction suppose que `tab_par` est initialisé avec les particules en positions
 *       valides dans l'espace de simulation.
 */
void calcul_forces(Particule tab_par[]) {

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
 * en fonction des forces appliquées.
 *
 * @param tab_par Tableau de structures `Particule`, contenant les particules du système 
 *                avec leurs positions, vitesses, et forces initiales.
 * @param delta_t Pas de temps pour l'intégration. Un petit pas améliore la précision mais
 *                augmente le temps de calcul.
 *
 * @note Avant d'appeler cette fonction, les forces sur chaque particule doivent avoir été
 *       calculées (par exemple, en utilisant `calcul_forces`), et les valeurs de `fx` et `fy`
 *       doivent être prêtes à l’emploi.
*/
void resolution_verlet(Particule tab_par[], double delta_t) 
{

    //calcul des force à t
    calcul_forces(tab_par);

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
    calcul_forces(tab_par);

    //on calcul ensuite les vitesses en utilisant toujours la méthode de Verlet
    for (int i = 0; i < nbx_particules; i++) 
    {    
        if (!tab_par[i].actif) continue ; //pour verifer qu'on fait pas la resolution pour un defaut si y'en a

        tab_par[i].vx += 0.5 *tab_par[i].fx * delta_t; // vitesse selon x
        tab_par[i].vy += 0.5 * tab_par[i].fy * delta_t; // vitesse selon y
    }
}
void resolution_euler(Particule tab_par[],double delta_t)
{
    //je commence par déterminer f_indice_x :
    calcul_forces(tab_par);
    for ( int i=0;i<nbx_particules;i++)
    {
        tab_par[i].vx = tab_par[i].vx + tab_par[i].fx*delta_t;
        tab_par[i].vy  = tab_par[i].vy + tab_par[i].fy*delta_t;
        tab_par[i].x =  tab_par[i].x + tab_par[i].vx*delta_t;
        tab_par[i].y = tab_par[i].y + tab_par[i].vy*delta_t;

        if(tab_par[i].x > Lmoitie){ tab_par[i].x-=L;}
        else if(tab_par[i].x < -Lmoitie){ tab_par[i].x+=L;}

        if(tab_par[i].y > Lmoitie) {tab_par[i].y-=L;}
        else if(tab_par[i].y < -Lmoitie) {tab_par[i].y+=L;}
    }


}

//methode des moindres carrés pour y=ax+b, n est la dimensuion de y et x (nombre de points)
double mc_ax_plus_b (double x[], double y[], double *a, double *b, double n)
{
    int i;
    double si, sx, sy, sxy, sxx;
    double delta;
    si = sx = sxy = sxx = 0;
    for (i = 0; i <= n - 1; i++)
      {
       if (x[i]>-9 && x[i]<-4)
        {
            si += 1;
            sxx += pow (x[i], 2);
            sxy += x[i] * y[i];
            sx += x[i];
            sy += y[i];
        }
      }

    delta = (sxx * si) - (sx * sx);
    *a = (1 / delta) * (sxy * si - sy * sx);
    *b = (1 / delta) * (sxx * sy - sxy * sx);
}

int main() {
    double max_energie_totale = -1.0; // Variable pour stocker le maximum de l'énergie totale
    double temps_ref = 0.0;

    double delta_t = 5.0e-4;

    Particule tab_particules[nbx_particules];

    init_Particule(&tab_particules[0], 0.369163, 7.001539, 2.461984, 6.004601);

    init_Particule(&tab_particules[1], 11.804992, 12.508632, 6.039470, 6.527302);

    init_Particule(&tab_particules[2], 0.720203, 14.251818, 1.136606, 7.847989);

    init_Particule(&tab_particules[3], 2.790903, 8.519731, 2.597047, 7.523878);

    init_Particule(&tab_particules[4], 5.316689, 9.872263, 0.398034, 5.516499);

    init_Particule(&tab_particules[5], 5.712756, 5.894922, 5.035953, 7.358547);

    init_Particule(&tab_particules[6],8.860042, 7.814479, 6.068676, 9.158643);

    init_Particule(&tab_particules[7], 12.048815, 1.838540, 1.266312, 8.278652);

    init_Particule(&tab_particules[8], 8.840079, 5.592444, 4.283253, 3.763380);

    init_Particule(&tab_particules[9],3.101077, 0.484085, 0.290682, 2.547519);


    double t=0.0,a,b;

    FILE *out;

    int i,j;
    out = fopen ("test1.txt", "w");

    if (!out) {
        perror("Erreur lors de l'ouverture du fichier");
        return EXIT_FAILURE;
    }

    for (j=0;j<3*1/delta_t;j++)
    {
        t+=delta_t;
        resolution_euler(tab_particules,delta_t);
        double energie_totale = Energie_totale(tab_particules);

        fprintf (out, "%22.8g  %22.8g \n",t,energie_totale);

        // Vérifier si delta_t est égal à 5.0e-4 et mettre à jour le maximum de l'énergie totale
        if (energie_totale > max_energie_totale) {
            max_energie_totale = energie_totale;
            temps_ref = t;
        }
    }

    fclose (out);

    if (max_energie_totale != -1.0) {
        printf("Le maximum de l'énergie totale pour delta_t = 5.0e-4 est %f\n et le temps_ref est %f\n", max_energie_totale, temps_ref);
    } else {
        printf("Aucune valeur de delta_t n'était égale à 5.0e-4\n");
    }

    double test[500],diviseur;
    const int dim = sizeof(test) / sizeof(test[0]);
    double x[dim],y[dim];

    for (int i = 0; i <= dim ;i++) {
        diviseur = (i+1)*(i+1)/10;
        test[i] = temps_ref/diviseur;
    }
    
    out= fopen("test.txt", "w");
    double E0=293.71623;
    for (int k=0; k<dim; ++k) 
    {

        delta_t = test[k];
        
        init_Particule(&tab_particules[0], 0.369163, 7.001539, 2.461984, 6.004601);
        init_Particule(&tab_particules[1], 11.804992, 12.508632, 6.039470, 6.527302);
        init_Particule(&tab_particules[2], 0.720203, 14.251818, 1.136606, 7.847989);
        init_Particule(&tab_particules[3], 2.790903, 8.519731, 2.597047, 7.523878);
        init_Particule(&tab_particules[4], 5.316689, 9.872263, 0.398034, 5.516499);
        init_Particule(&tab_particules[5], 5.712756, 5.894922, 5.035953, 7.358547);
        init_Particule(&tab_particules[6], 8.860042, 7.814479, 6.068676, 9.158643);
        init_Particule(&tab_particules[7], 12.048815, 1.838540, 1.266312, 8.278652);
        init_Particule(&tab_particules[8], 8.840079, 5.592444, 4.283253, 3.763380);
        init_Particule(&tab_particules[9], 3.101077, 0.484085, 0.290682, 2.547519);

        double t = 0.0;

        for (int j = 0; j < 3 * 1 / delta_t; j++) 
        {
            t += delta_t;
            resolution_verlet(tab_particules,delta_t);
            if (t >= temps_ref) 
            {
                 max_energie_totale=Energie_totale(tab_particules);
                 break;
                
            }
            if (isnan(log(max_energie_totale-E0))) continue;
        
            fprintf(out, "%22.8g  %22.8g \n",log(delta_t), log(max_energie_totale-E0));
            x[k]=log(delta_t);
            y[k]=log(max_energie_totale-E0);
         } 
    fclose(out); 
    mc_ax_plus_b (x, y,&a,&b,dim);
     FILE *gnuplot = popen("gnuplot -persistent", "w");

    if (gnuplot) {
        // Script Gnuplot pour afficher trois sous-graphiques séparés horizontalement avec des couleurs différentes
        fprintf(gnuplot, "set terminal pngcairo size 800,600\n");  // Largeur augmentée pour chaque graphique
        fprintf(gnuplot, "set output 'exposant.png'\n");
        fprintf(gnuplot, "set xlabel '{/Times-Italic ln(δt)}' font ',18'\n");

        // Sous-graphe pour l'énergie potentielle (bleu)
        fprintf(gnuplot, "set ylabel '{/Times-Italic ln(E(t^*)-E0)}' font ',18'\n");
        //fprintf(gnuplot, "plot 'test.txt' using 1:2 title 'Data from test.txt', 2*x+5.6 with lines title 'y = 2x + 4.6'\n");
        fprintf(gnuplot, "plot 'test.txt' using 1:2 title 'Data', %lf*x+%lf with lines title 'y = %.1lfx + %.1lf'\n", a, b, a, b);

        pclose(gnuplot);
    }
    return 0; 

}