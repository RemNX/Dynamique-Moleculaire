#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//-----------------------------{Définition des constantes}-----------------------------

#define nbx_particules   10
const int sigma=1;
#define Rc   3.0*sigma
#define pi   M_PI
#define MIN_DISTANCE 2.0
#define L 100 //Taille de la boîte 
const double Lmoitie = L*0.5;
const double rho = nbx_particules/(L*L);
const double epsilon = 1;

const double Rc2 = Rc*Rc;
const double inv_Rc_2 = 1./Rc2;
const double inv_Rc_4 = inv_Rc_2*inv_Rc_2;
const double inv_Rc_6 = inv_Rc_2*inv_Rc_2*inv_Rc_2;
const double inv_Rc_10 = inv_Rc_6*inv_Rc_4;
const double inv_Rc_12 = inv_Rc_6*inv_Rc_6;



//Si jamais à un moment on veut que sigma soit différent de 1 : 
//sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
//sigma12 = sigma6*sigma6;

const double U_tail = pi*epsilon*nbx_particules*rho*((2./5)*inv_Rc_10 - inv_Rc_4);

//-----------------------------{Définition de la structure pour la particule}-----------------------------

/*Definition Particule avec x,y,vitesse selon x  
, vitesse selon y*/ 
typedef struct
{
    double x;
    double y;
    double vx;
    double vy;
    double fx;
    double fy;
} Particule;

/*Constructeur*/
void initPoint(Particule *p, double x, double y, double vx, double vy)
{
    p->x = x;
    p->y = y;
    p->vx = vx;
    p->vy = vy;
}

//-----------------------------{Genere la particule}-----------------------------

Particule genererParticule(Particule *existingParticles, int count) {
    Particule p;
    int fac = 20; // Facteur maximum pour les valeurs aléatoires
    int isTooClose;
    
    do {
        isTooClose = 0; // Réinitialiser le flag
        p.x = (double)rand() / RAND_MAX * fac;
        p.y = (double)rand() / RAND_MAX * fac;
        p.vx = (double)rand() / RAND_MAX * 10;
        p.vy = (double)rand() / RAND_MAX * 10;

        // Vérifier la distance avec chaque particule existante
        for (int i = 0; i < count; i++) {
            double distance = sqrt(pow(p.x - existingParticles[i].x, 2) + pow(p.y - existingParticles[i].y, 2));
            if (distance < MIN_DISTANCE) {
                isTooClose = 1; // Une particule est trop proche
                break; // Pas besoin de vérifier plus
            }
        }
    } while (isTooClose); // Regénérer si trop proche

    return p;
}

//-----------------------------{Calcul de l'énergie totale}-----------------------------

double Energie_totale(Particule tab_par[])
{
    double E=0.0,U=0.0;
    double rij2, inv_r2, inv_r2x3, inv_r2x6;

    for (int i=0;i<nbx_particules;i++)
    {
        for (int j=i+1;j<nbx_particules;j++)
        {

            double x_sous = tab_par[j].x - tab_par[i].x;
            double y_sous = tab_par[j].y - tab_par[i].y;
            rij2 = x_sous*x_sous + y_sous*y_sous;
            inv_r2 = 1.0 / rij2;
            inv_r2x3 = inv_r2 * inv_r2 * inv_r2;
            inv_r2x6 = inv_r2x3 * inv_r2x3;
            if (rij2<Rc2){
                U+=4*(inv_r2x6-inv_r2x3) - 4*(inv_Rc_12 - inv_Rc_6);
            }

        }
        E+=(tab_par[i].vx*tab_par[i].vx + tab_par[i].vy*tab_par[i].vy)/2.; //Energie cinétique

    }

    return U+E+U_tail;
}

//-----------------------------{Calcul_forces}-----------------------------

/*Fonction qui resoud les equations differentielles d'une particule
designé par un "indice", et qui prend un tableau dans lequel sont
stockées les particules, elle ne retourne rien car tout est modifié
directement */
void calcul_forces (Particule tab_par[])
{
    double f_x, f_y;
    double rij2, x_sous, y_sous;
    double inv_r2,inv_r6,terme_commun;
    double abs_distx, abs_disty;

    for ( int i=0;i<nbx_particules;i++)
    {
        f_x=0.0;
        f_y=0.0;

        for ( int j=0;j<nbx_particules;j++)
        {
            if (j != i)
            {
                
                x_sous = tab_par[j].x - tab_par[i].x;
                y_sous = tab_par[j].y - tab_par[i].y;

                abs_distx = fabs(x_sous);
                abs_disty = fabs(y_sous); 

                if (abs_distx > Lmoitie) {
                    x_sous = -(x_sous/abs_distx)*(L-abs_distx);
                }

                if (abs_disty > Lmoitie) {
                    y_sous = -(y_sous/abs_disty)*(L-abs_disty);
                }

                rij2 = x_sous * x_sous + y_sous * y_sous;

                if (rij2<Rc2){

                    inv_r2 = 1.0 / rij2;
                    inv_r6 = inv_r2 * inv_r2 * inv_r2;
                    terme_commun = 24.0 * (inv_r6*inv_r2 - 2*inv_r6*inv_r6*inv_r2);

                    f_x += terme_commun*x_sous;
                    f_y += terme_commun*y_sous;
                }
            }
        }

        tab_par[i].fx=f_x;
        tab_par[i].fy=f_y;
    }

}

//-----------------------------{resolution}-----------------------------

void resolution(Particule tab_par[],double delta_t)
{
    //je commence par déterminer f_indice_x :
    Particule tab_par_temp[nbx_particules];

    calcul_forces(tab_par);
    for ( int i=0;i<nbx_particules;i++)
    {
        
        tab_par[i].x = tab_par[i].x + tab_par[i].vx*delta_t + 0.5*tab_par[i].fx*(delta_t*delta_t);
        tab_par[i].y = tab_par[i].y + tab_par[i].vy*delta_t + 0.5*tab_par[i].fy*(delta_t*delta_t);

        if(tab_par[i].x > Lmoitie){ tab_par[i].x-=L;}
        else if(tab_par[i].x < -Lmoitie){ tab_par[i].x+=L;}

        if(tab_par[i].y > Lmoitie) {tab_par[i].y-=L;}
        else if(tab_par[i].y < -Lmoitie) {tab_par[i].y+=L;}

        tab_par_temp[i].x=tab_par[i].x;
        tab_par_temp[i].y=tab_par[i].y;
        tab_par_temp[i].fx=tab_par[i].fx;
        tab_par_temp[i].fy=tab_par[i].fy;

    }
    calcul_forces(tab_par_temp);
    for ( int i=0;i<nbx_particules;i++)
    {
        tab_par[i].vx = tab_par[i].vx + 0.5*(tab_par_temp[i].fx+ tab_par[i].fx)*delta_t;
        tab_par[i].vy = tab_par[i].vy + 0.5*(tab_par_temp[i].fy + tab_par[i].fy)*delta_t;
    }
}

//-----------------------------{main}-----------------------------

int main() {
    double max_energie_totale = -1.0; // Variable pour stocker le maximum de l'énergie totale
    double temps_ref = 0.0;

    double delta_t = 5.0e-4;

    Particule tab_particules[nbx_particules];

    initPoint(&tab_particules[0], 0.369163, 7.001539, 2.461984, 6.004601);

    initPoint(&tab_particules[1], 11.804992, 12.508632, 6.039470, 6.527302);

    initPoint(&tab_particules[2], 0.720203, 14.251818, 1.136606, 7.847989);

    initPoint(&tab_particules[3], 2.790903, 8.519731, 2.597047, 7.523878);

    initPoint(&tab_particules[4], 5.316689, 9.872263, 0.398034, 5.516499);

    initPoint(&tab_particules[5], 5.712756, 5.894922, 5.035953, 7.358547);

    initPoint(&tab_particules[6],8.860042, 7.814479, 6.068676, 9.158643);

    initPoint(&tab_particules[7], 12.048815, 1.838540, 1.266312, 8.278652);

    initPoint(&tab_particules[8], 8.840079, 5.592444, 4.283253, 3.763380);

    initPoint(&tab_particules[9],3.101077, 0.484085, 0.290682, 2.547519);

    char command[500];

    double t=0.0;

    FILE *out;

    int i,j;
    out = fopen ("testfin.txt", "w");

    if (!out) {
        perror("Erreur lors de l'ouverture du fichier");
        return EXIT_FAILURE;
    }

    for (j=0;j<3*1/delta_t;j++)
    {
        t+=delta_t;
        resolution(tab_particules,delta_t);
        double energie_totale = Energie_totale(tab_particules);

        //printf("%.10f\n",Energie_totale(tab_particules)); //verification energie totale

        fprintf (out, "%22.8g : %22.8g \n",t,\
        energie_totale);//Potentiel_LJ(&tab_particules[0],&tab_particules[1]));

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

    const int dim=54;
    double test[dim];
    double diviseur[]={ 1, 2, 3, 4, 5, 6, 7, 8, 9,
                    10, 20, 30, 40, 50, 60, 70, 80, 90,
                    100, 200, 300, 400, 500, 600, 700, 800, 900,
                    1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                    10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
                    100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000
                    };

    for (int i=0; i<dim; ++i) {
        test[i] = temps_ref/diviseur[i];
        printf(" ,%.10lf ",test[i]);
        if(i%5==0 ) printf("\n");
    }

     char nom_fichier[25];
    for (int k=0; k<dim; ++k) {

        sprintf(nom_fichier, "Energies/Etot_%d.txt",k);
        delta_t = test[k];

        Particule tab_particules_bis [nbx_particules];

        initPoint(&tab_particules_bis[0], 0.369163, 7.001539, 2.461984, 6.004601);
        initPoint(&tab_particules_bis[1], 11.804992, 12.508632, 6.039470, 6.527302);
        initPoint(&tab_particules_bis[2], 0.720203, 14.251818, 1.136606, 7.847989);
        initPoint(&tab_particules_bis[3], 2.790903, 8.519731, 2.597047, 7.523878);
        initPoint(&tab_particules_bis[4], 5.316689, 9.872263, 0.398034, 5.516499);
        initPoint(&tab_particules_bis[5], 5.712756, 5.894922, 5.035953, 7.358547);
        initPoint(&tab_particules_bis[6], 8.860042, 7.814479, 6.068676, 9.158643);
        initPoint(&tab_particules_bis[7], 12.048815, 1.838540, 1.266312, 8.278652);
        initPoint(&tab_particules_bis[8], 8.840079, 5.592444, 4.283253, 3.763380);
        initPoint(&tab_particules_bis[9], 3.101077, 0.484085, 0.290682, 2.547519);

        char command;
        double t = 0.0;
        FILE *out ; 
        out= fopen(nom_fichier, "w");

        for (int j = 0; j < 3 * 1 / delta_t; j++) 
        {
            t += delta_t;
            resolution(tab_particules_bis,delta_t);
            fprintf(out, "%22.8g : %22.8g \n", t, Energie_totale(tab_particules_bis));
            if (t >= temps_ref) break;
            
        }
        fclose(out);
    } 
    return 0;

}
