#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//#define delta_t  10e-6
#define nbx_particules   10

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
/*Calcul distance entre 2 particules*/
double distance(Particule *p1, Particule *p2) {
    return sqrt(pow(p2->x - p1->x, 2) + pow(p2->y - p1->y, 2));
}

/*Fonction qui genere des positions et vitesse
aleatoires, pour une particule*/
Particule genererParticule()
{
    Particule p;
    p.x = (double)rand()/ RAND_MAX;
    p.y = (double)rand() / RAND_MAX;
    p.vx = (double)rand() / RAND_MAX;
    p.vy = (double)rand() / RAND_MAX;
    return p;
}

/*Fonction qui calcule la force selon x en Lennard Jones entre
deux particules*/
double force_selon_x(Particule *p2, Particule *p1)
{
    double rij = distance(p1, p2);
    double terme_commun = 24*(pow(1/rij,8)-2*pow(1/rij,14));
    return terme_commun*(p2->x-p1->x);
}

/*Fonction qui calcule la force selon y en Lennard Jones entre
deux particules*/
double force_selon_y(Particule *p2, Particule *p1)
{
    double rij = distance(p1, p2);
    double terme_commun = 24*(pow(1/rij,8)-2*pow(1/rij,14));
    return terme_commun*(p2->y-p1->y);
}

double Potentiel_LJ(Particule *p2, Particule *p1)
{
    double rij=distance(p1,p2);
    double U=4*(pow((1/rij),12)-pow((1/rij),6));
    return U;
}

double Energie_totale(Particule tab_par[])
{
    double E=0.0,U=0.0;

    for (int i=0;i<nbx_particules;i++)
    {
        for (int j=i+1;j<nbx_particules;j++)
        {

            U+=Potentiel_LJ(&tab_par[j],&tab_par[i]);

        }
        E+=(pow(tab_par[i].vx,2)+pow(tab_par[i].vy,2))/2;

    }
    return U+E;
}

/*Fonction qui resoud les equations differentielles d'une particule
designé par un "indice", et qui prend un tableau dans lequel sont
stockées les particules, elle ne retourne rien car tout est modifié
directement */
void calcul_forces (Particule tab_par[])
{
    //je commence par déterminer f_indice_x :
    double f_x;
    double f_y;

    for ( int i=0;i<nbx_particules;i++)
    {
        f_x=0.0;
        f_y=0.0;
        for ( int j=0;j<nbx_particules;j++)
        {
            if (j != i)
            {
                f_x += force_selon_x(&tab_par[j],&tab_par[i]);
                f_y += force_selon_y(&tab_par[j],&tab_par[i]);
            }
        }

        tab_par[i].fx=f_x;
        tab_par[i].fy=f_y;
    }

}

void resolution(Particule tab_par[],double delta_t)
{
    //je commence par déterminer f_indice_x :
    calcul_forces(tab_par);
    for ( int i=0;i<nbx_particules;i++)
    {
        tab_par[i].vx = tab_par[i].vx + tab_par[i].fx*delta_t;
        tab_par[i].vy  = tab_par[i].vy + tab_par[i].fy*delta_t;
        tab_par[i].x =  tab_par[i].x + tab_par[i].vx*delta_t;
        tab_par[i].y = tab_par[i].y + tab_par[i].vy*delta_t;


       /*  printf("%d %15.4f %15.4f %15.4f\n",i,(pow(tab_par[i].vx,2)+pow(tab_par[i].vy,2)) \
        ,tab_par[i].fx,tab_par[i].fy); */
    }


}

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
    out = fopen ("testou.txt", "w");

    for (j=0;j<3*1/delta_t;j++)
    {
        t+=delta_t;
        resolution(tab_particules,delta_t);
        double energie_totale = Energie_totale(tab_particules);

        //printf("%.10f\n",Energie_totale(tab_particules)); //verification energie totale

        fprintf (out, "%22.8g : %22.8g \n",t,\
        Energie_totale(tab_particules));//Potentiel_LJ(&tab_particules[0],&tab_particules[1]));

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

    const int dim=30;
    double test[dim];

    for (int i=0; i<dim; ++i) {
        double diviseur = 100*(i+1);
        test[i] = temps_ref/diviseur;
        printf(" ,%lf ",test[i]);
    }

    char nom_fichier[20];
    for (int k=0; k<dim; ++k) {

        sprintf(nom_fichier, "Etot_%d.txt",k);
        delta_t = test[k];

        Particule tab_particules_bis[nbx_particules];

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
        FILE *out;

        out = fopen(nom_fichier, "w");
        for (int j = 0; j < 3 * 1 / delta_t; j++) {
            t += delta_t;
            resolution(tab_particules_bis,delta_t);
            fprintf(out, "%22.8g : %22.8g \n", t, Energie_totale(tab_particules_bis));
        }
        fclose(out);
    }
    return 0;

}
