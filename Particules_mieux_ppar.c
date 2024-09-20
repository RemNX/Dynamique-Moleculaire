#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define delta_t  10e-4
#define nbx_particules   2

/*Definition Particule avec x,y,vitesse selon x
, vitesse selon y*/
typedef struct 
{
    double x;
    double y;
    double vx;
    double vy;
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
    return sqrt(pow(p1->x - p2->x, 2) + pow(p1->y - p2->y, 2));
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
double force_selon_x(Particule *p1, Particule *p2)
{
    double rij = distance(p1, p2);
    double terme_commun = (24/pow(rij,2))*pow((1/rij),6)*(2*pow((1/rij),6)-1);
    return terme_commun*(p2->x-p1->x);
}

/*Fonction qui calcule la force selon y en Lennard Jones entre 
deux particules*/
double force_selon_y(Particule *p1, Particule *p2)
{ 
    double rij = distance(p1, p2);
    double terme_commun = (24/pow(rij,2))*pow((1/rij),6)*(2*pow((1/rij),6)-1);
    return terme_commun*(p2->y-p1->y);
}

double Potentiel_LJ(Particule *p1, Particule *p2)
{
    double rij=distance(p1,p2);
    double U=4*(pow((1/rij),12)-pow((1/rij),6));
    return U;
}

double Energie_totale(Particule tab_par[])
{
    double E=0,U=0;
    for (int i=0;i<nbx_particules-1;i++)
    { 
        U+=Potentiel_LJ(&tab_par[i],&tab_par[i+1]);
        E+=(pow(tab_par[i].vx,2)+pow(tab_par[i].vy,2))/2 ;
    }
    E+=(pow(tab_par[nbx_particules-1].vx,2)+pow(tab_par[nbx_particules-1].vy,2))/2 ;
    return U+E;
}

/*Fonction qui resoud les equations differentielles d'une particule
designé par un "indice", et qui prend un tableau dans lequel sont 
stockées les particules, elle ne retourne rien car tout est modifié 
directement */
void resolution(Particule tab_par[], int indice)
{
    //je commence par déterminer f_indice_x :
    double f_x =0;
    double f_y =0;

    for ( int i=0;i<nbx_particules;i++)
    {
        if (indice != i)
        { 
            f_x += force_selon_x(&tab_par[indice],&tab_par[i]);
            f_y += force_selon_y(&tab_par[indice],&tab_par[i]);
        }
    }

    tab_par[indice].vx = tab_par[indice].vx + (f_x)*delta_t;
    tab_par[indice].vy  = tab_par[indice].vy + (f_y)*delta_t;
    tab_par[indice].x =  tab_par[indice].x + tab_par[indice].vx*delta_t;
    tab_par[indice].y = tab_par[indice].y + tab_par[indice].vy*delta_t; 
}  

int main() 
{
    //srand(time(NULL)); // pour les nombres aleatoires
    
    Particule tab_particules [nbx_particules];
    initPoint(&tab_particules[0],1.4,2,6,9);
    initPoint(&tab_particules[1],1.2,1.1,2,1);
    char command[500];
    FILE *out;

     int i,j;
   /*  for ( i=0;i<10;i++)
    {
        tab_particules[i]=genererParticule();
        //printf("%f\n",tab_particules[0].vx);
    }
 */
    out = fopen ("Ep_Lennard_Jones.txt", "w");
    for (j=0;j<20000;j++)
    {
        for ( i=0;i<nbx_particules;i++)
        {   
            resolution(tab_particules, i);
        }
        printf("%.10f\n",Energie_totale(tab_particules)); //verification energie totale
        fprintf (out, "%22.8g  %22.8g \n", \
        distance(&tab_particules[0],&tab_particules[1]),\
        Potentiel_LJ(&tab_particules[0],&tab_particules[1]));
    }
    fclose (out);
    //sprintf (command, "gnuplot -p -e \"plot 'potentiel.txt' w l , \" ");
   sprintf(command, 
        "gnuplot -p -e \"set title 'Simulation Potentiel Lennard Jones' font ',16'; \
        set xlabel 'r/d'; \
        set ylabel 'E/E0'; \
        set xrange [0.7:3]; \
        set yrange [-1.5:3.5]; \
        set grid;\
        plot 'Ep_Lennard_Jones.txt' w l\"");

    //system (command);

    return 0;
}