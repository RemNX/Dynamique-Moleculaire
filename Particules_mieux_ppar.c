#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#define delta_t  10e-6
double delta_t;
#define nbx_particules   5

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

void resolution(Particule tab_par[])
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

int main() 
{
    char nom_fichier[20];
    for (int h=-9; h<=2;h++){
        for (int k=1;k<=9;k++){
        sprintf(nom_fichier,"amp_plot/Etot_%de%d.txt",k,h);
        delta_t = k*pow(10,h);

        //srand(time(NULL)); // pour les nombres aleatoires
        
        Particule tab_particules [nbx_particules];

        initPoint(&tab_particules[0], 2, 2.3,  4,  3);

        initPoint(&tab_particules[1], 1.5,  3,  3,  2);

        initPoint(&tab_particules[2],  1,  1.5,  4,  5);

        initPoint(&tab_particules[3],  1.7,  1.6,  2,  3);

        initPoint(&tab_particules[4],  2.4,  3,  3,  5);
        
        char command[500];

        double t=0.0;
        
        FILE *out;

        int i,j;
    /*  for ( i=0;i<10;i++)
        {
            tab_particules[i]=genererParticule();
            //printf("%f\n",tab_particules[0].vx);
        }
    */
        //out = fopen ("Etot.txt", "w");
        out = fopen (nom_fichier, "w");
        for (j=0;j<5000;j++)
        {
            t+=delta_t;
            resolution(tab_particules);
            //printf("%.10f\n",Energie_totale(tab_particules)); //verification energie totale
            fprintf (out, "%22.8g  %22.8g \n",t,\
            Energie_totale(tab_particules));//Potentiel_LJ(&tab_particules[0],&tab_particules[1]));
        }
        fclose (out);
        //sprintf (command, "gnuplot -p -e \"plot 'potentiel.txt' w l , \" ");
    /* sprintf(command, 
            "gnuplot -p -e \"set title 'Etot_{12}' font ',16'; \
            set xlabel 'r/d'; \
            set ylabel 'E/E0'; \
            set grid;\
            plot Etot_5e-6' w l\"");

        system (command); */

        }
    }
    return 0;
}