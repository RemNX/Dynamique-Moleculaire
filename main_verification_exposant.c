#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//-----------------------------{Définition des constantes}-----------------------------

#define nbx_particules   10
const int sigma=1;
#define Rc   2.5*sigma
#define pi   M_PI
#define MIN_DISTANCE 2.0
#define L 50//Taille de la boîte 
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

void resolution_verlet(Particule tab_par[],double delta_t)
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
double mc_ax_plus_b (double x[], double y[], double *a, double *b, double n)
{
    int i;
    double si, sx, sy, sxy, sxx;
    double delta;
    si = sx = sxy = sxx = 0;
    for (i = 0; i <= n - 1; i++)
      {
        //if (x[i]>-7 && x[i]<-2  )// Euler
       if (x[i]>-8 && x[i]<-4)//Verlet
        {si += 1;
        sxx += pow (x[i], 2);
        sxy += x[i] * y[i];
        sx += x[i];
        sy += y[i];}
      }

    delta = (sxx * si) - (sx * sx);
    *a = (1 / delta) * (sxy * si - sy * sx);
    *b = (1 / delta) * (sxx * sy - sxy * sx);
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
        resolution_verlet(tab_particules,delta_t);
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
        printf(" ,%.10lf ",test[i]);
        if(i%5==0 ) printf("\n");
    }
    
    out= fopen("test.txt", "w");
    //double E0=293.67526; //Euler
    double E0=293.71623; //verlet
    for (int k=0; k<dim; ++k) {

        delta_t = test[k];
        
        initPoint(&tab_particules[0], 0.369163, 7.001539, 2.461984, 6.004601);
        initPoint(&tab_particules[1], 11.804992, 12.508632, 6.039470, 6.527302);
        initPoint(&tab_particules[2], 0.720203, 14.251818, 1.136606, 7.847989);
        initPoint(&tab_particules[3], 2.790903, 8.519731, 2.597047, 7.523878);
        initPoint(&tab_particules[4], 5.316689, 9.872263, 0.398034, 5.516499);
        initPoint(&tab_particules[5], 5.712756, 5.894922, 5.035953, 7.358547);
        initPoint(&tab_particules[6], 8.860042, 7.814479, 6.068676, 9.158643);
        initPoint(&tab_particules[7], 12.048815, 1.838540, 1.266312, 8.278652);
        initPoint(&tab_particules[8], 8.840079, 5.592444, 4.283253, 3.763380);
        initPoint(&tab_particules[9], 3.101077, 0.484085, 0.290682, 2.547519);

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
