#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>

// En caso de que queramos utilizar la fórmula completa
const double G = 6.674*1e-10;

typedef struct{
    double x, y, z;
    double vx, vy, vz;
    double m;
} Cuerpo;

void gravedad(int nCuerpos, double* Fx, double* Fy, double* Fz, Cuerpo* cuerpos, int ldn, double t, double* t_fuerzas, double* t_velocidades, double* t_posiciones){
    int i, j;
    double t_start, t_end;

    memset(Fx, 0, sizeof(double) * nCuerpos * nCuerpos);
    memset(Fy, 0, sizeof(double) * nCuerpos * nCuerpos);
    memset(Fz, 0, sizeof(double) * nCuerpos * nCuerpos);

    t_start = omp_get_wtime();
    for(i = 0; i < nCuerpos; i++){

	for(j = i+1; j < nCuerpos; j++){
	    double dx = cuerpos[j].x - cuerpos[i].x;
	    double dy = cuerpos[j].y - cuerpos[i].y;
	    double dz = cuerpos[j].z - cuerpos[i].z;

	    double softening = 2.0f;

	    double distancia = dx*dx + dy*dy + dz*dz + softening;
            
	    double F = G * (cuerpos[i].m * cuerpos[j].m) / pow(distancia, 1.5);

	    Fx[i*ldn+j] += dx * F;
            Fy[i*ldn+j] += dy * F;
            Fz[i*ldn+j] += dz * F;
	}
    }
    t_end = omp_get_wtime();
    *t_fuerzas += (t_end - t_start);

    t_start = omp_get_wtime();
    for(i = 0; i < nCuerpos; i++){
        // Todas las fuerzas de este cuerpo
        double FCx = 0, FCy = 0, FCz = 0;
        int k;
        for(k = 0; k < nCuerpos; k++){
            // Accedemos a los valores calculados en
            // otras iteraciones anteriores
            if(k < i){
	        // Cambiamos el signo, si k tiene fuerza hacia i
	        // i tiene fuerza opuesta hacia k
            	FCx -= Fx[k*ldn+i];
              	FCy -= Fy[k*ldn+i];
              	FCz -= Fz[k*ldn+i];
            }
            // Accedemos a los valores calculados en
            // esta iteración
            else if(k > i){
            	FCx += Fx[i*ldn+k];
            	FCy += Fy[i*ldn+k];
            	FCz += Fz[i*ldn+k];
            }
        }
        // Añadimos la velocidad al cuerpo
        cuerpos[i].vx += t * (FCx / cuerpos[i].m);
        cuerpos[i].vy += t * (FCy / cuerpos[i].m);
        cuerpos[i].vz += t * (FCz / cuerpos[i].m);
    }
    t_end = omp_get_wtime();
    *t_velocidades += (t_end - t_start);
    
    // Una vez hemos calculado la aceleración
    // de todos los cuerpos en un instante concreto
    // podemos moverlos todos a la vez
    t_start = omp_get_wtime();
    for(i = 0; i < nCuerpos; i++){
    	cuerpos[i].x += t*cuerpos[i].vx;
	cuerpos[i].y += t*cuerpos[i].vy;
	cuerpos[i].z += t*cuerpos[i].vz;
    }
    t_end = omp_get_wtime();
    *t_posiciones += (t_end - t_start);
}

void momentoLineal(int nCuerpos, Cuerpo* cuerpos) {
    double Px = 0.0, Py = 0.0, Pz = 0.0;
    
    for(int i = 0; i < nCuerpos; i++) {
        Px += (double)cuerpos[i].m * (double)cuerpos[i].vx;
        Py += (double)cuerpos[i].m * (double)cuerpos[i].vy;
        Pz += (double)cuerpos[i].m * (double)cuerpos[i].vz;
    }

    printf("Momento Total del Sistema -> X: %e, Y: %e, Z: %e\n", Px, Py, Pz);
}


int main(int argc,char *argv[]){

    if(argc < 3){
	printf("Error de uso: indica una cantidad de cuerpos.\n");
    	return -1;
    }

    int nCuerpos = atoi(argv[1]);
    int iteraciones = atoi(argv[2]);
    
    double t = 10.0f;

    Cuerpo* cuerpos = (Cuerpo*)malloc(nCuerpos * sizeof(Cuerpo));

    double *Fx, *Fy, *Fz;

    Fx = (double *) malloc(sizeof(double)*nCuerpos*nCuerpos);
    Fy = (double *) malloc(sizeof(double)*nCuerpos*nCuerpos);
    Fz = (double *) malloc(sizeof(double)*nCuerpos*nCuerpos);


    srand(time(NULL));
    int i;
    for(i = 0; i < nCuerpos; i++){
    	cuerpos[i].x = (double)rand() / (double)RAND_MAX * 20.0f - 10.0f ;
	cuerpos[i].y = (double)rand() / (double)RAND_MAX * 20.0f - 10.0f;
	cuerpos[i].z = (double)rand() / (double)RAND_MAX * 20.0f - 10.0f;
	cuerpos[i].vx = 0.0f; //(double)rand() / (double)RAND_MAX;
	cuerpos[i].vy = 0.0f; //(double)rand() / (double)RAND_MAX;
	cuerpos[i].vz = 0.0f; //(double)rand() / (double)RAND_MAX;
    	cuerpos[i].m = (double)rand() / (double)RAND_MAX * 500000.0f + 1.0f;
    }
    
    double t_fuerzas = 0.0;
    double t_velocidades = 0.0;
    double t_posiciones = 0.0;
    
    double t_total_start = omp_get_wtime();
    
    int step;
    gravedad(nCuerpos, Fx, Fy, Fz, cuerpos, nCuerpos, t, &t_fuerzas, &t_velocidades, &t_posiciones);
    momentoLineal(nCuerpos, cuerpos);
    for(step = 1; step < iteraciones; step++){
      gravedad(nCuerpos, Fx, Fy, Fz, cuerpos, nCuerpos, t, &t_fuerzas, &t_velocidades, &t_posiciones);
    }
    momentoLineal(nCuerpos, cuerpos);
    
    double t_total_end = omp_get_wtime();
    double t_total = t_total_end - t_total_start;

    printf("\n================ REPORTE DE RENDIMIENTO ================\n");
    printf("Tiempo Total de Bucle:   %f segundos\n", t_total);
    printf("--------------------------------------------------------\n");
    printf("Fase 1 (Calc Fuerzas):   %f seg (%.2f%%)\n", t_fuerzas, (t_fuerzas/t_total)*100);
    printf("Fase 2 (Sumar Vels):     %f seg (%.2f%%)\n", t_velocidades, (t_velocidades/t_total)*100);
    printf("Fase 3 (Mover Pos):      %f seg (%.2f%%)\n", t_posiciones, (t_posiciones/t_total)*100);
    printf("========================================================\n\n");

    free(cuerpos);
    free(Fx);
    free(Fy);
    free(Fz);

    return 0;
}























