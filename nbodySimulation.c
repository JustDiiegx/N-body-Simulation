#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <raylib.h>
#include <string.h>

// En caso de que queramos utilizar la fórmula completa
const double G = 6.674*1e-10;

typedef struct{
    double x, y, z;
    double vx, vy, vz;
    double m;
} Cuerpo;

void gravedad(int nCuerpos, double* Fx, double* Fy, double* Fz, Cuerpo* cuerpos, int ldn, double t){
    int i, j;
    
    memset(Fx, 0, sizeof(double) * nCuerpos * nCuerpos);
    memset(Fy, 0, sizeof(double) * nCuerpos * nCuerpos);
    memset(Fz, 0, sizeof(double) * nCuerpos * nCuerpos);

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
    
    // Una vez hemos calculado la aceleración
    // de todos los cuerpos en un instante concreto
    // podemos moverlos todos a la vez
    for(i = 0; i < nCuerpos; i++){
    	cuerpos[i].x += t*cuerpos[i].vx;
	cuerpos[i].y += t*cuerpos[i].vy;
	cuerpos[i].z += t*cuerpos[i].vz;
    }
}


int main(int argc,char *argv[]){

    if(argc < 2){
	printf("Error de uso: indica una cantidad de cuerpos.\n");
    	return -1;
    }

    int nCuerpos = atoi(argv[1]);
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

    // Apartado gráfico para la visualización

    const int screenWidth = 1000;
    const int screenHeight = 800;

    InitWindow(screenWidth, screenHeight, "Problema de los N-Cuerpos");

    // Define the camera to look into our 3d world
    Camera3D camera = { 0 };
    camera.position = (Vector3){ 0.0f, 10.0f, 30.0f };  // Camera position
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };      // Camera looking at point
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };          // Camera up vector
    camera.fovy = 45.0f;                                // Camera field-of-view Y
    camera.projection = CAMERA_PERSPECTIVE;             // Camera mode type

    SetTargetFPS(30);

    while(!WindowShouldClose()){
    
    	gravedad(nCuerpos, Fx, Fy, Fz, cuerpos, nCuerpos, t);

	BeginDrawing();

	    ClearBackground(BLACK);
	    BeginMode3D(camera);
	    DrawGrid(20, 8.0f);
            for(i = 0; i < nCuerpos; i++){
	    	Vector3 pos = { cuerpos[i].x, cuerpos[i].y, cuerpos[i].z };
		DrawSphere(pos, cuerpos[i].m * 0.000001, WHITE);
	    }
	    EndMode3D();
	    DrawFPS(10,10);
    	EndDrawing();
    }
    CloseWindow();
    
    free(cuerpos);
    free(Fx);
    free(Fy);
    free(Fz);

    return 0;
}























