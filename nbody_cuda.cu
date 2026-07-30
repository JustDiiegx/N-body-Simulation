#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cuda_runtime.h>

const double G = 6.674*1e-10;

typedef struct{
    double x, y, z;
    double vx, vy, vz;
    double m;
} Cuerpo;

void momentoLineal(int nCuerpos, Cuerpo* cuerpos) {
    double Px = 0.0, Py = 0.0, Pz = 0.0;
    
    int i;
    for(i = 0; i < nCuerpos; i++) {
        Px += cuerpos[i].m * cuerpos[i].vx;
        Py += cuerpos[i].m * cuerpos[i].vy;
        Pz += cuerpos[i].m * cuerpos[i].vz;
    }

    printf("Momento Total del Sistema -> X: %e, Y: %e, Z: %e\n", Px, Py, Pz);
}

__global__ void kernel_fuerzas_y_velocidades(Cuerpo* cuerpos, int nCuerpos, double t) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j;

    if (i < nCuerpos) {
        double FCx = 0.0, FCy = 0.0, FCz = 0.0;
        
        for (j = 0; j < nCuerpos; j++) {
            if (i != j) {
                double dx = cuerpos[j].x - cuerpos[i].x;
                double dy = cuerpos[j].y - cuerpos[i].y;
                double dz = cuerpos[j].z - cuerpos[i].z;
                
                double softening = 2.0;
                
                double distancia = dx*dx + dy*dy + dz*dz + softening;
                
                double F = G * (cuerpos[i].m * cuerpos[j].m) / pow(distancia, 1.5);
                
                FCx += dx * F;
                FCy += dy * F;
                FCz += dz * F;
            }
        }
        
        cuerpos[i].vx += t * (FCx / cuerpos[i].m);
        cuerpos[i].vy += t * (FCy / cuerpos[i].m);
        cuerpos[i].vz += t * (FCz / cuerpos[i].m);
    }
}

__global__ void kernel_posiciones(Cuerpo* cuerpos, int nCuerpos, double t) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nCuerpos) {
        cuerpos[i].x += t * cuerpos[i].vx;
        cuerpos[i].y += t * cuerpos[i].vy;
        cuerpos[i].z += t * cuerpos[i].vz;
    }
}

int main(int argc, char *argv[]) {

    if(argc < 3){
        printf("Error de uso: ./nbody_cuda <nCuerpos> <iteraciones>\n");
        return -1;
    }

    int nCuerpos = atoi(argv[1]);
    int iteraciones = atoi(argv[2]);
    
    double t = 10.0;
    cudaError_t err = cudaSuccess;

    size_t size = nCuerpos * sizeof(Cuerpo);
    Cuerpo* h_cuerpos = (Cuerpo*)malloc(size);

    srand(time(NULL));
    int i;
    for(i = 0; i < nCuerpos; i++){
        h_cuerpos[i].x = (double)rand() / RAND_MAX * 20.0 - 10.0;
        h_cuerpos[i].y = (double)rand() / RAND_MAX * 20.0 - 10.0;
        h_cuerpos[i].z = (double)rand() / RAND_MAX * 20.0 - 10.0;
        h_cuerpos[i].vx = 0.0;
        h_cuerpos[i].vy = 0.0;
        h_cuerpos[i].vz = 0.0; 
        h_cuerpos[i].m = (double)rand() / RAND_MAX * 500000.0 + 1.0;
    }
    
    // Reservar memoria y copiar a GPU
    Cuerpo* d_cuerpos;
    err = cudaMalloc((void**)&d_cuerpos, size);
    
     if (err != cudaSuccess)
    {
        fprintf(stderr, "Fallo al reservar memoria en la GPU: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    err = cudaMemcpy(d_cuerpos, h_cuerpos, size, cudaMemcpyHostToDevice);
    
     if (err != cudaSuccess)
    {
        fprintf(stderr, "Fallo al copiar memoria de la CPU a la GPU: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    int threadsPerBlock = 256;
    int blocksPerGrid = (nCuerpos + threadsPerBlock - 1) / threadsPerBlock;

    // Eventos para medir el tiempo
    cudaEvent_t start_total, stop_total;
    cudaEvent_t start_fase1, stop_fase1;
    cudaEvent_t start_fase3, stop_fase3;
    
    cudaEventCreate(&start_total);
    cudaEventCreate(&stop_total);
    
    cudaEventCreate(&start_fase1);
    cudaEventCreate(&stop_fase1);
    
    cudaEventCreate(&start_fase3);
    cudaEventCreate(&stop_fase3);

    double t_fuerzas_vel = 0.0;
    double t_posiciones = 0.0;
    float ms = 0.0f;

    // Comenzamos la ejecución, hacemos la primera iteración por separado
    // para calcular el momento lineal al comienzo.
    cudaEventRecord(start_total);
    
    // Fases 1 y 2 (ahora juntas)
    cudaEventRecord(start_fase1);
    kernel_fuerzas_y_velocidades<<<blocksPerGrid, threadsPerBlock>>>(d_cuerpos, nCuerpos, t);
    cudaEventRecord(stop_fase1);
    cudaEventSynchronize(stop_fase1);
    
    cudaEventElapsedTime(&ms, start_fase1, stop_fase1);
    t_fuerzas_vel += (ms / 1000.0);
    
    // Fase 3
    cudaEventRecord(start_fase3);
    kernel_posiciones<<<blocksPerGrid, threadsPerBlock>>>(d_cuerpos, nCuerpos, t);
    cudaEventRecord(stop_fase3);
    cudaEventSynchronize(stop_fase3);
        
     cudaEventElapsedTime(&ms, start_fase3, stop_fase3);
     t_posiciones += (ms / 1000.0);
    
    // Copiamos de vuelta a la CPU para calcular el momento lineal inicial
    err = cudaMemcpy(h_cuerpos, d_cuerpos, size, cudaMemcpyDeviceToHost);
    
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Fallo al copiar memoria de la GPU a la CPU: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    momentoLineal(nCuerpos, h_cuerpos);


    for(int step = 1; step < iteraciones; step++) {
        // Fases 1 y 2
        cudaEventRecord(start_fase1);
        kernel_fuerzas_y_velocidades<<<blocksPerGrid, threadsPerBlock>>>(d_cuerpos, nCuerpos, t);
        cudaEventRecord(stop_fase1);
        cudaEventSynchronize(stop_fase1); // Sincronizamos para medir exactamente este kernel
        
        cudaEventElapsedTime(&ms, start_fase1, stop_fase1);
        t_fuerzas_vel += (ms / 1000.0);

        // Fase 3
        cudaEventRecord(start_fase3);
        kernel_posiciones<<<blocksPerGrid, threadsPerBlock>>>(d_cuerpos, nCuerpos, t);
        cudaEventRecord(stop_fase3);
        cudaEventSynchronize(stop_fase3);
        
        cudaEventElapsedTime(&ms, start_fase3, stop_fase3);
        t_posiciones += (ms / 1000.0);
    }

    cudaEventRecord(stop_total);
    cudaEventSynchronize(stop_total);
    
    float total_ms = 0.0f;
    cudaEventElapsedTime(&total_ms, start_total, stop_total);
    double t_total = total_ms / 1000.0;

    // Volvemos a calcular el momento lineal tras la última iteración
    err = cudaMemcpy(h_cuerpos, d_cuerpos, size, cudaMemcpyDeviceToHost);
    
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Fallo al copiar memoria de la GPU a la CPU: %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    momentoLineal(nCuerpos, h_cuerpos);

    printf("\n================ REPORTE DE RENDIMIENTO (CUDA) ================\n");
    printf("Tiempo Total de Bucle:       %f segundos\n", t_total);
    printf("---------------------------------------------------------------\n");
    printf("Fases 1 y 2 (Fuerzas+Vels):  %f seg (%.2f%%)\n", t_fuerzas_vel, (t_fuerzas_vel/t_total)*100);
    printf("Fase 3 (Mover Posiciones):   %f seg (%.2f%%)\n", t_posiciones, (t_posiciones/t_total)*100);
    printf("===============================================================\n\n");

    cudaEventDestroy(start_total);
    cudaEventDestroy(stop_total);
    cudaEventDestroy(start_fase1);
    cudaEventDestroy(stop_fase1);
    cudaEventDestroy(start_fase3);
    cudaEventDestroy(stop_fase3);
    cudaFree(d_cuerpos);
    free(h_cuerpos);

    return 0;
}
