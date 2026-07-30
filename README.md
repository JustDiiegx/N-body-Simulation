# N-Body Gravitational Simulation (Sequential, OpenMP & CUDA)

An optimized N-Body physics simulator developed in **C** and **CUDA** that models gravitational interactions among $N$ celestial bodies in 3D space. The project explores algorithm optimization and parallel computing paradigms, scaling from a sequential CPU baseline to multi-core shared memory execution with **OpenMP** and massive parallelization on GPUs with **CUDA**.

Developed for the _Multicore Architecture Programming_ course at the University of Murcia.

## Features & Physical Model

- **Newtonian Gravitational Physics:** Solves 3D particle interactions using Newton's Law of Universal Gravitation and Newton's Second/Third Laws of Motion.
    
- **Singularity Softening:** Implements a softening parameter ($\epsilon$) to prevent division-by-zero errors and unrealistic velocity spikes during close-body collisions.
    
- **Numerical Integration:** Uses Euler's method for trajectory integration over discrete time steps ($\Delta t$).
    
- **Physical Correctness Validation:** Includes linear momentum conservation checks ($P = \sum m \cdot v$) to verify simulation accuracy across all versions.
    
- **3D Visualizer (Optional):** Interactive real-time 3D simulation renderer built with **Raylib**.
    

## Implementations & Architecture

### 1. Sequential Baseline ($O(n^2)$ Optimized)

- Leverages **Newton's Third Law** ($F_{ij} = -F_{ji}$) to compute forces using an upper triangular matrix, cutting the required force calculations and memory usage in half.
    

### 2. Multi-Core CPU Parallelization (OpenMP)

- Parallelizes the force calculation and velocity update loops across CPU threads using `#pragma omp parallel for schedule(dynamic, 1)` to handle dynamic loop workloads efficiently.
    
- **Speedup:** **2.53x** over sequential.
    

### 3. GPU Acceleration (NVIDIA CUDA)

- Redesigned to eliminate shared memory bottlenecks and triangular matrix lookup overhead.
    
- **Kernel Design:**
    
    - `kernel_fuerzas_y_velocidades`: Each GPU thread handles one body $i$, accumulating gravitational forces from all other bodies $j$ in local registers and updating velocity.
        
    - `kernel_posiciones`: Updates 3D spatial coordinates based on newly calculated velocities.
        
- **PCIe Bus Optimization:** Data transfers (`cudaMemcpy`) are isolated to initial setup and final teardown, keeping data residing entirely on device memory during loop iterations.
    
- **Speedup:** **4.61x** over sequential.
    

## Performance Benchmark

Tested with $N = 2000$ bodies over $500$ simulation steps:

|**Implementation**|**Execution Time**|**Speedup**|**Parallel Paradigm**|
|---|---|---|---|
|**Sequential Baseline**|54.18 s|1.00x|Single Thread (CPU)|
|**OpenMP**|21.39 s|**2.53x**|Multi-Core Shared Memory (CPU)|
|**CUDA**|11.75 s|**4.61x**|Massive Parallelism (GPU)|

## Prerequisites & Compilation

### Requirements

- `gcc` / `g++`
    
- `make`
    
- OpenMP runtime library
    
- NVIDIA CUDA Toolkit (`nvcc`)
    
- Raylib _(Optional, required only for the 3D visual version)_
    

### Build and Run

#### 1. OpenMP Version

```
gcc -O3 -fopenmp nbodyParallel.c -o nbodyParallel -lm
./nbodyParallel <number_of_bodies> <iterations>
# Example: ./nbodyParallel 2000 500
```

#### 2. CUDA Version

```
nvcc -O3 nbody_cuda.cu -o nbodyCuda
./nbodyCuda <number_of_bodies> <iterations>
# Example: ./nbodyCuda 2000 500
```
