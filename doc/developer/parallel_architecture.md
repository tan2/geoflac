# GeoFLAC Developer Guide: Parallelization & GPU Acceleration Architecture

This document describes the parallel execution model of GeoFLAC. It details the hybrid OpenMP (CPU) and OpenACC (GPU) implementation, memory synchronization strategies, and best practices for writing GPU-compatible solver routines.

---

## 1. Dual-Parallelization Model

To handle massive geodynamical simulations efficiently, GeoFLAC supports both multi-core CPU computing and graphics card (GPU) acceleration within the same code base.

```
                  +--------------------------------+
                  |      Main Loop (Host CPU)      |
                  +--------------------------------+
                                   |
                  +--------------------------------+
                  |  !$ACC DATA COPYIN / COPYOUT   |
                  +--------------------------------+
                                   |
                  +----------------+----------------+
                  |                                 |
         [ OpenMP Directives ]             [ OpenACC Directives ]
                  |                                 |
       Parallel CPU Threads (OMP)       Massive GPU Threads (ACC)
                  |                                 |
                  +----------------+----------------+
                                   |
                  +--------------------------------+
                  |       !$ACC END DATA           |
                  +--------------------------------+
```

*   **OpenMP (CPU)**: Parallelizes operations across multiple host CPU cores using shared memory.
*   **OpenACC (GPU)**: Offloads compute-intensive loops to NVIDIA GPUs using compiler directives (e.g., compiled with NVHPC compiler).
*   Both paradigms exist side-by-side. The compiler ignores OpenACC directives when compiling for CPU-only execution, and ignores OpenMP directives when building for GPU offloading.

---

## 2. Memory Management (OpenACC Data Clauses)

GPU devices have separate physical memory (VRAM) from the host system memory (RAM). Offloading execution requires copying data back and forth, which is the primary bottleneck in GPU computing.

### A. Structured Data Regions
To minimize data transfer overhead, GeoFLAC defines a unified **structured data region** wrapping the main simulation loop:
```fortran
!$acc data copy(cord, temp, phase, visc, stress, strain, aps, fmelt, fmagma, area)
do nloop = 1, nsteps
    ! Calculation loops occur on the GPU without transferring arrays back to CPU
enddo
!$acc end data
```
Inside this data region, variables remain resident in GPU VRAM across timesteps.

### B. Explicit Memory Synchronization
If a subroutine on the host CPU needs to access data calculated on the GPU (or vice-versa), developers must explicitly synchronize using `update` directives:
*   **Host-to-Device Sync**:
    ```fortran
    !$acc update device(cord, temp)
    ```
    Copies updated host arrays into GPU memory. Used during initialization, input parsing, or grid remeshing.
*   **Device-to-Host Sync**:
    ```fortran
    !$acc update host(phase, visc)
    ```
    Copies calculated results from GPU memory back to host CPU memory. Used right before writing output files.

---

## 3. GPU-Compatible Subroutines (`!$ACC routine`)

If a loop on the GPU calls an external subroutine or function, that subroutine/function must be compiled to run on the device. This is done by adding the `!$ACC routine` directive in the subroutine definition.

### A. Routine Parallelism Levels
You must specify the degree of parallelism inside the routine:
*   **`seq` (Sequential)**: The routine is run sequentially by a single GPU thread. Most helper math functions and single-marker update functions are sequential.
    ```fortran
    subroutine count_phase_ratio(j, i)
    !$ACC routine seq
    ! Logic executed by a single thread
    end subroutine count_phase_ratio
    ```
*   **`worker` / `vector`**: The routine is executed in parallel across a group of threads.
    ```fortran
    subroutine newphase2marker(n, iphase)
    !$ACC routine worker
    ! Parallelized work executed inside a worker team
    end subroutine newphase2marker
    ```

---

## 4. Thread Safety & Race Conditions

When thousands of GPU threads run simultaneously, concurrent writes to the same memory cell will cause race conditions.

### A. Atomic Updates
If multiple threads must add or write values to a shared array index (such as accumulating marker compositions into grid element phase ratios), developers **must** use atomic directives:
```fortran
!$acc parallel loop async(1)
do kk = 1, nmarkers
    ! get element indices (j, i) of the marker
    ...
    ! Atomic increment to prevent multiple threads from overwriting each other
    !$acc atomic update
    phase_count(j, i) = phase_count(j, i) + 1.0d0
enddo
```
*   **`!$acc atomic update`**: Ensures read-modify-write operations are performed atomically on the GPU.
*   **`!$acc atomic write`**: Ensures simple write operations are atomic.
*   **OpenMP Fallback**: Always write the corresponding OpenMP CPU atomic directive (`!$omp atomic update`) next to the ACC directive for CPU thread safety.
