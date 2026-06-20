# GeoFLAC Tutorial 04: Parallel Scaling & Thread Configuration

This tutorial demonstrates how to configure parallel execution in GeoFLAC using OpenMP, and how to measure the solver's parallel performance (scaling, speedup, and efficiency).

---

## 1. Introduction to OpenMP Parallelism in GeoFLAC

GeoFLAC is a parallelized geodynamic solver using **OpenMP (Open Multi-Processing)**. Critical loops (such as stress updates, coordinate changes, and thermal diffusion) are parallelized to run concurrently across multiple CPU cores.

The number of threads used by the solver is controlled by the standard environment variable **`OMP_NUM_THREADS`**. Configuring this parameter correctly is essential to optimize solver speed.

---

## 2. Model Setup (`plastic.inp`)

To get a measurable and meaningful execution time for parallel scaling, we reuse the simple elastoplastic bar compression test (`plastic.inp` from Tutorial 3) and scale the mesh resolution to:
* **Grid Resolution**: $50 \times 150$ elements ($7500$ elements total).
* **Grid Size**: $1000 \text{ m} \times 3000 \text{ m}$.
* **Rheology**: Elastoplastic Mohr-Coulomb bar (`irheol = 6`).
* **Boundary Conditions**: Horizontal shortening (compression in Z) from the top at $10^{-9} \text{ m/s}$, with the bottom boundary fixed in Z (free-slip).
* **Max Time**: $0.5 \text{ Kyr}$ (running for 450 steps).

With this grid resolution, a single-thread execution takes approximately $6.6\text{ seconds}$, providing a healthy baseline to evaluate parallel speedup.

---

## 3. How to Find the Number of CPU Cores

Before configuring threads, it is useful to know how many CPU cores (logical processors) are available on your system. Here is how to find out on different operating systems:

### Linux
In the terminal, run any of the following commands:
* **`nproc`**: Prints the number of processing units available (recommended).
* **`lscpu`**: Displays detailed CPU architecture information (look for the `CPU(s):` line).
* **`grep -c ^processor /proc/cpuinfo`**: Counts the logical cores listed in CPU info.

### macOS
In the terminal, run:
* **`sysctl -n hw.ncpu`**: Prints the number of logical CPU cores.

### Windows
* **PowerShell**: Run `(Get-CimInstance Win32_ComputerSystem).NumberOfLogicalProcessors`
* **Command Prompt**: Run `wmic cpu get NumberOfCores,NumberOfLogicalProcessors`
* **Task Manager**: Open Task Manager (`Ctrl+Shift+Esc`), select the **Performance** tab, choose **CPU**, and read the **Cores** and **Logical processors** counts at the bottom.

---

## 4. How to Set `OMP_NUM_THREADS`

You can configure the thread count in your terminal in two ways:

### Option 1: Temporary (Single Command)
To run the solver with a specific number of threads for a single run, prepend the environment variable directly to the command:
```bash
OMP_NUM_THREADS=4 ../../src/flac plastic.inp
```

### Option 2: Persistent (Shell Session)
To set the thread count for all subsequent runs in your current terminal session, use `export`:
```bash
export OMP_NUM_THREADS=4
../../src/flac plastic.inp
```

To verify the setting, you can print the environment variable:
```bash
echo $OMP_NUM_THREADS
```

---

## 5. Measuring Parallel Scaling

To measure the parallel scaling of the solver, we run the same simulation under different thread counts ($N = 1, 2, 4, 8, 12, 16$) and record the wall-clock execution time ($T_N$). From these, we calculate:
* **Speedup ($S_N$)**: How much faster the parallel run is compared to the single-thread baseline:
  $$S_N = \frac{T_1}{T_N}$$
* **Parallel Efficiency ($E_N$)**: How effectively the solver utilizes the allocated processor cores:
  $$E_N = \frac{S_N}{N} \times 100\%$$

### Automated Scaling Script (`measure_scaling.py`)
We have provided a helper Python script `measure_scaling.py` in this directory. It automatically cleans output files, runs the simulation across multiple thread counts, and generates a scaling report:
```bash
python3 measure_scaling.py
```

---

## 6. Parallel Benchmark Results

The following measurements were taken on a **16-core CPU**:

| Threads ($N$) | Time ($T_N$) | Speedup ($S_N$) | Parallel Efficiency ($E_N$) |
|:---:|:---:|:---:|:---:|
| 1 | 6.56 s | 1.00x | 100.0% |
| 2 | 3.43 s | 1.91x | 95.6% |
| 4 | 1.91 s | 3.44x | 86.1% |
| 8 | 1.29 s | 5.08x | 63.5% |
| 12 | 1.31 s | 5.01x | 41.8% |
| 16 | 1.48 s | 4.45x | 27.8% |

### Analysis of Scaling Behavior
* **Excellent Scaling ($N \le 4$)**: Near-linear speedup and high efficiency ($>85\%$). For small to moderate thread counts, the CPU cores are utilized very effectively.
* **Diminishing Returns ($N \ge 8$)**: Speedup saturates around $5.0\text{x}$ at 8 threads.
* **Performance Drop ($N = 16$)**: Running with 16 threads is slightly slower than 8 threads.
* **Reasons for Saturation**:
  1. **Memory Bandwidth Bottleneck**: Finite-difference/finite-element geodynamic solvers are highly memory-bound. As more threads request grid and marker data from the RAM simultaneously, memory bandwidth saturates.
  2. **Thread Management Overhead**: For smaller grid sizes ($50 \times 150$), the computational load per thread is small. The overhead of spawning, synchronizing, and joining OpenMP threads becomes comparable to the time saved by parallel computation.
  
> [!TIP]
> **Best Practice**: For maximum throughput on medium-sized models, setting `OMP_NUM_THREADS` to **4 or 8** often provides the best balance of speedup and CPU resource conservation.
