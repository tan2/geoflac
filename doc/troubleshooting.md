# Troubleshooting GeoFLAC Runtime Aborts

This guide lists the common runtime abort codes (`STOP` statements) that you may encounter when executing simulations, their causes, and how to resolve them.

---

## 1. Grid Geometry & Zone Configuration Aborts

These errors occur during initialization inside the coordinate setup module (`src/init_cord.f90`) and grid area calculator (`src/init_areas.f90`).

### `STOP 16` — Number of zones is not odd (X-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Number of zones is not odd.(X-direction)`
* **Cause:** The number of grid zones in the X-direction (`nzonx` parameter in the `.inp` file) is an even number (e.g., `2` or `4`).
* **Why this is a constraint:** The mesh generation algorithm (`mesh1`) requires an odd number of zones. The odd-numbered zones define regions of constant spacing, while the even-numbered zones act as transitional regions between them. Having an even number of zones causes memory out-of-bounds access.
* **Resolution:** Set `nzonx` to an odd number (e.g., `1`, `3`, `5`) or `0` for a regular grid.

### `STOP 17` — Nelem in zones is not correct (X-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Nelem in zones is not correct.(X-direction)`
* **Cause:** The sum of elements across all X-direction zones does not match the total elements specified on the mesh parameters line (`nex`).
* **Resolution:** Adjust the number of elements in each zone so that their sum exactly equals the first integer on the Mesh Parameters line (`nex = nx - 1`).

### `STOP 18` — Number of zones is not odd (Z-direction / Y-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Number of zones is not odd.(Y-direction)`
* **Cause:** The number of grid zones in the Z-direction (`nzony` parameter in the `.inp` file) is an even number (e.g., `2` or `4`).
* **Resolution:** Change `nzony` to an odd number (e.g., `1`, `3`, `5`) or `0` for a regular grid.

### `STOP 19` — Nelem in zones is not correct (Z-direction / Y-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Nelem in zones is not correct.(Y-direction)`
* **Cause:** The sum of elements across all Z-direction zones does not match the total elements specified on the mesh parameters line (`nez`).
* **Resolution:** Adjust the number of elements in each zone so that their sum exactly equals the second integer on the Mesh Parameters line (`nez = nz - 1`).

### `STOP 41` — Negative Element Area during Initialization
* **Error Message in `sys.msg`:** `INIT_AREAS: Negative area!`
* **Cause:** The coordinates specified for the grid result in a negative or zero-area element at the very beginning of the simulation. This is typically caused by:
  1. A bad custom coordinate input file where nodes overlap, intersect, or are defined in an inverted order.
  2. Extremely warped zone size scaling that forces nodes to cross each other.
* **Resolution:** Check your initial grid parameters or the custom coordinate input file. Ensure that nodes are defined in a clean, non-intersecting grid and that coordinate boundaries do not cross.

---

## 2. File & Path Aborts

These errors are related to external file access, parsing errors, or model coordinate and state mismatches.

### `STOP 21` — Coordinates/Temperature File Error
* **Error Message in `sys.msg`:** `INIT_CORD: Cannot open file with initial coordinates!` or `INIT_CORD: Error reading file with initial coordinates!` or similar error in `INIT_TEMP` for temperature profiles.
* **Cause:** The parameter `ircoord` is set to `1` (which tells the solver to read coordinates from an external file), but the file specified by `coordfile` (e.g., `coord.dat`) is missing or cannot be opened/read. Alternatively, a temperature profile file could not be read.
* **Resolution:** 
  1. Ensure that the coordinates file is in your running directory.
  2. Or, set `ircoord` to `0` to have the solver automatically generate the mesh coordinates.

### `STOP 11` — Input Parameter File Parsing Error
* **Error Message in stdout:** `Error reading file "<filename>" at line <line_number>`
* **Cause:** The solver encountered a parsing error (e.g., incorrect format, wrong data type, missing parameter, or trailing garbage text) while parsing the `.inp` parameter file at the specified line.
* **Resolution:** Open your `.inp` file, go to the reported line number, and check the format of the parameter value. Ensure it matches the expected type (integer, float, or string) and that no required entries on that line are missing.

### `STOP` — Unexpected End of Input File (EOF)
* **Error Message in stdout:** `AdvanceToNextInputLine: EOF reached!`
* **Cause:** The solver reached the end of the `.inp` file unexpectedly while looking for the next non-comment parameter line. This means the input file is incomplete or missing necessary parameters (e.g., output parameters or process control configurations at the end of the file).
* **Resolution:** Compare your `.inp` file with a complete, working example (e.g., from the benchmarks/examples). Add any missing lines at the end of the file.

### `STOP` — Custom Coordinate Grid Spacing Conflict
* **Error Message in stdout:** `ircoord is set, but nzonx is not 0!` or `ircoord is set, but nzony is not 0!`
* **Cause:** The parameter `ircoord` is set to `1` (or greater), which instructs the solver to read node coordinates from an external file (`coordfile`), but the grid zone count `nzonx` or `nzony` is set to a non-zero value.
* **Why this is a constraint:** When custom coordinates are loaded, the built-in zone spacing generator (`mesh1`) is bypassed. Specifying a non-zero number of zones is contradictory and invalid.
* **Resolution:** Edit your `.inp` file. If `ircoord` is set to `1` (or greater), make sure the subsequent lines for `nzonx` and `nzony` are both set to `0` (e.g., `0   ! nzonx` and `0   ! nzony`).

### `STOP 999` — Model Restart Dimension Mismatch
* **Error Message in stdout:** `Wrong marker ntriag ...`
* **Cause:** A restart was triggered because a binary save state (`_contents.rs` and matching `.rs` files) was found in the execution directory, but the grid dimensions (`nx`, `nz`) or element sizes specified in the `.inp` file do not match the binary layout.
* **Resolution:** Make sure that when restarting a run, you do not modify the grid dimensions or mesh layout in the `.inp` file. If you wish to start a new model, delete the `.rs` files and `_contents.rs` first.

---

## 3. Boundary & Configuration Aborts

These errors occur due to invalid, incompatible, or array-limit-exceeding configuration parameter inputs.

### `STOP 26` — Inhomogeneities Array Limit Exceeded
* **Error Message in `sys.msg`:** `Read_params: Increase arrays for inhomogenities`
* **Cause:** The number of inhomogeneities (`inhom` parameter in the `.inp` file) exceeds the maximum compile-time limit `maxinh` defined in the solver source code.
* **Resolution:** 
  1. Reduce the number of inhomogeneities specified in the `.inp` file.
  2. Or, increase the value of `maxinh` inside `src/params.f90` (and `src/arrays.f90` if applicable) and recompile the solver.

### `STOP 55` — Mid-Plate Push Boundary Mismatch
* **Cause:** The parameter `nofside` is set to `5` (mid-plate velocity push segment), but the corresponding segment boundary type `nbc` is not set to `10` (horizontal velocity `velx`).
* **Resolution:** Set `nbc` to `10` on the mid-plate boundary condition row.

### `STOP` — Boundary Side Mismatch for Viscous Velocity
* **Error Message in stdout:** `Wrong side for velbv_visc`
* **Cause:** A viscous-dependent velocity boundary condition (`velbc_visc`) is defined on a boundary side (`nofside`) other than side `1` (left) or side `3` (right).
* **Resolution:** Ensure that the boundary condition for viscous velocity is only assigned to the left or right boundaries.

### `STOP` — Illegal Remeshing Mode
* **Error Message in `sys.msg`:** `Illegal remeshing mode! Allowable - 1, 3, 4or 11`
* **Cause:** The parameter `mode_rem` under the REMESHING section in your `.inp` file is set to an unsupported value.
* **Resolution:** Edit the `.inp` file and change `mode_rem` to one of the allowed integer options:
  * `1`: Standard grid remeshing.
  * `3`: Remeshing with erosion/deposition.
  * `4`: Remeshing with topography smoothing.
  * `11`: Adaptive remeshing based on element deformation.

### `STOP` — Remeshing Test Frequency Mismatch
* **Error Message in `sys.msg`:** `ntest_rem must be multiples of process frequency`
* **Cause:** The remeshing test frequency `ntest_rem` is not a multiple of one of the process frequency parameters: `ifreq_rmasses`, `ifreq_imasses`, `ifreq_visc`, or `ifreq_avgsr`.
* **Why this is a constraint:** The solver must synchronize the remeshing checks with the solver steps where masses, viscosity, and average strain rates are recalculated.
* **Resolution:** Adjust `ntest_rem` in the `.inp` file so that it is exactly divisible by (or a multiple of) the process control frequencies (`ifreq_rmasses`, `ifreq_imasses`, `ifreq_visc`, `ifreq_avgsr`). For example, if all frequencies are set to `10`, then `ntest_rem` must be a multiple of 10 (e.g., `10`, `20`, `50`).

---

## 4. Runtime Physics & Stability Aborts

These errors occur during model calculations and are caused by numerical instability, massive deformation, or marker problems.

### `STOP 22` — Maxwell Viscoelastic Stability Limit
* **Cause:** The numerical time step `dt` is too large relative to the viscoelastic relaxation time of the deforming elements ($dt / (2 \cdot \eta / \mu) > 0.5$, where $\eta$ is viscosity and $\mu$ is shear modulus).
* **Resolution:** 
  1. Decrease `dt_scale` (e.g., to `0.5` or smaller) under Process Control in your `.inp` file to force a smaller time step.
  2. Increase the minimum viscosity limit (`vis_min`) to prevent highly fluid elements from destabilizing the solver.

### `STOP 40` — Grid Inversion (Excessive Element Deformation)
* **Error Message in `sys.msg`:** `Negative area!` during calculation steps.
* **Cause:** Shear deformation has caused the element grid to fold over or invert coordinates (negative area during move step), causing code to abort to prevent garbage output.
* **Resolution:** 
  1. Enable and configure remeshing (`ny_rem = 1`) to reset the grid when deformation is high.
  2. Reduce your grid spacing or boundary velocity to slow down deformation rates.

### `STOP 133` / `STOP 134` — Marker Tracking Out-of-Bounds
* **Cause:** Finding the host elements for markers has failed to converge after 100 iterations. This is usually caused by excessive shear motion of markers or grid boundaries collapsing.
* **Resolution:** Enable grid remeshing or decrease the execution time step to smooth deformation.

### `STOP` — Element Marker Depletion
* **Error Message in stdout:** `No markers in element`
* **Cause:** During the phase ratio calculation (`count_phase_ratio` in `src/marker_data.f90`), the solver detected that an element has no markers inside it (`nmark_elem(j,i) == 0`). This occurs when large material shear or grid deformation causes all markers to move out of a single element, or if the initial marker count was configured too low.
* **Resolution:** 
  1. Increase the initial marker density in the `.inp` file (e.g., the number of markers per element).
  2. If this occurs mid-simulation, configure grid remeshing (`ny_rem = 1` or lower threshold) to reset the grid before elements become too distorted.
  3. Reduce the time step `dt_scale` to stabilize marker advection.

---

## 5. General Execution & CLI Aborts

These errors relate to solver invocation, environment setup, or unsupported configuration boundaries.

### `STOP 1` — Command-line Usage or Geotherm Profile Configuration Error
* **Cause 1: Command-line Usage Error**
  * **Error Message in stdout:** `usage: flac inputfile`
  * **Triggering File:** `src/par.f90`
  * **Details:** You executed the solver without providing exactly one command-line argument for the input file.
  * **Resolution:** Execute the solver by passing the path to a valid `.inp` input file as the single argument (e.g., `./flac model.inp`).
* **Cause 2: Unsupported Geotherm Profile**
  * **Triggering File:** `src/init_temp.f90`
  * **Details:** The geothermal boundary configuration parameter `ictherm` in the `.inp` file is set to an unsupported integer code.
  * **Resolution:** Under the temperature profiles section in your `.inp` file, set `ictherm` to one of the supported options:
    * `1`: Half-space cooling model (oceanic geotherm).
    * `2`: Plate cooling model (oceanic geotherm).
    * `12`: Plate cooling model with radiogenic heating (continental geotherm).
    * `21`: Constant geotherm gradient at top layer, then constant `t_bot`.
    * `22`: Constant geotherm gradient at top two layers, then constant `t_bot`.

### Solver Hangs or Takes Too Long
* **Cause:** The grid size is very large (e.g. 500x200), and the solver is running on a single CPU core.
* **Resolution:** Ensure OpenMP is enabled by setting `OMP_NUM_THREADS` in your shell to use all available cores:
  ```bash
  export OMP_NUM_THREADS=$(nproc)
  ```
