# Troubleshooting GeoFLAC Runtime Aborts

This guide lists the common runtime abort codes (`STOP` statements) that you may encounter when executing simulations, their causes, and how to resolve them.

---

## 1. Grid Geometry & Zone Configuration Aborts

These errors occur during initialization inside the coordinate setup module (`src/init_cord.f90`).

### `STOP 16` — Number of zones is not odd (X-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Number of zones is not odd.(X-direction)`
* **Cause:** The number of grid zones in the X-direction (`nzonx` parameter in the `.inp` file) is an even number (e.g., `2` or `4`).
* **Why this is a constraint:** The mesh generation algorithm (`mesh1`) requires an odd number of zones. The odd-numbered zones (1, 3, 5) define regions of constant spacing, while the even-numbered zones (2, 4) act as transitional regions between them. Having an even number of zones causes memory out-of-bounds access.
* **Resolution:** Set `nzonx` to an odd number (e.g., `1`, `3`, `5`) or `0` for a regular grid.

### `STOP 17` — Nelem in zones is not correct (X-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Nelem in zones is not correct.(X-direction)`
* **Cause:** The sum of elements across all X-direction zones (`nelz-x`) does not match the total elements specified on the mesh parameters line (`nex`).
* **Resolution:** Adjust the number of elements in each zone so that their sum exactly equals the first integer on the Mesh Parameters line (`nex = nx - 1`).

### `STOP 18` — Number of zones is not odd (Z-direction / Y-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Number of zones is not odd.(Y-direction)`
* **Cause:** The number of grid zones in the Z-direction (`nzony` parameter in the `.inp` file) is an even number (e.g., `2` or `4`).
* **Resolution:** Change `nzony` to an odd number (e.g., `1`, `3`, `5`) or `0` for a regular grid.

### `STOP 19` — Nelem in zones is not correct (Z-direction / Y-direction)
* **Error Message in `sys.msg`:** `INIT_CORD: Nelem in zones is not correct.(Y-direction)`
* **Cause:** The sum of elements across all Z-direction zones (`nelz-y`) does not match the total elements specified on the mesh parameters line (`nez`).
* **Resolution:** Adjust the number of elements in each zone so that their sum exactly equals the second integer on the Mesh Parameters line (`nez = nz - 1`).

---

## 2. File & Path Aborts

### `STOP 21` — Coordinates File Error
* **Error Message in `sys.msg`:** `INIT_CORD: Cannot open file with initial coordinates!` or `INIT_CORD: Error reading file with initial coordinates!`
* **Cause:** The parameter `ircoord` is set to `1` (which tells the solver to read coordinates from an external file), but the file specified by `coordfile` (e.g., `coord.dat` or `points.xy`) is missing or cannot be opened.
* **Resolution:** 
  1. Ensure that the coordinates file is in your running directory.
  2. Or, set `ircoord` to `0` to have the solver automatically generate the mesh coordinates.

---

## 3. General Execution Issues

### Solver Hangs or Takes Too Long
* **Cause:** The grid size is very large (e.g. 500x200), and the solver is running on a single CPU core.
* **Resolution:** Ensure OpenMP is enabled by setting `OMP_NUM_THREADS` in your shell to use all available cores:
  ```bash
  export OMP_NUM_THREADS=$(nproc)
  ```
