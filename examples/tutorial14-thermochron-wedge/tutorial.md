# GeoFLAC Tutorial: Taiwan Orogenic Wedge & Offline Thermochronology

This tutorial explains the geodynamical principles, model setup, boundary conditions, and offline thermochronological modeling of a 2D **Orogenic Wedge** simulation in **GeoFLAC**, based on the mountain building process of the Taiwan orogeny (Tan et al., 2024, *Science Advances*).

---

## 1. Physical and Geodynamical Principles

Orogenic belts, such as the magmatically quiet Taiwan collision zone, are often studied using the critical wedge theory and bivergent wedge S-point models. In Taiwan, the collision between the Eurasian continental margin and the heading Luzon volcanic arc drives rapid deformation, crustal shortening, and exhumation.

This tutorial implements a thermomechanical wedge model that refines the classic S-point framework by incorporating:
1. **Brittle-Ductile Transition**: Varies with temperature and lithology at depth.
2. **Decollement Geometry**: Features a ramp-and-flat geometry that guides subduction and thrusting of the incoming crust.
3. **Lithology-dependent Surface Erosion**: Simulated via free-surface diffusion.
4. **Offline Thermochronology**: Marker-based tracking of cooling histories to predict low-temperature thermochronometer closure ages (ZFT, ZHe, AFT) offline, avoiding inline calculation overhead.

```
       Hsuehshan Range (HR)        Backbone Range (BR)
               |                            |
               v                            v
  ============+======+=====================+======+======
  \ Accretionary Wedge                     |  Backstop
   \                                       |  (Rigid)
    \  Basal Detachment (Ramp-and-Flat)    |
     \______________________/\____________/|
```

---

## 2. Model Setup

The model represents a regional cross-section of the orogenic wedge down to $30 depth.

### Geometry and Mesh
* **Dimensions**: $250\text{ km}$ wide ($X \in [0, 250]\text{ km}$), $30\text{ km}$ deep ($Z \in [-30, 0]\text{ km}$).
* **Grid Resolution**: $125 \times 30$ elements ($126 \times 31$ nodes), configured in [`thermochron_wedge.inp`](thermochron_wedge.inp).
* **Grid Type**: Uniform regular grid.

---

## 3. Initial Thermal and Phase Structure

Because the orogen grows by compressing a variable continental margin, the model incorporates pre-defined crustal layers and a transition of thermal ages:

### 1. Varying Thermal Age
We divide the domain into 18 zones using `nzone_age` to initialize a geotherm that matches the incoming continental margin and the backstop region.

### 2. Layers and Inhomogeneities
The layers of the continental crust (such as basalt, continental basement, metasediments, and weaker decollement layers) are defined using inhomogeneities with `geometry = 4` (diagonal slab ramp-and-flats) and `geometry = 0` (rectangular blocks).
* **Decollement Layer (Phase 18)**: Sets the low-friction basal detachment channel to guide the wedge growth.
* **Basement and Accretionary wedge (Phases 12, 15, 17)**: Stratigraphic layers undergoing compression.

---

## 4. Boundary Conditions

The mechanical boundary conditions simulate the tectonic convergence:
1. **Left Boundary (Side 1)**: Inflow velocity representing the incoming continental plate.
2. **Right Boundary (Side 3)**: Locked vertical/horizontal velocities representing the collision backstop.
3. **Bottom Boundary (Side 2)**: Hydrostatic pressure support to represent buoyancy.
4. **Top Boundary (Side 4)**: Free surface with topography diffusion (`topo_removal_rate = 1.0e-6`).

---

## 5. Offline Thermochronology

Unlike other setups where thermochronological ages are computed inline during the simulation, this tutorial computes closure ages **offline** in a post-processing step:
* During the solver run, the node temperatures and the cooling rates (`coolingrate.0`, representing $dT/dt$) are saved.
* A Python post-processing script, [`flac2vtk.py`](../../util/flac2vtk.py), reads the marker coordinate files, the cooling rates, and a reference parameter database [`thermo_chron.dat`](../../util/thermo_chron.dat).
* It computes the Zircon Fission Track (ZFT), Zircon (U-Th)/He (ZHe), and Apatite Fission Track (AFT) ages on the Lagrangian markers (which advect with the rock material) and interpolates them back to the grid nodes.
* The output is written directly into structured VTS files (`flac.*.vts`) for visualization.

---

## 6. Running the Simulation and Plotting

### Step 1: Run the Solver
Run the FLAC solver on the input file:
```bash
../../src/flac thermochron_wedge.inp
```
The solver will output binary files (such as `phase.0`, `temp.0`, `coolingrate.0`, `vel.0`) at regular time intervals.

### Step 2: Compute Thermochronological Ages
Run the post-processing script with the `-t` option to compute offline thermochronology and generate `.vts` visualization files:
```bash
python3 ../../util/flac2vtk.py -t .
```
This script reads `coolingrate.0`, the markers, and `util/thermo_chron.dat` to compute ZFT, ZHe, and AFT ages, writing them into `flac.*.vts`.

### Step 3: Plot results
Generate diagnostic figures showing the shear zones and phases:
```bash
python3 plot_thermochron_wedge.py
```
This script creates:
1. **`thermochron_wedge.png`**: A two-panel plot displaying:
   * **Top Panel**: Accumulated plastic strain (shear zones) and temperature isotherms.
   * **Bottom Panel**: Deformed lithological phases and grid mesh.
2. **`plots/evolution_thermochron_wedge.png`**: The evolution of shear localization and faults over time.

---

## 7. Reference Literature
For the physical details, parameter studies, and application to the Taiwan Orogeny, please refer to:
* **Tan et al. (2024)**, *Mountain building process of the Taiwan orogeny*, Science Advances, 10, eadp8056. [PDF reference](file:///home/tan2/Dropbox/Papers/2024/Science%20Advances/Tan%20et%20al-2024-Mountain%20building%20process%20of%20the%20Taiwan%20orogeny.pdf)
