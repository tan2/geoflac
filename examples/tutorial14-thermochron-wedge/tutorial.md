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

### Geometry, Mesh, and Time
* **Dimensions**: $250\text{ km}$ wide ($X \in [0, 250]\text{ km}$), $30\text{ km}$ deep ($Z \in [-30, 0]\text{ km}$).
* **Grid Resolution**: $125 \times 30$ elements ($126 \times 31$ nodes), configured in [`thermochron_wedge.inp`](thermochron_wedge.inp).
* **Grid Type**: Uniform regular grid.
* **Model Time**: $6.5\text{ Myr}$ ($6501.0\text{ Kyr}$ total duration).

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

## 5. Topographic Surface Processes (Erosion & Sedimentation)

During mountain building, rapid tectonic uplift creates steep topography. Without erosion, the growing wedge would eventually become too heavy, preventing further deformation. In this model, surface processes are simulated via a hillslope diffusion equation on the free surface:

$$\frac{\partial z_s}{\partial t} = \frac{\partial}{\partial x} \left( \kappa \frac{\partial z_s}{\partial x} \right)$$

Where:
* $z_s(x)$ is the surface node elevation (`cord(1, :, 2)`).
* $\kappa$ is the topography diffusivity (`topo_removal_rate = 1.0e-6` $\text{m}^2/\text{s}$ or $\approx 31.5\text{ m}^2/\text{yr}$ in [thermochron_wedge.inp](thermochron_wedge.inp)).

### Numerical Formulation
The second derivative of surface elevation is calculated using a central finite-difference scheme (discretized in [fl_move.f90](../../src/fl_move.f90)):

$$dtopo_i = 0.5 \times dt \times \frac{\kappa_{i+1} S_{i+1/2} - \kappa_{i-1} S_{i-1/2}}{x_{i+1} - x_{i-1}}$$

Where $S_{i+1/2}$ is the local slope. Mountain peaks (convex curvature) experience **erosion** ($dtopo_i < 0$), while valleys and basins (concave curvature) experience **sedimentation** ($dtopo_i > 0$). To prevent numerical instability, the erosion depth is capped at half the top-row element height.

### Lagrangian Marker Adjustment (Resurfacing)
Because material properties are carried by Lagrangian markers, the solver runs a `resurface` routine (every `ifreq_avgsr = 10` steps):
* **Erosion**: Markers that end up above the new eroded surface are deleted.
* **Sedimentation**: The new space created by deposition is filled with new markers assigned to the **unconsolidated sediment phase** (`ksed1 = 1`).

This mass transfer unloading at the surface drives rapid exhumation, which directly influences the cooling paths and closure ages tracked by the thermochronology module.

---

## 6. Offline Thermochronology

Unlike other setups where thermochronological ages are computed inline during the simulation, this tutorial computes closure ages **offline** in a post-processing step:
* During the solver run, the node temperatures and the cooling rates (`coolingrate.0`, representing $dT/dt$) are saved.
* A Python post-processing script, [`flacmarker2vtk.py`](../../util/flacmarker2vtk.py), reads the marker coordinate files, the cooling rates, and a reference parameter database [`thermo_chron.dat`](../../util/thermo_chron.dat).
* It computes the Zircon Fission Track (ZFT), Zircon (U-Th)/He (ZHe), Apatite Fission Track (AFT), and other closure ages directly on the Lagrangian markers (which advect with the rock material).
* The output thermochronology ages are saved directly as point data arrays (`age_ZFT`, `age_ZHe`, `age_AFT`, etc.) in the marker VTP files (`flacmarker.*.vtp`) for visualization.

---

## 7. Running the Simulation and Plotting

### Step 1: Run the Solver
Run the FLAC solver on the input file:
```bash
../../src/flac thermochron_wedge.inp
```
The solver will output binary files (such as `phase.0`, `temp.0`, `coolingrate.0`, `vel.0`) at regular time intervals.

### Step 2: Compute Thermochronological Ages
Run the marker post-processing script with the `-t` option to compute offline thermochronology and generate `.vtp` visualization files for the markers:
```bash
python3 ../../util/flacmarker2vtk.py -t .
```
This script reads `coolingrate.0`, the markers, and `util/thermo_chron.dat` to compute ZFT, ZHe, and AFT ages directly on the markers, writing them as point arrays (`age_ZFT`, `age_ZHe`, `age_AFT`, etc.) into `flacmarker.*.vtp`.

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

## 8. Reference Literature
For the physical details, parameter studies, and application to the Taiwan Orogeny, please refer to:
* **Tan et al. (2024)**, *Mountain building process of the Taiwan orogeny*, Science Advances, 10, eadp8056. [PDF reference](file:///home/tan2/Dropbox/Papers/2024/Science%20Advances/Tan%20et%20al-2024-Mountain%20building%20process%20of%20the%20Taiwan%20orogeny.pdf)
