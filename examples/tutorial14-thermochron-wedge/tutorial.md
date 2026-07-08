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
   \  (Crustal Shortening)                 |  (Rigid)
    \                                      |
     \___ (Flat 1)                         |
         \___ (Flat 2)                     |
             \___ (Flat 3)                 |
                 \________________________/|
             Basal Detachment (Ramp-and-Flat)
```

---

## 2. Model Setup

The model represents a regional cross-section of the orogenic wedge down to $30\text{ km}$ depth.

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

### 2. Horizontal Zoning and Varying Layer Thicknesses (`nzone_age`)
To represent the complex lateral transitions of a continental margin, the domain is divided horizontally into 18 zones using the `nzone_age` input block. In the input file, each zone is defined by three lines of parameters:

**Line 1: Thermal model and horizontal bounds**
```
ictherm   age_1   tp1   tp2   ixtb1   ixtb2
```
* **`ictherm`**: Geotherm model selector (e.g., `1` for oceanic half-space cooling, `12` for continental plate cooling).
* **`age_1`**: Thermal age of the lithosphere [Myr].
* **`tp1`, `tp2`**: Model-dependent constants (typically `0.0`).
* **`ixtb1`, `ixtb2`**: Horizontal node index range defining the zone bounds.

**Line 2: Layer thicknesses**
```
nph_layer   hc(1)   hc(2)   ...   hc(N-1)
```
* **`nph_layer`**: Total number of lithological layers in the column.
* **`hc(j)`**: Thickness of layer $j$ [km] (the bottom layer $N$'s thickness is auto-calculated to fit the domain depth).

**Line 3: Layer phases**
```
iph_col(1)   iph_col(2)   ...   iph_col(N)
```
* **`iph_col(j)`**: Lithological phase IDs of the layers from top to bottom.

#### Linear Interpolation of Transitional Thicknesses
A **smooth transition** between zone $k$ and zone $k+1$ is established by setting `ixtb2(k) = -1` and `ixtb1(k+1) = -1`. The solver will then linearly interpolate the layer boundary depths and thermal ages between node `ixtb1(k)` and node `ixtb2(k+1)`.

**Example: Tapering Continental margin**
```
1, 10.0, 0.0, 0.0, 20, -1 
   5, 4.0 4.0 4.0 5.0
   13 14 16 18 19
1, 10.0, 0.0, 0.0, -1, 25 
   5, 8.0 8.0 8.0 9.0
   13 14 16 18 19
```
* In this setup, the transition spans from horizontal node index 20 to 25.
* At node 20, the upper crustal layers (Phases 13, 14, 16) are each $4.0\text{ km}$ thick.
* At node 25, the layers thicken to $8.0\text{ km}$.
* The solver linearly interpolates the thickness of each layer for the intermediate nodes (21 to 24), simulating a smooth, realistic thickening of the continental crust margin.

### 3. Layers and Inhomogeneities
Initial crustal layers and local geometric structures (such as basal detachment ramps or sedimentary basins) are defined using inhomogeneities. In the input file [thermochron_wedge.inp](thermochron_wedge.inp), inhomogeneities are set using lines of numeric parameters:

```
ix1   ix2   iy1   iy2   inphase   igeom   xinitaps
```

Where the parameters represent:
* **`ix1`, `ix2`**: Element column index range in the horizontal direction ($X$).
* **`iy1`, `iy2`**: Element row index range in the vertical direction ($Z$).
* **`inphase`**: Rock type phase ID to assign. If `inphase > 0`, the elements are overwritten with this phase. If `inphase = -1`, the existing phase is preserved, allowing local weakening without changing the lithology.
* **`igeom`**: Geometry selector:
  * `0`: Rectangular block defined by the upper-left corner (`ix1`, `iy1`) and lower-right corner (`ix2`, `iy2`).
  * `3`: Diagonal line segment connecting (`ix1`, `iy1`) to (`ix2`, `iy2`).
  * `4`: Diagonal line segment with pre-existing plastic shear strain.
* **`xinitaps`**: Initial plastic strain value (only used when `igeom = 4` to set the weakness of the zone).

#### Concrete Setup Examples

1. **Defining a Pre-weakened Diagonal Ramp (Basal Detachment)**:
   ```
   19  24  4  8  -1  4  1.0
   ```
   * This defines a slanting line segment spanning columns 19 to 24, dipping from row 4 to row 8.
   * `inphase = -1` keeps the original rock phases.
   * `igeom = 4` and `xinitaps = 1.0` initializes the plastic strain in this zone to $1.0$, creating a pre-existing shear zone that acts as a weak thrust fault guide.

2. **Defining a Flat Crustal Layer Boundary**:
   ```
   46  115  12  12  15  4  0.0
   ```
   * This defines a horizontal segment (`iy1 = iy2 = 12`) from column 46 to 115.
   * `inphase = 15` changes the elements along this boundary to Phase 15 (metasediment).
   * `igeom = 4` specifies a line inclusion, and `xinitaps = 0.0` leaves the initial plastic strain at 0.


---

## 4. Boundary Conditions & Nodal Flexure

The mechanical boundary conditions simulate the tectonic convergence and lithospheric flexural loading:
1. **Left Boundary (Side 1)**: Inflow velocity representing the incoming continental plate.
2. **Right Boundary (Side 3)**: Locked vertical/horizontal velocities representing the collision backstop.
3. **Top Boundary (Side 4)**: Free surface with topography diffusion (`topo_removal_rate = 1.0e-6`).
4. **Bottom Boundary (Side 2)**: Flexural boundary condition (`nbc = 200`) simulating the elastic bending of the underlying lithosphere under tectonic loads.

### Multi-Segment Velocity BC (Side 3)
To simulate bivergent continental convergence, the right boundary (Side 3) is divided into multiple velocity segments that define an S-point (velocity discontinuity):
*   **Nodes 1 to 20 (Upper segment)**: Moves horizontally to the left at a constant velocity of $-4.65 \times 10^{-10}\text{ m/s}$ ($\approx -1.47\text{ cm/yr}$), representing the tectonic push of the overriding plate.
*   **Nodes 23 to 31 (Lower segment)**: Locked horizontally ($v_x = 0.0\text{ m/s}$), representing the stationary underthrusting lithosphere.
*   **Nodes 21 to 22 (Transition zone)**: The horizontal velocity is linearly interpolated ($v_x = a + bx$) from the upper segment velocity to $0.0\text{ m/s}$. This transitional segment prevents numerical stress singularities at the velocity discontinuity and guides the initiation of the basal detachment thrust.

### Nodal Flexure (Buoyancy Support)
As mountains grow, their massive weight pushes down on the Earth's outermost rigid shell (the lithosphere). To represent this in the model, the bottom boundary acts as a flexible, elastic sheet floating on a thick, buoyant fluid (representing the mantle asthenosphere). 

When the tectonic wedge compresses and thickens, the weight of the crust increases. The bottom boundary responds by bending downwards under this extra load. Conversely, if mountain peaks erode and lighten, the boundary rebounds upwards. 

Crucially, when the bottom boundary bends, the entire vertical crustal column above it translates downwards or upwards together with it. This ensures that the crust behaves as a coherent, elastic sheet and prevents the grid elements from stretching or distorting artificially.

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

### What is Thermochronology? (Non-Technical Explanation)
For geologists, **thermochronology** is a method used to reconstruct the cooling histories of rocks as they travel through the Earth's crust. Unlike traditional geochronology (which measures when a mineral first crystallized), thermochronology measures the time elapsed since a rock cooled past a specific threshold temperature called the **closure temperature** ($T_c$). 

Inside the Earth, temperature increases with depth (geotherm). When a rock is deep and hot ($T > T_c$), isotopic decay products (such as helium gas or fission tracks in apatite and zircon crystals) diffuse out of the crystal structure immediately, meaning the "cooling clock" is open and cannot tick. As tectonic forces uplift the rock and erosion strips away the surface, the rock moves upwards into cooler crustal levels. Once the rock cools past the closure temperature ($T \le T_c$), the crystal lattice locks, traps the decay products, and the cooling clock starts ticking. 

Thus, younger thermochronological ages at the surface indicate rapid tectonic uplift and exhumation, while older ages indicate slow, quiet geological cooling.

### Running the Offline Thermochronology Script
To avoid adding computational overhead to the main Fortran solver, thermochronological ages are tracked on Lagrangian markers and computed in a post-processing step using the following commands:

1. **Convert the structural grid to VTK format** (generates `.vts` grid files):
   ```bash
   python3 ../../util/flac2vtk.py .
   ```
2. **Compute thermochronology ages on the markers** (generates `.vtp` marker files):
   ```bash
   python3 ../../util/flacmarker2vtk.py -t .
   ```
   * *Note*: The `-t` flag instructs the script to calculate offline cooling closure ages by reading the cooling rates (`coolingrate.0`) and the reference closure database [`thermo_chron.dat`](../../util/thermo_chron.dat).
3. **Interpolate marker ages onto the structural grid** (adds age fields to `.vts` files):
   ```bash
   python3 ../../util/tcages2grid.py .
   ```
   * *Note*: This maps the tracked marker ages onto the Lagrangian grid element cells, making them available as element datasets in the VTS structural grid files.

### Understanding Special Age Values in the Output
When visualizing the marker age datasets in VisIt or ParaView, you will find three distinct types of ages:

1. **Positive Ages ($Age \ge 0.0\text{ Myr}$)**:
   * **Geological Meaning**: The rock was once buried deeply (heated past its closure temperature $T_c$ to reset the clock) and has since cooled below $T_c$ due to exhumation. The value represents the number of millions of years ($Myr$) elapsed since the rock crossed the closure temperature boundary on its path to the surface.
2. **`-1.0` (Open/Unclosed System)**:
   * **Geological Meaning**: The marker is currently still at depth and is **hotter than the closure temperature** ($T > T_c$). Helium gas is still diffusing out of the crystal, so the cooling clock is open and has not started ticking.
3. **`NaN` (Unreset System)**:
   * **Geological Meaning**: The marker has **never been heated hot enough** during the experiment to exceed the closure temperature (always $T \le T_c$). Consequently, its pre-depositional geologic age has never been reset. This typically occurs for surface sediments or shallow crustal rocks that were never buried deeply.

---

## 7. Running the Simulation and Plotting

### Step 1: Run the Solver
Run the FLAC solver on the input file:
```bash
../../src/flac thermochron_wedge.inp
```
The solver will output binary files (such as `phase.0`, `temp.0`, `coolingrate.0`, `vel.0`) at regular time intervals.

### Step 2: Compute and Interpolate Thermochronological Ages
1. Run the marker post-processing script with the `-t` option to compute offline thermochronology and generate `.vtp` files for the markers:
   ```bash
   python3 ../../util/flacmarker2vtk.py -t .
   ```
   This script reads `coolingrate.0`, the markers, and `util/thermo_chron.dat` to compute ZFT, ZHe, and AFT ages directly on the markers, writing them as point arrays (`age_ZFT`, `age_ZHe`, `age_AFT`, etc.) into `flacmarker.*.vtp`.
2. Interpolate the computed marker ages onto the structural grid to update `.vts` files:
   ```bash
   python3 ../../util/tcages2grid.py .
   ```
   This maps the cooling ages from the marker points onto the Lagrangian mesh elements, allowing you to plot and analyze thermochronology as grid-wide fields.

### Step 3: Plot results
Generate diagnostic figures showing the shear zones and phases:
```bash
python3 plot_thermochron_wedge.py
```
This script creates:
1. **`thermochron_wedge.png`**: A two-panel plot displaying:
   * **Top Panel**: Accumulated plastic strain (shear zones) and temperature isotherms.
   * **Bottom Panel**: Deformed lithological phases and grid mesh.
2. **`images/evolution_thermochron_wedge.png`**: The evolution of shear localization and faults over time.

---

## 8. Reference Literature
For the physical details, parameter studies, and application to the Taiwan Orogeny, please refer to:
* **Tan et al. (2024)**, *Mountain building process of the Taiwan orogeny*, Science Advances, 10, eadp8056. [PDF reference](file:///home/tan2/Dropbox/Papers/2024/Science%20Advances/Tan%20et%20al-2024-Mountain%20building%20process%20of%20the%20Taiwan%20orogeny.pdf)

---

## Appendix: Mathematical and Numerical Formulation of the Flexural BC

### Governing Equation
The vertical deflection $w(x)$ of the bottom boundary under a vertical load $q(x)$ is governed by the 1D thin elastic plate flexure equation:

$$\frac{d^2}{dx^2} \left( D(x) \frac{d^2 w}{dx^2} \right) + \rho_m g w = q(x) - q_{init}(x)$$

Where:
* $w(x)$ is the downward vertical deflection [m].
* $D(x) = \frac{E T_e^3}{12 (1 - \nu^2)}$ is the flexural rigidity [N$\cdot$m], calculated dynamically using the Lamé parameters of the bottom elements.
* $T_e$ is the elastic thickness of the plate (`Te = bca = 10000.0` meters, or $10\text{ km}$).
* $\rho_m$ is the underlying asthenosphere density (`rho_m = bcb = 3300.0` $\text{kg/m}^3$).
* $g$ is gravitational acceleration ($10.0\text{ m/s}^2$).
* $\rho_m g w$ is the Winkler foundation buoyancy restoring force representing asthenosphere displacement.
* $q(x) - q_{init}(x)$ is the excess load [Pa] representing the change in vertical lithostatic column weight relative to the initial state.

### Numerical Implementation (Tridiagonal Solver & Whole Column Coupling)
1. **Block Tridiagonal Formulation**: The fourth-order flexure equation is split into two coupled second-order differential equations in terms of vertical deflection $w$ and bending moment $M = D \frac{d^2 w}{dx^2}$. These equations are discretized in [fl_flexure.f90](../../src/fl_flexure.f90) using a second-order non-uniform finite-difference scheme, yielding a block tridiagonal linear system solved on the GPU.
2. **Elastic Column Translation**: In each step, if a column $i$ undergoes an incremental vertical deflection $\Delta w_i = w_i(t) - w_i(t - dt)$:
   * The coordinates of all nodes in that vertical column ($j = 1$ to $nz$) are shifted downward by $\Delta w_i$ to keep element thicknesses unchanged:
     $$z_{j,i} \leftarrow z_{j,i} - \Delta w_i$$
   * The vertical velocities of the entire column are updated by the flexural deflection velocity:
     $$v_{z,\, j,i} \leftarrow v_{z,\, j,i} - \frac{\Delta w_i}{dt}$$
   * To prevent velocity accumulation/inflation across timesteps, the previous step's flexural velocity is subtracted at the start of the next flexure step.
