# GeoFLAC Tutorial: Ocean-Ocean Subduction & Mantle Wedge Melting

This tutorial explains the geodynamical principles, model setup, boundary conditions, dynamic mineral phase transitions, magma generation/migration, and visualization of a 2D **Ocean-Ocean Subduction** simulation in **GeoFLAC**.

---

## 1. Physical and Geodynamical Principles

Subduction zones are the primary drivers of plate tectonics and volcanic arc systems. An oceanic plate subducts into the mantle due to negative buoyancy (slab pull). As the slab descends, it undergoes substantial thermomechanical and mineralogical changes:

```
                  Volcanic Arc
                      |
                      v
      Incoming        |   Overriding Plate
    Oceanic Plate     |       |
  =================   v       v
  \   \   \   \   \=======+==============
   \   \ Basalt    \  Serp.|  Mantle Wedge
    \   \ (Phase 3) \(Ph 9)| (Phase 4/8)
     \   \           \~~~~~+
      \   v           \ H.M.|  Magma Migration
       \ Eclogite      \(Ph16) (Diking to Arc)
        \(Phase 13)     \ ^ |  
         \               \| | 
          v               v v 
```

### Key Processes Modeled:
1. **Slab Pull (Eclogitization)**: As the subducting oceanic crust (basalt, Phase 3/7) goes deeper, pressure and temperature increase. This induces a phase transition to **Eclogite** (Phase 13), which is significantly denser ($\rho = 3480\text{ kg/m}^3$ vs $2880\text{ kg/m}^3$), providing the slab pull force that stabilizes subduction.
2. **Slab Dehydration & Mantle Hydration**: Dehydrating subducted crust releases water into the overlying mantle wedge (olivine, Phase 4/8). At shallow depths ($< 65\text{ km}$), this hydrates the mantle wedge to form **Serpentinite** (Phase 9), a very weak rock ($\phi = 3^\circ$, $c = 4\text{ MPa}$) that decouples the plates.
3. **Mantle Wedge Melting**: At greater depths ($> 65\text{ km}$), water hydrates the hot mantle wedge to form **Hydrated Mantle** (Phase 16). The presence of water lowers the solidus (melting temperature). When temperatures exceed this wet solidus, the hydrated mantle wedge undergoes partial melting.
4. **Magma Diking & Volcanism**: The generated magma collects and migrates vertically through the crust (via a parameterized diking channel, Phase 14) to build a volcanic arc at the surface.

---

## 2. Model Setup

The model represents a regional 2D cross-section of the upper mantle and crust down to $300\text{ km}$ depth.

### Geometry and Mesh
* **Dimensions**: $960\text{ km}$ wide ($X \in [0, 960]\text{ km}$), $300\text{ km}$ deep ($Z \in [-300, 0]\text{ km}$).
* **Grid Resolution**: $185 \times 64$ elements ($186 \times 65$ nodes), configured in [`subduction.inp`](subduction.inp).
* **Non-Uniform Grid Zones**:
  To resolve the subduction zone and arc with high detail while saving computation time elsewhere, a variable grid is defined:
  * **X direction (5 zones)**: High resolution ($2\text{ km}$ elements) in the center ($X \in [300, 600]\text{ km}$) where subduction occurs, and coarser resolution ($6\text{ km}$ elements) near the left and right boundaries.
  * **Z direction (5 zones)**: High resolution ($1.5\text{ km}$ elements) near the surface ($Z \in [-50, 0]\text{ km}$) to capture crustal faults and the volcanic arc, and coarser resolution ($9\text{ km}$ elements) at depth ($Z < -150\text{ km}$).
  
  Review the zone configurations under the mesh section in [`subduction.inp`](subduction.inp) and compare them with the [`nzonx` and `nzony` descriptions](../doc/input_description.md#mesh-parameters).

---

## 3. Initial Thermal and Phase Structure

Because subduction requires a pre-existing slab to start, we initialize the model with a dipping slab and a realistic thermal profile:

### 1. Lithospheric Geotherm
We divide the domain into 5 thermal zones using [`nzone_age`](../doc/input_description.md#initial-structure):
* **Left Plate (Nodes 1 to 58)**: Young oceanic plate (thermal age $100\text{ Myr}$) initialized with a cooling half-space geotherm.
* **Wedge and Right Plate (Nodes 120 to 186)**: Older, stable plate (thermal age $200\text{ Myr}$).
* **Transition Zone (Nodes 59 to 119)**: Linear transition of thermal age to smooth the geotherm.

### 2. Dipping Slab Initialization (Inhomogeneities)
We insert slanting inclusions to model the initial dipping subducted oceanic crust and slab core:
* **Slab Crust (Phase 3)**: A $7.5\text{ km}$ thick slanting basalt layer dipping at $\sim 45^\circ$, initialized using slanting inhomogeneities (`geometry = 4`).
* **Weak Seed (Phase 16)**: A weak zone immediately above the slab to decouple the plates and initialize the shear zone (`init.pl.strain = 1.0`).
* **Cold Slab Core**: A slanting temperature anomaly (`geometry = 13`, `amp = -500°C`) is applied to model the cold core of the subducting lithosphere, preventing the slab from warming up and weakening too quickly.

Review these settings under `inhom` in [`subduction.inp`](subduction.inp) and check their syntax in the [Initial Inhomogeneities section](../doc/input_description.md#initial-inhomogeneities).

---

## 4. Boundary Conditions

The mechanical boundary conditions compress the domain to drive subduction:

1. **Left Boundary (Side 1)**: Inflow velocity of $V_x = 4.9\text{ cm/yr} \approx 1.56 \times 10^{-9}\text{ m/s}$ (`nbc = 10`), driving the subducting plate.
2. **Right Boundary (Side 3)**: Locked horizontal velocity $V_x = 0.0$ (`nbc = 10`), representing a rigid backstop for the overriding plate.
3. **Bottom Boundary (Side 2)**: Hydrostatic pressure support (`nyhydro = 2`, `iphsub = 4`) balancing buoyancy and allowing material to flow in/out of the bottom boundary as the slab subducts.
4. **Top Boundary (Side 4)**: Free surface with topography diffusion (`topo_removal_rate = 1.0e-6`), simulating surface erosion and sedimentation.

---

## 5. Mantle Wedge Melting & Magma Diking

Magma melting and transport are enabled in [`subduction.inp`](subduction.inp) using the following parameters:
* **`itype_melting = 1`**: Enables wet melting of the hydrated mantle wedge.
* **`nelem_serp = 2`**: Specifies a $2$-element-thick serpentinized layer above the subducting slab.
* **`nelem_dike = 1`**: Specifies that magma migrates upward through a $1$-element-wide vertical conduit.
* **`prod_magma = 2e-15`**: Magma extraction rate.
* **`weaken_ratio_plastic = 1.0`, `weaken_ratio_viscous = 1.0`**: Stresses and viscosities within magma-filled cells are weakened proportionally to the magma fraction.

For details on melting mathematics, see the [Magma Parameters](../doc/input_description.md#misc-parameters) description.

---

## 6. Remeshing

As the slab subducts deep into the mantle, elements around the subduction interface undergo severe shear strain and flatten. To prevent the simulation from crashing due to grid distortion:
* **`ny_rem = 1`**: Activates the automatic remeshing engine.
* **`mode_rem = 3`**: Restores the Left, Right, and Bottom boundaries to their initial vertical/horizontal walls while preserving the top free-surface topography.
* **`ntest_rem = 500`**: Checks the grid angles every 500 time steps.
* **`angle_rem = 5.`**: Triggers a rebuild if any grid element is distorted by more than $5^\circ$.

---

## 7. Running the Simulation and Plotting

### Step 1: Run the Solver
Because the 24 Myr subduction simulation is computationally heavy ($\sim 4.8$ million steps), we run a short validation case of **5.0 Myr**:
```bash
export OMP_NUM_THREADS=15
../../src/flac subduction.inp
```
The solver will output binary files (e.g. `phase.0`, `temp.0`, `fmagma.0`, `vel.0`) every 200 kyrs.

### Step 2: Generate the Visualizations
Run the provided Python plotting script:
```bash
python3 plot_subduction.py
```
This script reads the binary outputs and generates three publication-ready figures:
1. **`subduction_full_zone.png`**: A full 2D cross-section showing the subducting slab, phase distribution, temperature isotherms, and plate velocity field.
2. **`subduction_mantle_wedge.png`**: A two-panel zoomed-in view of the mantle wedge ($X \in [300, 700]\text{ km}$, $Z \in [-150, 0]\text{ km}$) displaying:
   * **Top Panel**: Active phase boundaries (serpentinite decoupling, hydrated mantle wedge) and velocity vectors.
   * **Bottom Panel**: Magma/melt generation zone in the mantle wedge beneath the volcanic arc.
3. **`plots/subduction_evolution.png`**: A three-panel evolutionary sequence showing how the slab sinks over time.

---

## 8. Analyzing the Results

### Slab Eclogitization and Angle
In the full profile (`subduction_full_zone.png`), watch how the subducting slab (green olivine mantle + blue oceanic crust) sinks. As the blue oceanic crust descends below $\sim 60\text{ km}$, it transforms into dark red **Eclogite** (Phase 13). The high density of eclogite drives the slab down, steepening the subduction angle.

### Mantle Wedge Flow and Hydration
Observe the velocity vectors in the mantle wedge. The descending slab drags the adjacent wedge mantle downwards, establishing a **corner flow** convection cell. 
Directly above the slab, you will see a thin yellow layer of **Serpentinite** (Phase 9) at shallow depths, which transitions into **Hydrated Mantle** (Phase 16) at greater depths.

### Arc Magmatism
Look at the zoomed-in plot (`subduction_mantle_wedge.png`). The hydrated mantle (Phase 16) sits in the hot core of the mantle wedge ($T > 800^\circ\text{C}$). Because of hydration, this region starts melting, visible as a red-colored anomaly of high **magma fraction**. The extracted magma travels vertically through the diking channel (Phase 14) directly above the melting zone to form the volcanic arc at the surface.
