# GeoFLAC Tutorial: Ocean-Ocean Subduction & Mantle Wedge Melting

This tutorial explains the geodynamical principles, model setup, boundary conditions, dynamic mineral phase transitions, magma generation/migration, and visualization of a 2D **Ocean-Ocean Subduction** simulation in **GeoFLAC**.

---

## 1. Physical and Geodynamical Principles

Subduction zones are the primary drivers of plate tectonics and volcanic arc systems. An oceanic plate subducts into the mantle due to negative buoyancy (slab pull). As the slab descends, it undergoes substantial thermomechanical and mineralogical changes:

```
                                Trench                  Arc
                                  |                     Volcano
                                  |
  Subducting                      V                    +---+   Overriding 
  Plate (Oceanic)                                     /     \  Plate (Oceanic)
  =======================================================================
  Crust (Basalt, Ph 3)             \ |   Crust (Basalt, Ph 3)
  -----------------------+          \+---------------------------
  Lithospheric Mantle     \ Basalt   \  \<-- Serp. (Ph 9)  
  (Phase 4)                \ (Ph 3)   \  \ 
                            \          \  \           Mantle Wedge
                             \        ..\  \         (Phase 4/8)
                              \    ...   \ |
                               \...       \|       Magma Migration
                                \          \             ^
                                 \ Eclogite \            |
                                  \ (Ph 13)  \           |
                                   v          v
```

### Key Processes Modeled:
1. **Slab Pull (Eclogitization)**: As the subducting oceanic crust (basalt, Phase 3/7) goes deeper, pressure and temperature increase. This induces a phase transition to **Eclogite** (Phase 13), which is significantly denser ($\rho = 3480\text{ kg/m}^3$ vs $2880\text{ kg/m}^3$), providing the slab pull force that stabilizes subduction.
2. **Slab Dehydration & Mantle Hydration**: Dehydrating subducted crust releases water into the overlying mantle wedge (olivine, Phase 4/8). At shallow depths ($< 65\text{ km}$), this hydrates the mantle wedge to form **Serpentinite** (Phase 9), a very weak rock ($\phi = 3^\circ$, $c = 4\text{ MPa}$) that decouples the plates.
3. **Mantle Wedge Melting**: Phase 16 represents chlorite-bearing hydrated mantle in the subducting lithosphere. As the slab descends, Phase 16 undergoes dehydration breakdown. The released water migrates upward to wet the hot mantle wedge, significantly lowering its solidus (melting temperature) and inducing wet partial melting.
4. **Magma Dynamics & Crustal Processing**:
   * **a. Magma Migration Modeling**: When partial melting occurs in the mantle wedge, it will produce `prod-magma` amount of magma per unit of time per unit of partial melting. Parts (`ratio-mantle-mzone`) of the magma is extracted and evenly distributed in an invertical triangle, whose tip is in the partial melting element and base on the Moho. The angle of the triangle is controlled by the parameter `angle-mzone`. The rest of the magma is distributed within the crust, forming vertical dikes with `nelem-dike` element thick. The volume of these magma will form arc crust (Phase 14) at the surface.
   * **b. Weakening Effects of Magma**: The presence of magma reduces local rock strength. Magma-bearing elements undergo severe rheological weakening (lowered viscosity `weaken-ratio-viscous`, reduced yield stress `weaken-ratio-plastic,`).
   * **c. Magma Freezing (Crystallization) Modeling**: As ascending magma enters colder environments ($T < T_{\text{solidus}}$), it cools and freezes. The freezing rate is controlled with `lambda-freeze` and its temperature dependence is `lambda-freeze-tdep`. During freezing, the latant heat is released, which will increase the local tempeature.

---

## 2. Model Setup

The model represents a regional 2D cross-section of the upper mantle and crust down to $200\text{ km}$ depth and $600\text{ km}$ wide.

### Geometry and Mesh
* **Dimensions**: $600\text{ km}$ wide ($X \in [0, 600]\text{ km}$), $200\text{ km}$ deep ($Z \in [-200, 0]\text{ km}$).
* **Grid Resolution**: $100 \times 40$ elements ($101 \times 41$ nodes), defined in [`subduction.inp`](subduction.inp).
* **Non-Uniform Grid Zones**:
  To resolve the subduction zone and arc with high detail while saving computation time elsewhere, a variable grid is defined:
  * **X direction (5 zones)**: High resolution ($\sim 2.7\text{ km}$ elements) in the center ($X \in [200, 450]\text{ km}$) where subduction occurs, and coarser resolution ($\sim 11\text{ km}$ elements) near the left and right boundaries.
  * **Z direction (3 zones)**: High resolution near the surface ($Z \in [-50, 0]\text{ km}$) to capture crustal faults and the volcanic arc, and coarser resolution at depth.

---

## 3. Initial Thermal and Phase Structure

Because subduction requires a pre-existing slab to start, we initialize the model with a dipping slab extending down to **$50\text{ km}$ depth** and distinct plate thermal profiles. Both the subducting (left) and overriding (right) plates are modeled as oceanic lithosphere:

### 1. Thermal Age Contrast (`nzone_age`)
We partition the domain horizontally into 2 primary plate thermal age segments using the `nzone_age` input block:
* **Subducting Plate (Nodes 1 to 31)**: Older, colder oceanic plate (thermal age **$140\text{ Myr}$**) initialized with a cooling half-space geotherm (`ictherm = 1`).
  * **Layer Structure**: $1.5\text{ km}$ Sediment (Phase 1), $6.0\text{ km}$ Basaltic Crust (Phase 3), $10.0\text{ km}$ Hydrated Slab Mantle (Phase 16), and Olivine Lithospheric Mantle (Phase 4) below.
* **Overriding Plate (Nodes 32 to 101)**: Younger, warmer oceanic plate (thermal age **$60\text{ Myr}$**).
  * **Layer Structure**: $1.5\text{ km}$ Sediment (Phase 1), $6.0\text{ km}$ Basaltic Crust (Phase 3), and Olivine Mantle (Phase 4) below.

### 2. Pre-existing Slab Inclusions (Extended to $50\text{ km}$ Depth)
To initiate subduction smoothly at $X \approx 200-280\text{ km}$ (node 30), initial heterogeneities are inserted dipping at $\sim 45^\circ$ to the right:
* **Dipping Basaltic Slab Crust (Phase 3)**: Inserted between nodes 30-42 dipping from surface down to `iy2 = 16` (**$50\text{ km}$ depth**).
* **Weak Shear Zone (Phase 16 & -1)**: Hydrated weak shear zone immediately above the dipping slab (nodes 27-45, `iy2 = 16`).
* **Cold Slab Core**: Temperature anomaly (`geometry = 13`, `amp = -500°C`) extending down to `iy2 = 17` to preserve the cold thermal core of the subducting slab.

---

## 4. Mechanical Boundary Conditions

* **Left Boundary (Side 1)**: Pushed rightward at a constant convergence velocity $V_x = +1.56 \times 10^{-9}\text{ m/s}$ ($\sim 5\text{ cm/yr}$).
* **Right Boundary (Side 3)**: Fixed stationary boundary ($V_x = 0$).
* **Bottom Boundary (Side 2)**: Hydrostatic pressure boundary with basal temperature fixed at $1330^\circ\text{C}$.

---

## 5. Running the Simulation

1. **Compile the Solver** (with OpenMP support for speed):
   ```bash
   cd src/
   make clean && make omp=1
   cd ../examples/tutorial13-subduction/
   ```

2. **Execute the Code**:
   ```bash
   ../../src/flac subduction.inp
   ```

---

## 6. Visualizing the Evolution

A python script [`plot_subduction.py`](plot_subduction.py) is provided to visualize the lithological evolution and temperature geotherms of the subduction process.

To generate the plots, execute:
```bash
./plot_subduction.py
```
This reads the simulation output files and saves the lithological phase and temperature profile sequence to `images/subduction_evolution.png`.
