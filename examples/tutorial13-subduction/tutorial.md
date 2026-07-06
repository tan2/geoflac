# GeoFLAC Tutorial: Ocean-Ocean Subduction & Mantle Wedge Melting

This tutorial explains the geodynamical principles, model setup, boundary conditions, dynamic mineral phase transitions, magma generation/migration, and visualization of a 2D **Ocean-Ocean Subduction** simulation in **GeoFLAC**.

---

## 1. Physical and Geodynamical Principles

Subduction zones are the primary drivers of plate tectonics and volcanic arc systems. An oceanic plate subducts into the mantle due to negative buoyancy (slab pull). As the slab descends, it undergoes substantial thermomechanical and mineralogical changes:

```
  Subducting Plate (Oceanic)                   Overriding Plate (Oceanic)
  =======================================================================
  Crust (Basalt, Ph 3)   |                   | Crust (Basalt, Ph 3)
  -----------------------+                   +---------------------------
  Lithospheric Mantle     \  Serp. (Ph 9)   /
  (Phase 4)                \~~~~~~.        /
                            \ H.M. `·.    /         Mantle Wedge
                             \        `·./           (Phase 4/8)
                              \ Basalt   \
                               \ (Ph 3)   \       Magma Migration
                                \          \             ^
                                 \ Eclogite \            |
                                  \ (Ph 13)  \           |
                                   v          v
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

---

## 3. Initial Thermal and Phase Structure

Because subduction requires a pre-existing slab to start, we initialize the model with a dipping slab and a realistic thermal profile. Both the subducting (left) and overriding (right) plates are modeled as oceanic lithosphere:

### 1. Lithospheric Geotherm
We divide the domain into 5 thermal zones using [`nzone_age`](../doc/input_description.md#initial-structure):
* **Subducting Plate (Nodes 1 to 58)**: Young oceanic plate (thermal age $100\text{ Myr}$) initialized with a cooling half-space geotherm.
* **Wedge and Overriding Plate (Nodes 59 to 186)**: Older, stable oceanic plate (thermal age $200\text{ Myr}$), showing a thicker lithospheric geotherm.

### 2. Layer structure
Both plates have a standard 3-layer oceanic crust structure:
* **Sediment (Phase 1)**: $1.5\text{ km}$ thick.
* **Basaltic Crust (Phase 3)**: $6.0\text{ km}$ thick.
* **Hydrated Slab Mantle (Phase 16)**: $10.0\text{ km}$ thick.
* **Lithospheric Mantle (Phase 4)**: Olivine below $17.5\text{ km}$ depth.

### 3. Dipping Slab Initialization (Inhomogeneities)
We insert slanting inclusions to model the initial dipping subducted oceanic crust and slab core:
* **Slab Crust (Phase 3)**: A $7.5\text{ km}$ thick slanting basalt layer dipping at $\sim 45^\circ$, initialized using slanting inhomogeneities (`geometry = 4`).
* **Weak Seed (Phase 16)**: A weak zone immediately above the slab to decouple the plates and initialize the shear zone (`init.pl.strain = 1.0`).
* **Cold Slab Core**: A slanting temperature anomaly (`geometry = 13`, `amp = -500°C`) is applied to model the cold core of the subducting lithosphere, preventing the slab from warming up and weakening too quickly.

---

## 4. Boundary Conditions

The mechanical boundary conditions compress the domain to drive subduction:

1. **Left Boundary (Side 1)**: Inflow velocity of $V_x = 4.9\text{ cm/yr} \approx 1.56 \times 10^{-9}\text{ m/s}$ (`nbc = 10`), driving the subducting plate.
2. **Right Boundary (Side 3)**: Locked horizontal velocity $V_x = 0.0$ (`nbc = 10`), representing a rigid backstop for the overriding plate.
3. **Bottom Boundary (Side 2)**: Hydrostatic pressure support (`nyhydro = 2`, `iphsub = 4`) balancing buoyancy and allowing material to flow in/out of the bottom boundary as the slab subducts.
4. **Top Boundary (Side 4)**: Free surface with topography diffusion (`topo_removal_rate = 1.0e-6`), simulating surface erosion and sedimentation.

---

## 5. Mantle Wedge Melting, Magma Distribution & Freezing

The magma module in GeoFLAC models wet melting, transport, distribution, and freezing of magma:

### 1. Magma Generation & Wet Melting
When `itype_melting = 1` is set in [`subduction.inp`](subduction.inp), wet melting calculations are activated. Melting is computed for cells in the mantle wedge containing water/hydration (specifically **Phase 16, Hydrated Mantle**). The water-saturated solidus is calculated based on Grove et al. (2009) Nature:
- For depths $d < 80\text{ km}$:
  $$T_{solidus} = 800^\circ\text{C}$$
- For depths $d \ge 80\text{ km}$:
  $$T_{solidus} = 800 + 6.2 \times 10^{-8} \cdot (d - 80000)^2\text{ } (^\circ\text{C})$$

If the local temperature $T$ exceeds $T_{solidus}$, the partial melt fraction (magma fraction increment) $F$ generated in that step is:
$$F = \min\left(\frac{C_p \cdot (T - T_{solidus})}{L_{fusion}}, 0.10\right)$$
where $C_p$ is heat capacity and $L_{fusion}$ is latent heat of fusion (`latent_heat_magma = 420000.0` J/kg). The melt generation is capped at a maximum of 10% per step.
To simulate the energy consumed by the phase transition, the local temperature is pulled back to $T_{solidus}$ (latent heat buffering), subtracting the overshoot $T - T_{solidus}$ from the element's corner nodes.

### 2. Magma Distribution & Transport (Diking vs Percolation)
Once generated in the mantle wedge ($j > j_{moho}$), the melt is instantaneously extracted and transported upwards. It is distributed between the mantle wedge and the crust:

#### Mantle Wedge Percolation (Porosity Flow)
* In the mantle wedge (depths below the Moho, $j_{moho} < jj \le j$), magma ascends by porous percolation through a slanted cone centered at the melting point, with a dip defined by `angle_mzone = 30.0` degrees.
* A fraction of the extracted magma (`ratio_mantle_mzone = 0.1`, i.e. 10%) is distributed to the mantle wedge elements within this cone:
  $$\Delta f_{magma}^{mantle} = f_{melt} \cdot \text{ratio\_mantle\_mzone} \cdot A_{ratio} \cdot \text{prod\_magma} \cdot dt$$
  where $A_{ratio}$ is the area ratio of the mantle cone to the melting element, and `prod_magma = 2e-15` is the magma extraction rate.

#### Crustal Diking Conduit
* In the overriding crust (depths above the Moho, $jj \le j_{moho}$), magma ascends by fast vertical diking through a conduit of width `nelem_dike = 1` element, centered horizontally above the melting source.
* The remaining fraction of extracted magma ($1 - \text{ratio\_mantle\_mzone} = 0.9$, i.e. 90%) is transported directly to the crust:
  $$\Delta f_{magma}^{crust} = \frac{f_{melt} \cdot (1 - \text{ratio\_mantle\_mzone}) \cdot A_{ratio} \cdot \text{prod\_magma} \cdot dt}{n_{col}}$$
  where $n_{col}$ is the number of elements in the diking channel.
* The total accumulated magma fraction $f_{magma}$ in any element is capped at `fmagma_max = 0.1` (10%).

### 3. Freezing and Crystallization
As the magma ascends to shallower, colder parts of the crust, it crystallizes (freezes) dynamically:
* The crystallization rate is temperature-dependent:
  $$\lambda_{eff} = \lambda_{freeze} \cdot e^{-\lambda_{freeze\_tdep} \cdot (T - T_{top})}$$
  where `lambda_freeze = 1e-13` s$^{-1}$ and `lambda_freeze_tdep = 0.002` K$^{-1}$.
* The magma fraction is reduced by $\Delta f_{magma} = \min(f_{magma}, f_{magma} \cdot dt \cdot \lambda_{eff})$.
* Freezing releases latent heat back into the node coordinates, raising their temperatures by:
  $$\Delta T = \frac{\Delta f_{magma} \cdot L_{fusion}}{4 \cdot C_{p\_eff}}$$
* Elements within the diking conduit that accumulate frozen magma are dynamically converted into arc crust (Phase 14).

### 4. Feedbacks on Rheology
The presence of magma significantly weakens both the viscous and plastic strength of the host rock:
* **Viscous Weakening**:
  $$\eta_{eff} = \eta_0 \cdot e^{-\text{weaken\_ratio\_viscous} \cdot f_{magma} / f_{magma\_max}}$$
* **Plastic Weakening**: Cohesion and friction angle are weakened linearly:
  $$c = c_0 \cdot \left(1 - (1 - \text{weaken\_ratio\_plastic}) \cdot \frac{f_{magma}}{f_{magma\_max}}\right)$$
  $$\phi = \phi_0 \cdot \left(1 - (1 - \text{weaken\_ratio\_plastic}) \cdot \frac{f_{magma}}{f_{magma\_max}}\right)$$
  Since `weaken_ratio_plastic = 1.0` and `weaken_ratio_viscous = 1.0`, the presence of magma reduces the viscosity to $1/e \approx 37\%$ of its original value, and cohesion/friction angle to 0!

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
This script reads the binary outputs and generates three publication-ready figures inside the `images/` directory:
1. **`images/subduction_full_zone.png`**: A full 2D cross-section showing the subducting slab, phase distribution, temperature isotherms, and plate velocity field.
2. **`images/subduction_mantle_wedge.png`**: A two-panel zoomed-in view of the mantle wedge ($X \in [300, 700]\text{ km}$, $Z \in [-150, 0]\text{ km}$) displaying:
   * **Top Panel**: Active phase boundaries (serpentinite decoupling, hydrated mantle wedge) and velocity vectors.
   * **Bottom Panel**: Magma/melt generation zone in the mantle wedge beneath the volcanic arc.
3. **`images/subduction_evolution.png`**: A three-panel evolutionary sequence showing how the slab sinks over time.

---

## 8. Analyzing the Results

### Slab Eclogitization and Angle
In the full profile (`images/subduction_full_zone.png`), watch how the subducting slab (green olivine mantle + blue oceanic crust) sinks. As the blue oceanic crust descends below $\sim 60\text{ km}$, it transforms into dark red **Eclogite** (Phase 13). The high density of eclogite drives the slab down, steepening the subduction angle.

### Mantle Wedge Flow and Hydration
Observe the velocity vectors in the mantle wedge. The descending slab drags the adjacent wedge mantle downwards, establishing a **corner flow`** convection cell. 
Directly above the slab, you will see a thin yellow layer of **Serpentinite** (Phase 9) at shallow depths, which transitions into **Hydrated Mantle** (Phase 16) at greater depths.

### Arc Magmatism
Look at the zoomed-in plot (`images/subduction_mantle_wedge.png`). The hydrated mantle (Phase 16) sits in the hot core of the mantle wedge ($T > 800^\circ\text{C}$). Because of hydration, this region starts melting, visible as a red-colored anomaly of high **magma fraction**. The extracted magma travels vertically through the diking channel (Phase 14) directly above the melting zone to form the volcanic arc at the surface.
