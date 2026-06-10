# GeoFLAC Tutorial: Convergent Wedge with Basal Weak Detachment

This tutorial explains the physical principles, model configuration, custom boundary conditions, and visualization of an **Orogenic / Accretionary Convergent Wedge** simulation in **GeoFLAC**.

---

## 1. Physical Principles (Critical Taper Theory)

A **convergent wedge** (also known as an accretionary prism or fold-and-thrust belt) forms when a layered rock sequence is compressed horizontally over a weak basal detachment layer against a rigid backstop. This simulates tectonic settings like the Himalayas, the Swiss Alps, or accretionary prisms at subduction zones.

The structural geometry and evolution of convergent wedges are governed by **Critical Taper Theory** (Chapple, 1978; Davis, Suppe, & Dahlen, 1983). A wedge is considered "critical" when it is at the verge of sliding along its basal detachment and yielding internally everywhere due to shear failure.

```
                  Wedge Topographic Slope (alpha)
                           /
                          v             +---- Rigid Backstop
                       .·`|             |
                   .·`    |             v
  Vx (<-)      .·`        |             |   /|
    \      .·`            |             |  / |
     v .·` Accretionary   |             | /  |
  ====+  Wedge (Phase 1)  |             |/   |
  ===========================================+  <-- Basal Detachment (Phase 2)
  -------------------------------------------+  <-- Base (fixed vertical, free-slip)
      \_________________________________/
             Basal Slope (beta)
```

The critical taper angle ($\alpha + \beta$, where $\alpha$ is the topographic slope and $\beta$ is the basal dip) is a balance between the tectonic push and gravitational collapse, mathematically described as:

$$\alpha + \beta \approx \frac{(1 - \lambda)\mu_b + \beta}{(1 - \lambda)\mu + 1}$$

Where:
* $\mu = \tan \phi$ is the internal friction coefficient of the wedge ($\phi$ is internal friction angle).
* $\mu_b = \tan \phi_b$ is the basal detachment friction coefficient ($\phi_b$ is basal friction angle).
* $\lambda$ is the pore-fluid pressure ratio.

In this tutorial:
- **Accretionary Wedge (Phase 1)** is stiff with high friction ($\phi = 30^\circ$) and high cohesion ($c = 20 \text{ MPa}$).
- **Basal Detachment (Phase 2)** is extremely weak with low friction ($\phi_b = 5^\circ$) and low cohesion ($c = 1 \text{ MPa}$).

Because the basal detachment is weak, compression from the left allows shear to localize immediately along the base, letting the upper layer slide. When the sliding layer hits the rigid backstop (right boundary), the rock accumulates stress, yields internally, and deforms upwards, building a beautiful accretionary wedge with active conjugate **pro-thrusts** (dipping toward the push) and **retro-thrusts** (dipping toward the backstop).

---

## 2. Model Setup

The model represents a vertical 2D cross-section of a 10 km thick crustal column undergoing horizontal convergence.

### Geometry and Mesh
* **Dimensions**: $100 \text{ km}$ wide ($X \in [0, 100] \text{ km}$), $10 \text{ km}$ deep ($Z \in [-10, 0] \text{ km}$).
* **Grid Resolution**: $50 \times 20$ elements, producing a regular grid of square-like elements ($2 \text{ km} \times 0.5 \text{ km}$).
* **Mesh Coordinates**: Defined in `convergent_wedge.inp`:
  ```fortran
  50,20            number of _elements_ in X and Z directions: (nx-1),(nz-1)
  0.e+3,0.           x0,z0 begin.coord
  100.e+3,-10.e+3    rxbo,rzbo (100 km wide, 10 km deep)
  ```

### Two-Layer Phase Distribution
We partition the domain into a stiff upper crust and a weak basal detachment:
* **Upper layer (Phase 1, 0–9 km depth)**: Stiff elastoplastic rock ($c = 20 \text{ MPa}$, $\phi = 30^\circ$).
* **Basal layer (Phase 2, 9–10 km depth)**: A $1 \text{ km}$ thick weak basal detachment layer ($c = 1 \text{ MPa}$, $\phi = 5^\circ$).
* **Parsing**: Configured in `convergent_wedge.inp` using `nzone_age` layer boundaries:
  ```fortran
  1              - nzone_age
   1, 100., 0, 0, 1, 51
      2, 9.0
      1, 2
  ```
  This places the interface at $9 \text{ km}$ depth, defining Phase 1 in the upper 18 element rows and Phase 2 in the bottom 2 element rows.

### Gravity and Thermal Settings
* **Gravity**: Set to $10.0 \text{ m/s}^2$ to balance tectonic push against gravitational loading.
* **Thermal Geotherm**: Uniform temperature $T = 0^\circ\text{C}$ everywhere (`t_top = 0.0`, `t_bot = 0.0`). This suppresses any temperature-dependent viscous creep, allowing the model to act as a purely mechanical, friction-controlled accretionary sandbox.

---

## 3. Boundary Conditions

The mechanical boundary conditions compress the column horizontally and allow easy basal sliding:

### Left Boundary (Side 1) - Piece-wise Linear Shear Inflow
We apply compression from the left boundary. To avoid numerical stress concentrations at the bottom corner and model a realistic shear zone in the basal detachment:
* **Top 9 km (nodes 1 to 19)**: Constant horizontal velocity $V_x = 1.0 \text{ cm/yr} \approx 3.1688 \times 10^{-10} \text{ m/s}$.
* **Bottom 1 km weak layer (nodes 19 to 21)**: Horizontal velocity decreases linearly from $V_x = 1.0 \text{ cm/yr}$ at the top of the layer to $V_x = 0$ at the bottom boundary ($V_x = z \times 1 \text{ cm/km/yr}$, where $z$ is distance from bottom).

In GeoFLAC, this is natively achieved by splitting Side 1 into two boundary segments:
1. **BC 1 (nodes 1 to 19)**: Constant function ($a = 3.1688\times 10^{-10}$):
   `1         1    19    10    3.1688e-10  0.   0.    0.    0.    0.     0.    0.    0.`
2. **BC 2 (nodes 19 to 21)**: Linear function ($a = 3.1688\times 10^{-10}$, $b = -3.1688\times 10^{-10}$):
   `1         19   21    10    3.1688e-10 -3.1688e-10 0. 0.   0.    0.     0.    0.    0.`
   *Here, $x$ is non-dimensionalized along the segment from 0 (at node 19) to 1 (at node 21), yielding a perfect linear gradient:*
   $$V_x(x) = 3.1688 \times 10^{-10} \times (1 - x) \text{ m/s}$$

### Other Boundaries
* **Right Boundary (Side 3, nodes 1 to 21)**: Horizontal velocity $V_x = 0.0$ (`nbc = 10`), acting as a rigid vertical backstop.
* **Bottom Boundary (Side 2, nodes 1 to 51)**: Vertical velocity $V_z = 0.0$ (`nbc = 01`). Horizontal velocity is free-slip, permitting the basal detachment to slide.

---

## 4. Mohr-Coulomb Rheology with Strain Softening

We use purely elastoplastic Mohr-Coulomb rheology (`irheol = 6`). Under stress, elements behave elastically until satisfying the Mohr-Coulomb yield criterion, after which they yield plastically.

* **Wedge (Phase 1) softening**:
  - Peak Cohesion $c_1 = 20 \text{ MPa}$, Residual Cohesion $c_2 = 2 \text{ MPa}$ (`coh1 = 2.0e7`, `coh2 = 2.0e6`).
  - Peak Friction Angle $\phi_1 = 30^\circ$, Residual Friction Angle $\phi_2 = 20^\circ$ (`fric1 = 30.0`, `fric2 = 20.0`).
  - Weakening occurs over plastic strain $\epsilon_p \in [0.0, 0.1]$ (`pls1 = 0.0`, `pls2 = 0.1`).
* **Basal Detachment (Phase 2) softening**:
  - Peak Cohesion $c_1 = 1 \text{ MPa}$, Residual Cohesion $c_2 = 0.1 \text{ MPa}$ (`coh1 = 1.0e6`, `coh2 = 1.0e5`).
  - Constant Friction Angle $\phi = 5^\circ$ (`fric1 = 5.0`, `fric2 = 5.0`).
  - Weakening occurs over plastic strain $\epsilon_p \in [0.0, 0.1]$ (`pls1 = 0.0`, `pls2 = 0.1`).

---

## 5. Running the Simulation and Plotting

### Step 1: Run the GeoFLAC Solver
Clean any old outputs and run the solver in this directory:
```bash
rm -f *.0 *.rs *.vts _contents.* _markers.* pisos.rs time.rs vbc.s output.asc sys.msg
../../src/flac convergent_wedge.inp
```
The solver will run to **1.5 Myr** (`1500.0` Kyr) of total convergence, outputting data frames every **0.15 Myr** (`150.0` Kyr). Because it is a purely mechanical elastoplastic grid of $50 \times 20$ elements, it executes extremely fast.

### Step 2: Plot the Accretionary Wedge
Run the provided Python plotting script:
```bash
python3 plot_convergent_wedge.py
```
This script reads the binary outputs and generates two premium visualizations:
1. **`convergent_wedge.png`**: A detailed two-panel plot of the final deformed state showing:
   - **Top Panel**: Accumulated Plastic Strain (`aps`) highlighting localized thrust faults (conjugate shear bands), overlaid with velocity vectors showing wedge uplift.
   - **Bottom Panel**: Layer phases (`phase`) showing accretionary topography and basal slip, overlaid with a deformed grid mesh.
2. **`plots/evolution_convergent_wedge.png`**: A three-panel evolutionary sequence showing how the thrust sheets nucleate and step outwards over time.

---

## 6. Analysis of Results

### Fault Localization and Thrust Imbricates
In the upper panel of `convergent_wedge.png`, you will observe that strain localizes into distinct narrow bands of high plastic strain (`aps`). These represent **thrust faults**!
Because the wedge material is compressed against the rigid backstop:
- Stress builds up near the backstop and exceeds the yield strength.
- Conjugate thrust faults (slanted shear bands) form: **pro-thrusts** (dipping to the left) and **retro-thrusts** (dipping to the right).
- As convergence continues, old thrust faults lock up (due to isostatic gravity load) and new thrust faults step outwards (propagate to the left) into the undeformed crust, forming a classic **imbricate thrust fan**.

### Topographic Wedge and Basal Slip
In the lower panel, the beige upper crust (Phase 1) thickens and deforms upwards, building a beautiful wedge-shaped topography (accretionary prism) sloping away from the rigid backstop.
The blue-grey basal weak layer (Phase 2) remains relatively flat at the bottom, absorbing most of the shear strain and acting as a perfect detachment layer (décollement) that lets the wedge slide and stack dynamically.
The deformed grid lines show strong vertical thickening and horizontal shortening within the wedge, while showing simple shear within the basal detachment layer.
