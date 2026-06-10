# GeoFLAC Tutorial: Metamorphic Core Complex (MCC) Extension

This tutorial guides you through the setup, running, physical principles, and visualization of a **Metamorphic Core Complex (MCC)** crustal extension simulation using **GeoFLAC**.

---

## 1. Running the Simulation and Plotting

### Step 1: Run the GeoFLAC Solver
Clean old outputs and execute the compiled `flac` binary in the directory:
```bash
rm -f *.0 *.rs *.vts _contents.* _markers.* pisos.rs time.rs vbc.s output.asc sys.msg
../../src/flac core_complex.inp
```
The solver will run for approximately 38,500 steps, simulating a total time of **2.0 Myr** and writing output frames every **0.2 Myr**.

### Step 2: Generate the Visualization
Run the provided Python plotting script:
```bash
python3 plot_core_complex.py
```
This script reads the output files, extracts the grid coordinates, accumulated plastic strain (`aps`), lithology phases (`phase`), velocity field (`vx`, `vz`), and temperature (`temp`), and creates two premium visualizations:
1. **`images/core_complex.png`**: A detailed, publication-ready two-panel plot of the final state (1.8 Myr) showing detachment faulting, exhumation velocity vectors, and temperature isotherms advection.
2. **`images/evolution_core_complex.png`**: An evolutionary chart showing how the shear bands and detachment faults develop over time.

### Step 3: Convert Output to VTK Format (VTS)
To visualize the spatial distribution of stresses, strain rates, temperature, and material phases in ParaView or VisIt, convert the binary output files to `.vts` structured grid files using the provided utility:
```bash
python3 ../../util/flac2vtk.py .
```
This will generate `flac.000001.vts` to `flac.000010.vts` in your current directory.

---

## 2. Physical Background

A **Metamorphic Core Complex (MCC)** is a spectacular tectonic feature formed during high-magnitude crustal extension. It is characterized by the exhumation of middle-to-deep crustal rocks (the "metamorphic core") through a low-angle detachment fault system, exposing dome-like structures at the surface.

The formation of an MCC requires specific thermomechanical conditions:
1. **A strong, brittle upper crust** that can localize strain into major detachment faults.
2. **A warm, highly ductile lower crust** that can flow laterally and vertically to fill the space opened by upper crustal extension, preventing the crust from breaking completely apart.
3. **Visco-Elasto-Plastic (VEP) Rheology**, which allows materials to dynamically transition from elastic deformation to brittle faulting (plastic yielding) or ductile flow (viscous relaxation) based on temperature, pressure, and strain rate.

```
       Extensional Velocity (<- Vx)                      Extensional Velocity (Vx ->)
            \                                                               /
             v                                                             v
             +----------------------- Detachment Fault --------------------+
             |   Brittle Upper Crust (Mohr-Coulomb Yielding & Softening)   |
             +~~~~~~~~~~~~~~~~~~~~ Brittle-Ductile Transition ~~~~~~~~~~~~~+
             |                                                             |
             |           Ductile Lower Crust (Viscous Creep Flow)          |
             |                       ^                                     |
             +---------------------- | ------------------------------------+
                                     |
                             Asthenospheric Uplift
```

---

## 3. Model Setup

The model represents a vertical 2D cross-section of a continental lithospheric crust undergoing horizontal extension.

### Geometry and Mesh
* **Dimensions**: $100 \text{ km}$ wide ($X \in [0, 100] \text{ km}$), $30 \text{ km}$ deep ($Z \in [-30, 0] \text{ km}$).
* **Grid Resolution**: $50 \times 15$ elements, creating a regular grid where each element is $2 \text{ km} \times 2 \text{ km}$.
* **Mesh Coordinates**: Defined in `core_complex.inp`:
  ```fortran
  50,15            number of _elements_ in X and Z directions
  0.e+3,0.           x0,z0 begin.coord
  100.e+3,-30.e+3    rxbo,rzbo (100 km wide, 30 km deep)
  ```

### Geotherm and Initial Thermal State
We apply a linear geothermal gradient (`ictherm = 21`) starting from the surface ($T_{\text{top}} = 0^\circ\text{C}$) to the bottom of the crust ($T_{\text{bot}} = 800^\circ\text{C}$). This corresponds to a hot geothermal gradient of $\sim 26.7^\circ\text{C/km}$, typical of active rift systems.
* **Surface Temperature**: $0^\circ\text{C}$ (`t_top = 0.0`)
* **Bottom Temperature**: $800^\circ\text{C}$ (`t_bot = 800.0` and bottom heat boundary `1 800.`)

### Crustal Phases and Rheology
We partition the continental crust into two distinct lithological layers:
1. **Upper Crust (Phase 1, 0–12 km depth)**: Wet Quartz VEP with strain softening.
   * **Rheology Type (`irheol`)**: Set to `12` (Visco-Elasto-Plastic) in the input file.
   * **Strain Softening**: As plastic strain accumulates, the material weakens.
     * Cohesion decreases from $c_1 = 20 \text{ MPa}$ to $c_2 = 2 \text{ MPa}$ (`coh1 = 2.0e7`, `coh2 = 2.0e6`).
     * Friction angle decreases from $\phi_1 = 30^\circ$ to $\phi_2 = 15^\circ$ (`fric1 = 30.0`, `fric2 = 15.0`).
     * Softening interval occurs over plastic strain $\epsilon_p \in [0.0, 0.1]$ (`pls1 = 0.0`, `pls2 = 0.1`).
2. **Lower Crust (Phase 2, 12–30 km depth)**: Wet Quartz VEP with constant material parameters.
   * **Rheology Type (`irheol`)**: Set to `12` (Visco-Elasto-Plastic) in the input file.
   * Cohesion and friction angle remain high ($c = 20 \text{ MPa}$, $\phi = 30^\circ$) to prevent shallow brittle failure, but because of high temperatures ($> 320^\circ\text{C}$), its dislocation creep viscosity is very low, making it behave like a highly ductile fluid.

---

## 4. Boundary Conditions

The model is extended horizontally by applying constant velocities on the left and right boundaries:

1. **Left Boundary (Side 1)**: Pulling leftward at $V_x = -0.5 \text{ cm/yr} = -1.58438 \times 10^{-10} \text{ m/s}$ across all 16 nodes.
2. **Right Boundary (Side 3)**: Pulling rightward at $V_x = +0.5 \text{ cm/yr} = +1.58438 \times 10^{-10} \text{ m/s}$ across all 16 nodes.
3. **Bottom Boundary (Side 2)**: Constrained to allow asthenospheric pressure support (`nyhydro = 2`) and mantle return flow (`iphsub = 2`). This allows the bottom of the crust to flex and dome upward under isostatic balance.
4. **Top Boundary (Side 4)**: Free surface with topography diffusion (`topo_removal_rate = 1.0e-6`), which simulates surface erosion and sedimentation in the rift basins.

---

## 5. Inhomogeneities (The Weak Seed)

To localize strain and nucleate a major detachment fault in the center of the domain, we insert a small **weak seed** at the surface:
* **Location**: $X \in [48, 52] \text{ km}$ (indices $25 \to 27$), $Z \in [-2, 0] \text{ km}$ (indices $1 \to 2$).
* **Mechanism**: Set the initial accumulated plastic strain to `1.0` (`init.pl.strain = 1.0`). This immediately softens the cohesion and friction angle of this zone to their residual values ($c = 2 \text{ MPa}$, $\phi = 15^\circ$), forcing extensional strain to localize there.

---

## 6. VEP Rheological Formulation

In Visco-Elasto-Plastic (VEP) rheology (`irheol = 12`), the solver computes the local stresses for both **viscoelastic** Maxwell relaxation and **elasto-plastic** Mohr-Coulomb yielding, and applies the one that results in a lower stress state:

$$\sigma = \min\left( \sigma_{\text{viscoelastic}}, \sigma_{\text{plastic}} \right)$$

### 1. Viscoelastic (Maxwell Creep)
Ductile flow is modeled using power-law dislocation creep for Wet Quartz:

$$\dot{\epsilon}_v = A \cdot \sigma^n \cdot \exp\left( -\frac{E + P V}{R T} \right)$$

Where:
* $A = 1.25 \times 10^{-1}$ $\text{Pa}^{-n}\text{s}^{-1}$ (Pre-exponential factor)
* $n = 3.0$ (Stress exponent)
* $E = 2.76 \times 10^5$ $\text{J/mol}$ (Activation energy)
* $R = 8.314$ $\text{J/(mol}\cdot\text{K)}$ (Gas constant)
* $T$ is temperature in Kelvin.

### 2. Elasto-Plastic (Mohr-Coulomb Yielding)
Brittle failure occurs when stresses satisfy the Mohr-Coulomb yield criterion:

$$f(\sigma_1, \sigma_3) = \sigma_3 - \sigma_1 N_\phi + 2 c \sqrt{N_\phi} = 0$$

Where $N_\phi = \frac{1 + \sin \phi}{1 - \sin \phi}$, $c$ is cohesion, and $\sigma_1 \ge \sigma_3$ are the principal stresses (tension positive).

---

## 7. Analysis of Results

### Detachment Faulting and Shear Localization
In the upper panel of `images/core_complex.png`, you will observe that strain localizes into a low-angle shear band originating from the central weak seed and cutting through the brittle upper crust. This is the **detachment fault**. 
Because of the horizontal extension, the upper crust is pulled apart, forming an asymmetrical rift valley (graben) bounded by this active detachment fault.

![Metamorphic Core Complex final State](images/core_complex.png)

### Ductile Exhumation and Domelike Uplift
In the lower panel, you can see that the boundary between Phase 1 (Upper Crust) and Phase 2 (Lower Crust) domes upward significantly in the center. 
The highly ductile lower crust flows horizontally from the sides and vertically upward into the rift core, filling the space vacated by the extending upper crust. The velocity vectors (arrows) show a strong upward exhumation flow beneath the detachment fault.

### Thermal Advection (Isotherms Pull-up)
The red dashed lines represent the temperature isotherms ($200^\circ\text{C}$ to $700^\circ\text{C}$). You will notice that the isotherms are significantly pulled upward in the center of the domain. This represents **advective heat transfer**: the rapid upward flow of hot lower crust physically carries heat toward the surface, steepening the geothermal gradient locally and helping to keep the exhuming core hot and ductile!

The evolutionary progress of the core complex is shown below:

![Metamorphic Core Complex Strain Evolution](images/evolution_core_complex.png)

---

> [!TIP]
> **Key takeaways for core complex physics:**
> * MCCs cannot form if the lower crust is cold (e.g. $T_{\text{bot}} < 500^\circ\text{C}$). A cold lower crust would be brittle and would break along narrow faults rather than flowing into a dome.
> * Strain softening in the upper crust is critical for localizing shear into a single major detachment fault rather than numerous distributed minor faults.
