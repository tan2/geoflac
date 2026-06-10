# GeoFLAC Tutorial: Simple Elastic Bar Compression

This tutorial explains the setup, boundary conditions, physical results, and analytical verification of the simple elastic bar compression benchmark in **GeoFLAC**.

---

## 1. Running the Simulation and Plotting

### Step 1: Run the GeoFLAC Solver
Execute the compiled `flac` binary using the input configuration file:
```bash
../../src/flac bar.inp
```
The solver will run for 3,000 steps, simulating a total time of $0.0108$ Myr ($10,800$ years) and creating binary files `ezz.0`, `szz.0`, and `pres.0`.

### Step 2: Plot the Stress-Strain Curve
Run the provided Python plotting script:
```bash
python3 plot_elastic.py
```
This script reads the binary files, averages the strain and total vertical stress across the domain for each output frame, plots them against the analytical line ($E_{eff} = 80$ GPa), and saves the figure to `stress_strain.png`.

---

## 2. Model Setup

The model represents a vertical two-dimensional column (bar) of homogenous elastic rock undergoing vertical compression. 

### Geometry and Mesh
* **Dimensions**: 1000 meters wide ($X \in [0, 1000]$ m), 3000 meters high ($Z \in [-3000, 0]$ m).
* **Grid Resolution**: $3 \times 10$ elements in the $X$ and $Z$ directions, yielding a grid of elements ($333.3 \times 300$ m).

![GeoFLAC grid structure and node/element relationship](images/mesh.png)

#### Grid Elements and Nodes:
GeoFLAC uses a structured 2D quadrilateral mesh. It distinguishes between **nodes** and **elements**:
* **Nodes**: Located at element corners. Nodal variables include coordinates (`x`, `z`), velocities (`vx`, `vz`), and temperature (`temperature`). For a grid of $3 \times 10$ elements, there are $4 \times 11 = 44$ nodes. Node indices range from `(1, 1)` (top-left) to `(nx, nz)` (bottom-right).
* **Elements**: The quadrilaterals bounded by nodes. Elemental variables include stresses (`sxx`, `szz`, `sxz`), pressure (`pres`), phases (`phase`), and accumulated plastic strain (`aps`). For a $3 \times 10$ element grid, there are $3 \times 10 = 30$ elements. Element indices range from `(1, 1)` to `(nex, nez)`.

### Material Properties
The material is a linear isotropic elastic rock defined in the input file `bar.inp` with the following parameters:
* **Rheology Type (`irheol`)**: Set to `1` (Elastic) for linear isotropic elastic behavior. Refer to the [Rheology Types table](../../doc/input_description.md#phases--rheology) for other options.
* **Lamé constant ($\lambda$)**: $3.0 \times 10^{10} \text{ Pa}$ (`Lame:rl`)
* **Shear Modulus ($\mu$)**: $3.0 \times 10^{10} \text{ Pa}$ (`Lame:rm`)
* **Density ($\rho$)**: $2700 \text{ kg/m}^3$ (`den`)
* **Gravity ($g$)**: $0.0 \text{ m/s}^2$ (pure compression without gravity-induced pre-stress or hydrostatic gradients).

---

## 3. Boundary Conditions

The mechanical boundary conditions are configured in [`bar.inp`](bar.inp) to produce a **uniaxial stress** state along the vertical axis:

```fortran
;nofside  nbc1 nbc2  nbc   a       b    c     d     e     f      g     h     i 
2         1    4     01    0.0     0.   0.    0.    0.    0.     0.    0.    0.  ; Bottom fixed in Z
4         1    4     01   -1.e-9   0.   0.    0.    0.    0.     0.    0.    0.  ; Top compression in Z
```

### Syntax and Parameters Breakdown:
* **Side Selection (`nofside`)**:
  * **`2`**: Represents the **bottom boundary** ($Z = -3000\text{ m}$).
  * **`4`**: Represents the **top boundary** ($Z = 0\text{ m}$).
* **Node Range (`nbc1` to `nbc2`)**:
  * The mesh is set to $3 \times 10$ elements, giving 4 nodes along the horizontal axis. Specifying nodes **`1` to `4`** applies the boundary condition continuously across the entire width of the domain.
* **Boundary Condition Type (`nbc`)**:
  * **`01`**: Specifies a **vertical velocity ($V_z$) constraint** in meters/second. Refer to the [Boundary Condition Types table](../../doc/input_description.md#mechanical-conditions) for other options.
* **Spatial Profile Coefficient (`a`)**:
  * The velocity profile is governed by a spatial function, where `a` is the constant coefficient:
    * **Bottom (`nofside = 2`)**: Setting `a = 0.0` m/s restricts vertical movement, fixing the base vertically ($V_z = 0$).
    * **Top (`nofside = 4`)**: Setting `a = -1.e-9` m/s moves the boundary vertically downward ($V_z = -1.0 \times 10^{-9}\text{ m/s}$), driving the compression.
* **Horizontal & Lateral Conditions (Implicit)**:
  * Since no horizontal velocity constraints (`nbc = 10`) are defined on the top or bottom, they default to **free-slip** (free to expand sideways).
  * Since no boundary conditions are defined for the left (`nofside = 1`) and right (`nofside = 3`) boundaries, they default to completely **free/traction-free surfaces** ($\sigma_{xx} = 0$, $\sigma_{xz} = 0$).

---

## 4. Analytical Formulation (2D Plane Strain)

### What is Plane Strain?
**Plane Strain** is a simplified state of strain in solid mechanics where the deformation in one direction (usually the out-of-plane direction, $Y$, or axis 3 in GeoFLAC) is constrained to be zero. Specifically:
$$\epsilon_{yy} = 0, \quad \gamma_{xy} = 0, \quad \gamma_{zy} = 0$$

This assumption is physically applicable to structures that are very long in one dimension relative to their cross-section (such as tunnels, retaining walls, or regional tectonic columns). 

Because the out-of-plane expansion or contraction is prevented, a normal stress is induced in that direction:
$$\sigma_{yy} = \nu (\sigma_{xx} + \sigma_{zz})$$
This out-of-plane stress acts as a confinement, which increases the effective stiffness of the rock column, leading to the **effective vertical modulus** $E_{eff} = E / (1 - \nu^2)$ instead of the standard Young's Modulus $E$.

### Mathematical Derivation
Using the 3D isotropic Hooke's Law under these boundary conditions:
1. **Zero Out-of-Plane Strain ($\epsilon_{yy} = 0$)**:
   $$\epsilon_{yy} = \frac{1}{E} \left[ \sigma_{yy} - \nu (\sigma_{xx} + \sigma_{zz}) \right] = 0 \implies \sigma_{yy} = \nu (\sigma_{xx} + \sigma_{zz})$$
2. **Zero Horizontal Confinement ($\sigma_{xx} = 0$)**:
   $$\sigma_{yy} = \nu \sigma_{zz}$$
3. **Vertical Strain ($\epsilon_{zz}$)**:
   $$\epsilon_{zz} = \frac{1}{E} \left[ \sigma_{zz} - \nu (\sigma_{xx} + \sigma_{yy}) \right] = \frac{1 - \nu^2}{E} \sigma_{zz}$$

Therefore, the vertical stress-strain relationship is governed by the **effective vertical modulus** $E_{eff}$ :
$$\sigma_{zz} = E_{eff} \epsilon_{zz}$$
where:
$$E_{eff} = \frac{E}{1 - \nu^2} = \frac{4\mu(\lambda + \mu)}{\lambda + 2\mu}$$

### Physical Values:
* **Poisson's Ratio ($\nu$)**:
  $$\nu = \frac{\lambda}{2(\lambda + \mu)} = \frac{3 \times 10^{10}}{2(3 \times 10^{10} + 3 \times 10^{10})} = 0.25$$
* **Young's Modulus ($E$)**:
  $$E = 2\mu(1 + \nu) = 2 \times (3 \times 10^{10}) \times 1.25 = 75 \text{ GPa} = 75,000 \text{ MPa}$$
* **Effective Modulus ($E_{eff}$)**:
  $$E_{eff} = \frac{75 \text{ GPa}}{1 - 0.25^2} = 80 \text{ GPa} = 80,000 \text{ MPa}$$

Under pure compression, the Poisson expansion in the horizontal direction is:
$$\epsilon_{xx} = -\frac{\lambda}{\lambda + 2\mu} \epsilon_{zz} = -\frac{1}{3} \epsilon_{zz}$$

### Configuration of Out-of-Plane Stress in GeoFLAC
While the elastic equations in GeoFLAC natively assume plane strain ($\epsilon_{yy} = 0$), calculating and tracking the out-of-plane stress component ($\sigma_{yy}$) for complex materials (like Maxwell creep and plastic yielding) is optional. 

This behavior is controlled by the **`ndim`** parameter in the [`PROCESS CONTROL`](../../doc/input_description.md#process-control) section of the input file:
* **`ndim = 2` (Default)**: Normal 2D mode. Only the in-plane stress components ($\sigma_{xx}$, $\sigma_{zz}$, $\sigma_{xz}$) are used to check for plastic yielding or to compute Maxwell relaxation. The deviatoric out-of-plane stress ($\sigma'_{yy}$) is assumed to be zero to simplify and speed up calculations.
* **`ndim = 3`**: In/out-of-plane stress mode (3D mode). The out-of-plane stress ($\sigma_{yy}$) is computed dynamically and included in plastic yield criteria (e.g., Mohr-Coulomb in 3D principal stress space) and viscoelastic creep flow.

---

## 5. Understanding GeoFLAC Stress Output

> [!IMPORTANT]
> **Stress Representation in GeoFLAC binary files:**
> * `szz.0` stores the **deviatoric vertical stress** ($\sigma'_{zz} = \sigma_{zz} - P$), not the total stress.
> * `pres.0` stores the **mean stress / pressure** ($P$).
> * Both stress and pressure variables are output in **Kilobars** ($1 \text{ kb} = 100 \text{ MPa} = 10^8 \text{ Pa}$).

To reconstruct the total physical vertical stress $\sigma_{zz}$, you must sum the deviatoric and mean stress components:
$$\sigma_{zz} = \sigma'_{zz} + P$$
Multiplying by $100.0$ converts the result from Kilobars to standard MegaPascals (MPa).

---

## 6. Simulation Results

### Stress and Strain Evolution
The simulation outputs 8 frames, showing excellent linear scaling matching the analytical slope.

### Verification Chart
The generated plot is saved to `images/stress_strain.png`:

![Stress-Strain Verification Chart](images/stress_strain.png)

1. **Blue Line**: The averaged simulation stress vs strain.
2. **Red Dashed Line**: The theoretical analytical solution slope of $80$ GPa.

The two lines align with exceptional accuracy, confirming the physical validity and high numerical precision of **GeoFLAC**'s elastic mechanics solver and boundary condition implementation.
