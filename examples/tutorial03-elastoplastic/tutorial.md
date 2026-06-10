# GeoFLAC Tutorial: Mohr-Coulomb Elastoplastic Compression

This tutorial explains the setup, boundary conditions, physical results, and analytical verification of the Mohr-Coulomb elastoplastic bar compression benchmark in **GeoFLAC**.

---

## 1. Running the Simulation and Plotting

### Step 1: Run the GeoFLAC Solver
Execute the compiled `flac` binary in the directory:
```bash
rm -f *.0 *.rs *.vts _contents.* _markers.* pisos.rs time.rs vbc.s output.asc sys.msg
../../src/flac plastic.inp
```
The solver will run for 500 steps, simulating a total time of $0.5$ Kyr and writing output files `ezz.0`, `szz.0`, and `pres.0` every $0.05$ Kyr.

### Step 2: Plot the Stress-Strain Curve
Run the provided Python plotting script:
```bash
python3 plot_elastoplastic.py
```
This script reads the binary files, averages the vertical strain and total vertical stress across the domain for each output frame, plots them against the analytical Mohr-Coulomb yield path ($\sigma_C = 69.28$ MPa), and saves the figure to `images/stress_strain_yield.png`.

### Step 3: Convert Output to VTK Format (VTS)
To visualize the spatial distribution of stresses, strains, and yielding in ParaView or VisIt, convert the binary output files to `.vts` structured grid files using the provided utility:
```bash
python3 ../../util/flac2vtk.py .
```
This will generate `flac.000001.vts` to `flac.000010.vts` in your current directory.

---

## 2. Model Setup

The model represents a vertical two-dimensional column (bar) of homogenous elastoplastic rock undergoing vertical compression. Under compression, the material behaves elastically until its state of stress satisfies the **Mohr-Coulomb yield criterion**, after which it deforms plastically at a constant yield stress.

### Geometry and Mesh
* **Dimensions**: 1000 meters wide ($X \in [0, 1000]$ m), 3000 meters high ($Z \in [-3000, 0]$ m).
* **Grid Resolution**: $10 \times 30$ elements in the $X$ and $Z$ directions, yielding a regular grid of square elements ($100 \times 100$ m).

### Material Properties
The material is a Mohr-Coulomb elastoplastic rock with strain softening defined in the input file `plastic.inp` with the following parameters:
* **Rheology Type (`irheol`)**: Set to `6` (Elasto-plastic, Mohr-Coulomb) in the input file. Refer to the [Rheology Types table](../../doc/input_description.md#phases--rheology) for other options.
* **Lamé constant ($\lambda$)**: $3.0 \times 10^{10} \text{ Pa}$ (`Lame:rl`)
* **Shear Modulus ($\mu$)**: $3.0 \times 10^{10} \text{ Pa}$ (`Lame:rm`)
* **Peak Cohesion ($c_1$)**: $2.0 \times 10^7 \text{ Pa} = 20 \text{ MPa}$ ([`coh1`](../../doc/input_description.md#phases--rheology))
* **Residual Cohesion ($c_2$)**: $2.0 \times 10^6 \text{ Pa} = 2 \text{ MPa}$ ([`coh2`](../../doc/input_description.md#phases--rheology) - defining a 90% cohesion weakening)
* **Plastic Strain weakening limits**: $[0.0, 0.1]$ ([`pls1`, `pls2`](../../doc/input_description.md#phases--rheology) - softening occurs over 10% plastic strain)
* **Friction Angle ($\phi$)**: $30.0^\circ$ ([`fric1`, `fric2`](../../doc/input_description.md#phases--rheology) - constant friction angle)
* **Dilatancy Angle ($\psi$)**: $0.0^\circ$ ([`dilat1`, `dilat2`](../../doc/input_description.md#phases--rheology) - non-associated plastic flow, no plastic volume change)
* **Density ($\rho$)**: $2700 \text{ kg/m}^3$ (`den`)
* **Gravity ($g$)**: $0.0 \text{ m/s}^2$ (pure compression without gravity-induced pre-stress or hydrostatic gradients).

---

## 3. Boundary Conditions

The mechanical boundary conditions are configured in [`plastic.inp`](plastic.inp) to compress the column vertically:

```fortran
;nofside  nbc1 nbc2  nbc   a       b    c     d     e     f      g     h     i 
2         1    11    01    0.0     0.   0.    0.    0.    0.     0.    0.    0.  ; Bottom fixed in Z (free-slip)
4         1    11    01   -1.e-9   0.   0.    0.    0.    0.     0.    0.    0.  ; Top compression in Z
4         1    11    10    0.0     0.   0.    0.    0.    0.     0.    0.    0.  ; Top fixed in X (no-slip)
```

### Parameters Breakdown:
1. **Side Selection (`nofside`)**:
   * **`2`**: Represents the **bottom boundary** ($Z = -3000\text{ m}$).
   * **`4`**: Represents the **top boundary** ($Z = 0\text{ m}$).
2. **Node Range (`nbc1` to `nbc2`)**:
   * The grid is $10 \times 30$ elements, giving 11 nodes along the horizontal axis. Specifying nodes **`1` to `11`** applies the boundary condition continuously across the entire width of the domain.
3. **Boundary Condition Type (`nbc`)**:
   * **`01`**: Specifies a **vertical velocity ($V_z$) constraint** in meters/second.
   * **`10`**: Specifies a **horizontal velocity ($V_x$) constraint** in meters/second.
4. **Boundary Condition Profiles**:
   * **Bottom Boundary**: Vertically fixed ($V_z = 0$) but free to slide horizontally (free-slip).
   * **Top Boundary**: Compressing vertically downward at a constant velocity of $V_z = -1.0 \times 10^{-9}\text{ m/s}$ and fixed horizontally at $V_x = 0$ (no-slip condition).
   * **Lateral Boundaries ($X = 0$ m and $X = 1000$ m, Sides 1 and 3)**: Completely free of traction ($\sigma_{xx} = 0.0$, $\sigma_{xz} = 0.0$).

---

## 4. Analytical Formulation (Mohr-Coulomb Yield Criterion)

The Mohr-Coulomb yield criterion describes the shear strength of rocks and soils in terms of normal and shear stresses on the failure plane:

$$\tau = c + \sigma'_n \tan \phi$$

In terms of principal stresses $\sigma_1 \ge \sigma_2 \ge \sigma_3$ (using the classical soil mechanics convention where compressive stresses are positive, so $p_1 \ge p_2 \ge p_3 \ge 0$):

$$p_1 = p_3 N_\phi + 2 c \sqrt{N_\phi}$$

where $N_\phi$ is the passive pressure coefficient:

$$N_\phi = \tan^2\left( 45^\circ + \frac{\phi}{2} \right) = \frac{1 + \sin \phi}{1 - \sin \phi}$$

### Conversion to GeoFLAC Sign Convention (Tension Positive)
In GeoFLAC, compressive stresses are negative. Let $\sigma_1$ be the maximum principal stress (least compressive, horizontal stress $\sigma_{xx}$) and $\sigma_3$ be the minimum principal stress (most compressive, vertical stress $\sigma_{zz}$). 

Substituting $p_1 = -\sigma_3$ (most compressive) and $p_3 = -\sigma_1$ (least compressive) yields the yield function $f = 0$ in principal stress space:

$$\sigma_3 = \sigma_1 N_\phi - 2 c \sqrt{N_\phi}$$

### Unconfined Compressive Strength (UCS)
Since the lateral boundaries are completely free and unconfined, the horizontal stress remains zero throughout the test ($\sigma_{xx} = \sigma_1 = 0.0$). Substituting $\sigma_1 = 0$ into the yield equation gives the **Unconfined Compressive Strength (UCS)** of the rock, denoted as $\sigma_C$:

$$\sigma_{zz}^{yield} = \sigma_3 = -2 c \sqrt{N_\phi} = -\sigma_C$$

where:
$$\sigma_C = 2 c \sqrt{\frac{1 + \sin \phi}{1 - \sin \phi}} = 2 c \tan\left( 45^\circ + \frac{\phi}{2} \right)$$

### Strain Softening Analytical path
With strain softening, cohesion decreases linearly from its peak value $c_1$ to its residual value $c_2$ as the accumulated plastic strain ($\epsilon_p$) increases from $\text{pls1}$ to $\text{pls2}$:
$$c(\epsilon_p) = c_1 + (c_2 - c_1) \frac{\epsilon_p}{\text{pls2} - \text{pls1}}$$

This yields a peak and residual unconfined compressive strength:
* **Peak Strength ($\sigma_C^{peak}$)**:
  $$\sigma_C^{peak} = 2 c_1 \tan\left( 45^\circ + \frac{\phi}{2} \right) = 2 \times (20 \text{ MPa}) \times \sqrt{3} \approx 69.282 \text{ MPa}$$
* **Residual Strength ($\sigma_C^{residual}$)**:
  $$\sigma_C^{residual} = 2 c_2 \tan\left( 45^\circ + \frac{\phi}{2} \right) = 2 \times (2 \text{ MPa}) \times \sqrt{3} \approx 6.928 \text{ MPa}$$

During yielding, the total vertical strain $\epsilon_{zz}$ is the sum of elastic and plastic components ($\epsilon_{zz} = \epsilon_{zz}^e - \epsilon_p$). Since elastic strain is $\epsilon_{zz}^e = -\sigma_C(\epsilon_p) / E_{eff}$, we can derive the homogeneous stress-strain softening path analytically:
$$\epsilon_{zz} = -\frac{\sigma_C(\epsilon_p)}{E_{eff}} - \epsilon_p = -\frac{\sigma_C^{peak}}{E_{eff}} - \epsilon_p \left[ 1 + \frac{\sigma_C^{residual} - \sigma_C^{peak}}{E_{eff} \cdot \text{pls2}} \right]$$

Thus, the vertical column deforms elastically (with effective modulus $E_{eff} = 80$ GPa) until reaching the peak stress of **$-69.282$ MPa** at strain $\epsilon_{zz}^{yield} \approx -0.000866$, after which the stress weakens progressively toward the residual plateau.

*Note: Total vertical stress is reconstructed by summing deviatoric stress and pressure: $\sigma_{zz} = \sigma'_{zz} + P$. Please refer to the Elastic tutorial for a detailed breakdown of stress decomposition and Kilobar-to-MPa conversion in GeoFLAC.*

---

## 5. Simulation Results

### Stress and Strain Evolution
The simulation captures the transition from elastic loading to progressive plastic softening with exceptional detail.

* **Confinement Effect**:
  Due to the horizontal no-slip condition ($V_x = 0$) at the top boundary, the material is restricted from expanding horizontally at the contact surface. This restriction acts as a local confinement, which increases the effective compressive strength. As a result, the peak vertical stress in the simulation increases to **$-74.49$ MPa** (Frame 7), exceeding the analytical unconfined compressive strength (UCS) of **$-69.282$ MPa** (which assumes completely free lateral boundaries without any horizontal restriction).
  
* **Yielding/Failure Nucleation**:
  Under vertical compression, the material naturally expands horizontally due to Poisson's effect.
  * At the bottom boundary, nodes are free to slide horizontally (free-slip, $V_x$ is unconstrained), so the base expands freely without stress concentration.
  * At the top boundary, however, the no-slip condition ($V_x = 0$) prevents this lateral expansion. This constraint forces the material near the top contact to remain stationary in X while being pushed downward in Z, creating a severe shear stress gradient and concentration at the top corners.
  * Consequently, yielding (plastic strain `aps`) initiates in the top corner elements (Row 1, columns 1 and 10) rather than the bottom, and then propagates downward as localized shear bands.

* **Stress Softening and Localization**: 
  After reaching the peak of **$-74.49$ MPa** (Frame 7), the average vertical stress softens progressively to **$-66.20$ MPa** (Frame 10). This structural softening is driven by **strain localization and shear banding**: plastic strain concentrates in narrow localized bands. Elements within the shear bands accumulate plastic strain much faster than the homogeneous average, which accelerates their softening and reduces the overall load-bearing capacity of the vertical column.

### Verification Chart
The generated plot is saved to `images/stress_strain_yield.png`:

![Stress-Strain Yielding Verification Chart](images/stress_strain_yield.png)

1. **Blue Line**: The simulation stress-strain loading path with strain softening under top confinement.
2. **Red Dashed Line**: The theoretical analytical Mohr-Coulomb path showing homogeneous elastic loading and softening (unconfined).

The simulation accurately captures the transition from elastic loading to progressive plastic softening, verifying both the correct implementation and numerical stability of GeoFLAC's strain-softening mechanics.
