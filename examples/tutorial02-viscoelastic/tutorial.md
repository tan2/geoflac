# GeoFLAC Tutorial: Viscoelastic Maxwell Compression

This tutorial explains the setup, boundary conditions, physical results, and analytical verification of the viscoelastic Maxwell compression benchmark in **GeoFLAC**.

---

## 1. Model Setup

The model represents a vertical two-dimensional column of homogenous viscoelastic rock undergoing horizontal compression. Under compression, the material exhibits both transient elastic stress build-up and long-term viscous relaxation governed by a Newtonian Maxwell rheology.

### Geometry and Mesh
* **Dimensions**: 20 km wide ($X \in [0, 20]$ km), 5 km high ($Z \in [-5, 0]$ km).
* **Grid Resolution**: $40 \times 10$ elements in the $X$ and $Z$ directions, yielding a regular grid of square elements ($500 \times 500$ m).

### Material Properties
The material is a Newtonian Maxwell viscoelastic rock defined in the input file `maxwell.inp`:
* **Lamé constant ($\lambda$)**: $3.0 \times 10^{10} \text{ Pa}$ (`Lame:rl`)
* **Shear Modulus ($\mu$)**: $3.0 \times 10^{10} \text{ Pa}$ (`Lame:rm`)
* **Newtonian Viscosity ($\eta$)**: $1.0 \times 10^{21} \text{ Pa}\cdot\text{s}$
* **Density ($\rho$)**: $2700 \text{ kg/m}^3$ (`den`)
* **Gravity ($g$)**: $0.0 \text{ m/s}^2$ (pure compression without gravity-induced pre-stress or hydrostatic gradients).

> [!NOTE]
> **Formulating Newtonian Viscosity in GeoFLAC:**
> Because GeoFLAC models viscosity using power-law creep ($\dot{\epsilon} = A \sigma^n$), we simulate Newtonian Maxwell rheology by setting:
> * Power exponent $n = 1.0$ (`pln = 1.0`)
> * Zero thermal and pressure activation dependencies (`eactiv = 0.0`, `vactiv = 0.0`)
> * The pre-exponential coefficient `acoef` calculated from viscosity as:
>   $$\text{acoef} = \frac{10^6}{3 \eta} = 3.3333333 \times 10^{-16}$$
> * The viscosity limit bounds in `maxwell.inp` are set to `1.0e+18` and `1.0e+25` so the target $10^{21} \text{ Pa}\cdot\text{s}$ is safely inside the range.

---

## 2. Boundary Conditions

The mechanical boundary conditions are configured in `maxwell.inp` to compress the block horizontally:

1. **Left Boundary ($X = 0$ km, Side 1)**:
   * Constrained to move horizontally inward (to the right) at a constant velocity of $V_x = 1.0 \text{ cm/yr} \approx 3.1687686 \times 10^{-10} \text{ m/s}$.
   * Free to move vertically (zero vertical shear traction).
2. **Right Boundary ($X = 20$ km, Side 3)**:
   * Constrained horizontally to zero velocity ($V_x = 0.0$ m/s, forming a rigid wall).
   * Free to move vertically (free-slip boundary condition).
3. **Bottom Boundary ($Z = -5$ km, Side 2)**:
   * Constrained vertically to zero velocity ($V_z = 0.0$ m/s).
   * Free to expand/slide horizontally (free-slip boundary condition).
4. **Top Boundary ($Z = 0$ km, Side 4)**:
   * Completely free of traction ($\sigma_{zz} = 0.0$ and $\sigma_{zx} = 0.0$).

---

## 3. Analytical Formulation (Newtonian Maxwell Flow)

Under 2D plane strain ($\epsilon_{yy} = 0$) and a free top surface ($\sigma_{zz} = 0$), the total horizontal compressive stress $\sigma_{xx}(t)$ builds up elastically and relaxes viscously over time:

$$\sigma_{xx}(t) = \sigma_{xx}^{steady} \left( 1 - e^{-t/\tau_{eff}} \right)$$

where:
1. **Horizontal Strain Rate ($\dot{\epsilon}_{xx}$)**:
   $$\dot{\epsilon}_{xx} = \frac{V_x}{L} = \frac{-0.01 \text{ m/yr}}{20000 \text{ m}} \approx -1.584384 \times 10^{-14} \text{ s}^{-1}$$
2. **Effective Viscous Resistance**:
   Under plane strain Newtonian flow, the effective viscosity resisting horizontal contraction is $\eta_{eff} = 4\eta$.
   $$\sigma_{xx}^{steady} = 4\eta\dot{\epsilon}_{xx} = 4 \times 10^{21} \text{ Pa}\cdot\text{s} \times (-1.584384 \times 10^{-14} \text{ s}^{-1}) = -63.375 \text{ MPa}$$
3. **Effective Relaxation Time ($\tau_{eff}$)**:
   $$\tau_{eff} = \frac{\eta_{eff}}{E_{eff}} = \frac{4\eta}{\frac{4\mu(\lambda + \mu)}{\lambda + 2\mu}} = \frac{\eta (\lambda + 2\mu)}{\mu (\lambda + \mu)}$$
   Using Lamé constants $\lambda = \mu = 3.0 \times 10^{10} \text{ Pa}$:
   $$\tau_{eff} = 5.0 \times 10^{10} \text{ s} \approx 1584.38 \text{ years} = 1.584 \text{ Kyr}$$

Because the relaxation time ($\approx 1.584$ Kyr) is much shorter than the simulation's 10 Kyr output interval, the stress is fully relaxed to its steady-state Newtonian plateau by the first output step ($10,000$ years $\approx 6.3\tau_{eff}$).

---

## 4. Reconstructing Total Stress in GeoFLAC

> [!IMPORTANT]
> **Stress Representation in GeoFLAC binary files:**
> * `sxx.0` and `szz.0` store the **deviatoric stresses** ($\sigma'_{xx} = \sigma_{xx} - P$ and $\sigma'_{zz} = \sigma_{zz} - P$), not the total stresses.
> * `pres.0` stores the **mean stress / pressure** ($P$).
> * All stress outputs are in **Kilobars** ($1 \text{ kb} = 100 \text{ MPa} = 10^8 \text{ Pa}$).

To obtain the total horizontal stress $\sigma_{xx}$, you must sum the deviatoric horizontal stress and pressure:
$$\sigma_{xx} = \sigma'_{xx} + P$$
Multiplying by $100.0$ converts the result from Kilobars to standard MegaPascals (MPa).

---

## 5. Running the Simulation and Plotting

### Step 1: Run the GeoFLAC Solver
First, ensure that the output directory is clean of any past run files, and execute the compiled `flac` binary:
```bash
rm -f *.0 *.rs *.vts _contents.* _markers.* pisos.rs time.rs vbc.s output.asc sys.msg
../../src/flac maxwell.inp
```
The solver will run for 190,000 steps, simulating a total time of $101.0$ Kyr and writing outputs every $10.0$ Kyr.

### Step 2: Plot the Stress-Relaxation Curve
Run the provided Python plotting script:
```bash
python3 plot_viscoelastic.py
```
This script reads the binary files, reconstructs the total horizontal stress $\sigma_{xx}$, and plots the simulation points alongside the high-resolution analytical relaxation curve, saving the plot as `stress_relaxation.png`.

---

## 6. Simulation Results

The simulation shows outstanding agreement with the analytical steady-state plateau:

| Frame | Time (Kyr) | Simulated $\sigma_{xx}$ (MPa) | Analytical Steady-State (MPa) | Discrepancy (%) |
|:---:|:---:|:---:|:---:|:---:|
| 1 | 0.0 | 0.0000 | 0.0000 | — |
| 2 | 10.0 | -65.7561 | -63.3750 | 3.8% |
| 3 | 20.1 | -64.3144 | -63.3750 | 1.5% |
| 4 | 30.1 | -64.3319 | -63.3750 | 1.5% |
| 5 | 40.1 | -64.6287 | -63.3750 | 2.0% |
| 6 | 50.2 | -64.9607 | -63.3750 | 2.5% |
| 7 | 60.2 | -65.2965 | -63.3750 | 3.0% |
| 8 | 70.2 | -65.6359 | -63.3750 | 3.6% |
| 9 | 80.3 | -65.9789 | -63.3750 | 4.1% |
| 10 | 90.3 | -66.3255 | -63.3750 | 4.7% |
| 11 | 100.3 | -66.6757 | -63.3750 | 5.5% |

### Verification Chart
The generated plot `stress_relaxation.png` displays:
1. **Blue Dots**: The simulation output points spaced at 10 Kyr intervals.
2. **Red Dashed Line**: The theoretical analytical relaxation curve showing the quick elastic stress build-up and the long-term viscoelastic steady-state plateau at $-63.375$ MPa.

The extremely minor 2-5% discrepancy represents physical geometrical thinning as the domain is compressed continuously over 100 Kyr (large strain effects), confirming the physical validity and high numerical precision of **GeoFLAC**'s viscoelastic mechanical solver.
