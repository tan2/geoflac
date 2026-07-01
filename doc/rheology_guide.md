# GeoFLAC User Guide: Rheology Reference

This guide provides a comprehensive reference explaining how physical rock properties are represented in GeoFLAC's input files (`.inp`) and how the solver calculates elastic, plastic, and viscous stresses.

---

## 1. Elasticity

GeoFLAC uses isotropic linear elasticity governed by Hooke's Law. In the input file, elastic properties are defined using the **Lamé parameters**:

*   **`rl` ($\lambda$)**: The first Lamé parameter (in Pa).
*   **`rm` ($\mu$)**: The shear modulus (in Pa).

### Relationship to Engineering Moduli
In laboratory measurements or rock mechanics, elastic properties are more commonly specified by **Young's Modulus ($E$)** and **Poisson's Ratio ($\nu$)**, or the **Bulk Modulus ($K$)**. You can convert these properties to GeoFLAC's Lamé parameters using the following formulas:

$$\lambda = \frac{E \nu}{(1 + \nu)(1 - 2\nu)}$$

$$\mu = \frac{E}{2(1 + \nu)}$$

Conversely, to calculate Young's Modulus and Poisson's Ratio from GeoFLAC's input parameters:

$$E = \frac{\mu (3\lambda + 2\mu)}{\lambda + \mu}$$

$$\nu = \frac{\lambda}{2(\lambda + \mu)}$$

The **Bulk Modulus ($K$)** (which controls volume changes under pressure) is calculated as:

$$K = \lambda + \frac{2}{3}\mu$$

> [!NOTE]
> Under geological time scales, effective elastic parameters are often set lower (by up to an order of magnitude) than those derived from short-term seismic wave propagation, to account for long-term crack damage and micro-fractures in the crust.

---

## 2. Mohr-Coulomb Plasticity & Strain Softening

Plastic deformation is modeled using a non-associated Mohr-Coulomb flow rule with linear strain softening (`irheol = 6` or `12`). 

### The Yield Criterion
An element yields plastically when its stress state reaches the Mohr-Coulomb failure envelope:

$$\sigma_s = C + \sigma_n \tan \phi$$

Where:
*   $\sigma_s$ is the shear stress.
*   $\sigma_n$ is the normal stress (confining pressure, where compression is negative).
*   $C$ is the **cohesion** (in Pa).
*   $\phi$ is the **friction angle** (in degrees).

### Dilation and Plastic Flow
The direction of plastic deformation is governed by the **dilation angle ($\psi$)**:
*   $\psi = 0^\circ$: Zero plastic volume change during shear.
*   $\psi > 0^\circ$: The material expands (dilates) during shear.
*   **Non-associated flow rule**: In GeoFLAC, $\psi$ is typically set much smaller than the friction angle $\phi$ (e.g., $\psi \le 5^\circ$ while $\phi = 30^\circ$) to reflect realistic rock behavior, preventing unphysical runaway volume expansion.

### Linear Strain Softening
As a fault zone accommodates slip, the material weakens. This is simulated by linearly interpolating the cohesion, friction, and dilation angles between their peak (pre-yield) values and their residual values based on the **accumulated plastic strain ($\epsilon_p$ or `aps`)**:

| Parameter in `.inp` | Description |
|:--------------------|:------------|
| **`plstrain1`** ($\epsilon_{p1}$) | Plastic strain at which softening **begins**. |
| **`plstrain2`** ($\epsilon_{p2}$) | Plastic strain at which softening **saturates** (reaches residual). |
| **`cohesion1`, `cohesion2`** | Peak and residual Cohesion ($C_1, C_2$ in Pa). |
| **`fric1`, `fric2`** | Peak and residual Friction Angle ($\phi_1, \phi_2$ in degrees). |
| **`dilat1`, `dilat2`** | Peak and residual Dilation Angle ($\psi_1, \psi_2$ in degrees). |

#### Softening Interpolation Logic:
For any element with an accumulated plastic strain $\epsilon_p$:
1.  **Before Softening** ($\epsilon_p \le \epsilon_{p1}$): Properties are at peak values ($C_1, \phi_1, \psi_1$).
2.  **During Softening** ($\epsilon_{p1} < \epsilon_p < \epsilon_{p2}$): Properties are linearly interpolated:
    $$C(\epsilon_p) = C_1 + (C_2 - C_1) \times \frac{\epsilon_p - \epsilon_{p1}}{\epsilon_{p2} - \epsilon_{p1}}$$
3.  **Saturated / Residual** ($\epsilon_p \ge \epsilon_{p2}$): Properties remain constant at their residual values ($C_2, \phi_2, \psi_2$).

---

## 3. Viscous Creep (Non-Newtonian Power-Law)

High-temperature deformation in the lower crust and mantle is dominated by ductile viscous flow (`irheol = 12`). GeoFLAC calculates the effective non-Newtonian viscosity ($\eta$) at each element using the standard dislocation creep formulation:

$$\eta = 0.25 \times 10^6 \times (0.75 \cdot A)^{-1/n} \cdot \dot{\epsilon}_{II}^{(1/n - 1)} \cdot \exp \left( \frac{E + P \cdot V}{n R (T + 273.15)} \right)$$

Where:

| Input Parameter | Physical Quantity | Standard Units |
|:----------------|:------------------|:---------------|
| **`acoef`** ($A$) | Power-law pre-exponential coefficient | $\text{Pa}^{-n} \text{ s}^{-1}$ |
| **`pln`** ($n$) | Power-law stress exponent (typically 3 to 4) | Dimensionless |
| **`eactiv`** ($E$) | Activation energy | $\text{J} / \text{mol}$ |
| **`vactiv`** ($V$) | Activation volume | $\text{m}^3 / \text{mol}$ (set to `0.0` if pressure-independent) |

Other variables calculated by the solver at runtime:
*   $\dot{\epsilon}_{II}$: Second invariant of the deviatoric strain rate tensor ($\text{s}^{-1}$).
*   $P$: Hydrostatic pressure ($\text{Pa}$).
*   $T$: Local temperature ($\text{°C}$).
*   $R$: Universal gas constant ($8.3144 \text{ J}/(\text{mol}\cdot\text{K})$).

### Viscosity Cutoffs (`v-min`, `v-max`)
To ensure numerical stability and prevent infinitely small timesteps in an explicit solver, viscosity is capped:
*   **`v-min`**: Minimum viscosity cutoff (typically $10^{19} - 10^{20} \text{ Pa}\cdot\text{s}$). If viscosity becomes too low, the Maxwell time step limit ($dt_{\text{maxwell}} = \eta / \mu$) forces the solver's timestep to become extremely small, slowing down the execution.
*   **`v-max`**: Maximum viscosity cutoff (typically $10^{26} - 10^{27} \text{ Pa}\cdot\text{s}$). Represents the transition to rigid elastic behavior.

---

## 4. Visco-Elasto-Plastic (VEP) Integration

In a VEP formulation (`irheol = 12`), the total strain rate tensor $\dot{\boldsymbol{\epsilon}}$ is decomposed into elastic, viscous, and plastic parts:

$$\dot{\boldsymbol{\epsilon}} = \dot{\boldsymbol{\epsilon}}_e + \dot{\boldsymbol{\epsilon}}_v + \dot{\boldsymbol{\epsilon}}_p$$

During each timestep:
1.  **Viscoelastic Trial**: Stresses are incremented using a viscoelastic Maxwell model (coupling elastic Lamé constants with the calculated non-Newtonian viscosity $\eta$).
2.  **Plastic Correction**: If the trial stress exceeds the Mohr-Coulomb yield stress ($\sigma_s > \sigma_s^{\text{yield}}$), the stress is projected back to the yield envelope, and the excess energy is recorded as a plastic strain increment ($d\epsilon_p$), which is added to the element's total `aps` accumulator.
