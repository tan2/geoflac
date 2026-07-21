# Material Phase Determination & Blending in geoflac

This document describes how material phases (lithologies) are initialized, dynamically updated, and blended to calculate effective physical properties in the `geoflac` engine.

---

## 1. Core Architecture: Markers vs. Elements

`geoflac` uses a hybrid Lagrangian-Eulerian approach (Active Marker-in-Cell method) to track material phases and historical properties:
*   **Lagrangian Markers (`mark_phase`)**: Individual tracking points carrying a specific phase ID. They move continuously with the velocity field.
*   **Eulerian Elements (`iphase` & `phase_ratio`)**: The grid cells where mechanical and thermal equations are solved. An element's material properties are determined by the collection of markers currently residing within it.

### Core Phase Arrays
| Array | Shape | Meaning |
| :--- | :--- | :--- |
| `mark_phase(max_markers)` | `(max_markers)` | Phase ID assigned to each individual marker. |
| `phase_ratio(nphase, nz-1, nx-1)` | `(nphase, nz-1, nx-1)` | Volume (number) fraction of each phase present in a grid element. |
| `iphase(nz-1, nx-1)` | `(nz-1, nx-1)` | The dominant material phase ID of a grid element. |

---

## 2. Element Phase & Ratio Calculation (`count_phase_ratio`)

The mapping from marker phases to grid element phases is performed by the subroutine `count_phase_ratio(j, i)` (implemented in `src/marker_data.f90`).

### Step-by-Step Logic
1.  **Safety Check**: Verifies that the element contains at least one marker (`nmark_elem(j,i) > 0`). If it is empty, the code halts with an error (`STOP 'No markers in element'`).
2.  **Phase Counting**: Counts the occurrences of each phase among the markers located inside the element:
    $$\text{ncounters}(k) = \sum_{n=1}^{N_{\text{markers}}} [ \text{mark\_phase}(\text{mark\_id\_elem}(n,j,i)) == k ]$$
3.  **Phase Ratio Calculation**: Computes the ratio (fraction) of each phase within the element:
    $$\text{phase\_ratio}(k, j, i) = \frac{\text{ncounters}(k)}{N_{\text{markers}}}$$
4.  **Dominant Phase Determination**: Iterates through all active phases (`1` to `nphase`) and selects the phase with the highest count to be `iphase(j,i)`. In case of a tie, the phase with the lower index is selected.

### Subroutine Implementation
```fortran
subroutine count_phase_ratio(j, i)
  !$ACC routine seq
  use params
  use arrays

  integer, intent(in) :: j, i
  integer :: n, kk, ncounters(maxph), iph, kph, nm

  if (nmark_elem(j,i) == 0) stop 'No markers in element'

  ncounters = 0
  do n = 1, nmark_elem(j,i)
    kk = mark_id_elem(n,j,i)
    iph = mark_phase(kk)
    ncounters(iph) = ncounters(iph) + 1
  enddo

  do kk = 1, nphase
    phase_ratio(kk,j,i) = ncounters(kk) / (nmark_elem(j,i) * 1.0d0)
  enddo

  ! The phase of this element is the most abundant marker phase
  iph = 1
  nm = 0
  do kk = 1, nphase
     if (ncounters(kk) > nm) then
         nm = ncounters(kk)
         iph = kk
     endif
  enddo
  iphase(j,i) = iph
end subroutine count_phase_ratio
```

---

## 3. Initial Phase Assignment

At the start of a simulation (in `src/init_marker.f90`), markers are uniformly distributed inside each element (typically 9 markers per element) and assigned initial phases based on the model geometry.

### 3.1. Layered Column Structure (`nzone_age`)
The model domain is divided horizontally into one or more columns/zones defined by grid range boundaries (`ixtb1` and `ixtb2`).
*   Within each zone `n`, layer boundaries are defined in depth kilometers (`hc(n, :)`).
*   For intermediate columns, layer depths can transition/interpolate horizontally:
    $$h_{\text{layer}}(x) = h_{\text{layer}}(n) + (h_{\text{layer}}(n+1) - h_{\text{layer}}(n)) \times \frac{x - x(n)}{x(n+1) - x(n)}$$
*   A marker at coordinates $(x, y)$ is assigned to a phase `kph` depending on which depth layer it falls into.

### 3.2. Initial Heterogeneities / Inclusions (`inhom`)
Initial heterogeneities (such as weak zones or inclusions) can override the layered setup using shapes:
*   **Rectangular Geometry (`igeom == 0`)**: Calls `newphase2marker` to assign all markers within a grid range $[j_1, j_2] \times [i_1, i_2]$ to the inclusion phase `inphase`.
*   **45-Degree Weak Zones (`igeom == 3` or `igeom == 4`)**: Resets markers along a diagonal band.

---

## 4. Dynamic Phase Changes (`change_phase`)

During a simulation, markers undergo physical and chemical phase transitions in `src/change_phase.f90` based on local temperature ($T$), pressure ($P$), depth, and structural configurations.

### 4.1. Serpentine Hydration / Dehydration
When subducting oceanic crust releases fluids, it hydrates the overlying mantle wedge to form serpentinite:
*   **Hydration Transition (`kmant1` / `kmant2` $\rightarrow$ `kserp`)**:
    *   Happens within a specific $P$-$T$ stability field taken from *Ulmer and Trommsdorff (Science, 1995)*.
    *   Occurs if the element is situated directly above subducting oceanic crust/sediment (`phase_ratio > 0.8`).
*   **Dehydration Transition (`kserp` $\rightarrow$ `khydmant`)**:
    *   Occurs if serpentinite is subducted beyond its pressure/temperature stability boundaries, dehydrating back into a hydrated mantle phase.

### 4.2. Basalt-to-Eclogite Transition (`koceandry`/`kocean1`/`kocean2` $\rightarrow$ `keclg`)
Oceanic basalt transitions to dense eclogite at depth:
*   Uses a parameterized stability boundary derived from *Hacker (JGR, 2003)*:
    $$P_{\text{transition}} = f(T)$$
*   Occurs when temperature $T > 400^\circ\text{C}$ and depth $z > 20\text{ km}$ if the pressure exceeds Hacker's transition pressure boundaries.

### 4.3. Sediment Subduction & Melting
*   Subducted continental sediments (`ksed1`, `ksed2`) can transform into `kmetased` and `kschist`.
*   If they subduct to hot melting zones ($T > T_{\text{solidus}}$), they generate dynamic melt fractions (`fmelt` and `fmagma`).

---

## 5. Blended Material Property Calculations

When solving mechanical and thermal equations, properties are often blended using either **Arithmetic Means** or **Harmonic Means** of all phases present in the element, weighted by `phase_ratio(k, j, i)`. This creates smooth material transitions across interfaces and prevents grid-step locking.

### 5.1. Effective Density (`Eff_dens` in `src/matprops.f90`)
Effective density is computed using an **Arithmetic Mean**:
$$\rho_{\text{eff}} = \sum_{k=1}^{N_{\text{phase}}} \left[ \text{phase\_ratio}(k, j, i) \times \rho_k(T, P) \right]$$
*   Where $\rho_k(T, P) = \rho_{0,k} \cdot \left( 1 - \alpha_k T + \beta_k P \right)$.
*   Phases with `phase_ratio < 0.01` are skipped.
*   The effective density is subsequently modified for magma presence using `fmagma`.

### 5.2. Effective Viscosity (`Eff_visc` in `src/matprops.f90`)
Effective viscosity is computed using a **Weighted Harmonic Mean** (which is physically appropriate for layers shearing in parallel):
$$\frac{1}{\eta_{\text{eff}}} = \sum_{k=1}^{N_{\text{phase}}} \frac{\text{phase\_ratio}(k, j, i)}{\eta_k(T, P, \dot{\epsilon}_{\text{II}})}$$
*   Where $\eta_k$ is the local power-law effective viscosity calculated for phase $k$.

### 5.3. Effective Plasticity Parameters (`src/rh_plastic.f90`)
For Mohr-Coulomb yielding, material parameters are blended differently based on their mechanical behaviors:
*   **Weighted Harmonic Mean** (ensures weak phases dominate yielding):
    *   **Friction Angle ($\phi$)**: $\frac{1}{\phi_{\text{eff}}} = \sum \frac{\text{phase\_ratio}(k)}{\phi_k}$
    *   **Cohesion ($C$)**: $\frac{1}{C_{\text{eff}}} = \sum \frac{\text{phase\_ratio}(k)}{C_k}$
*   **Weighted Arithmetic Mean**:
    *   **Dilation Angle ($\psi$)**: $\psi_{\text{eff}} = \sum \text{phase\_ratio}(k) \times \psi_k$
    *   **Hardening Parameter ($H$)**: $H_{\text{eff}} = \sum \text{phase\_ratio}(k) \times H_k$
