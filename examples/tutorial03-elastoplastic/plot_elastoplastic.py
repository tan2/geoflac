#!/usr/bin/env python3
import sys
import os

# Add geoflac utilities path
sys.path.append('../../util')
import flac
import numpy as np

try:
    import matplotlib
    # Use non-interactive Agg backend to avoid GUI window popup issues in terminal environments
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is required to plot. Please install it.")
    sys.exit(1)

def main():
    # Instantiate Flac reader
    fl = flac.Flac()
    
    strains = []
    stresses = []
    
    print("Frame\tTime (Kyr)\tStrain (zz)\tStress (zz, MPa)")
    for i in range(1, fl.nrec + 1):
        t_kyr = fl.time[i-1] * 1000.0
        
        # Read strain components: exx, ezz, exz
        _, ezz, _ = fl.read_strain(i)
        
        # GeoFLAC outputs deviatoric stresses to szz.0, and pressure (mean stress) to pres.0.
        # Both are stored in Kilobars (1 kb = 100 MPa = 10^8 Pa).
        szz_dev = fl.read_szz(i)
        pres = fl.read_pres(i)
        
        # Total vertical stress = deviatoric vertical stress + mean stress (pressure)
        szz_total = szz_dev + pres
        
        # Average over the entire model domain
        mean_ezz = np.mean(ezz)
        mean_szz_kb = np.mean(szz_total)
        
        # Convert stress to MPa (1 Kilobar = 100 MPa)
        mean_szz_mpa = mean_szz_kb * 100.0
        
        strains.append(mean_ezz)
        stresses.append(mean_szz_mpa)
        
        print(f"{i}\t{t_kyr:.2f}\t\t{mean_ezz:g}\t\t{mean_szz_mpa:.2f}")

    strains = np.array(strains)
    stresses = np.array(stresses)

    # Physical parameters for Mohr-Coulomb model with strain weakening
    cohesion1 = 20.0e6      # Peak cohesion in Pa (20 MPa)
    cohesion2 = 2.0e6       # Residual cohesion in Pa (2 MPa)
    friction_angle_deg = 30.0  # Friction angle in degrees
    friction_angle_rad = np.radians(friction_angle_deg)
    
    plstrain1 = 0.0
    plstrain2 = 0.1
    
    # Unconfined Compressive Strength (UCS) analytical solution:
    # N_phi = (1 + sin(phi)) / (1 - sin(phi)) = tan^2(45 + phi/2)
    N_phi = (1.0 + np.sin(friction_angle_rad)) / (1.0 - np.sin(friction_angle_rad))
    
    # Peak and Residual UCS:
    ucs_peak_mpa = 2.0 * cohesion1 * np.sqrt(N_phi) / 1.0e6       # 69.282 MPa
    ucs_residual_mpa = 2.0 * cohesion2 * np.sqrt(N_phi) / 1.0e6   # 6.928 MPa
    
    # Effective vertical modulus under plane strain (80 GPa)
    E_eff = 80000.0  # in MPa
    
    # Generate analytical stress-strain path under strain weakening
    # Elastic yield strain limit:
    yield_strain = -ucs_peak_mpa / E_eff  # approx -0.000866
    
    # Slope factor S for plastic strain vs total strain:
    # eps_zz = eps_zz_e + eps_zz_p = -sigma_C(aps)/E_eff - aps
    # S = 1 + (sigma_C_residual - sigma_C_peak) / (E_eff * (plstrain2 - plstrain1))
    S = 1.0 + (ucs_residual_mpa - ucs_peak_mpa) / (E_eff * (plstrain2 - plstrain1))
    
    analytical_strains = np.linspace(min(strains), 0.0, 1000)
    analytical_stresses = []
    
    for eps in analytical_strains:
        if eps > yield_strain:
            # Elastic phase
            analytical_stresses.append(E_eff * eps)
        else:
            # Homogeneous plastic softening phase
            # Solve for aps: eps = yield_strain - aps * S  => aps = (eps - yield_strain) / -S
            aps = (eps - yield_strain) / -S
            if aps < plstrain2:
                ucs_val = ucs_peak_mpa + (ucs_residual_mpa - ucs_peak_mpa) * (aps / (plstrain2 - plstrain1))
                analytical_stresses.append(-ucs_val)
            else:
                analytical_stresses.append(-ucs_residual_mpa)
                
    analytical_stresses = np.array(analytical_stresses)
 
    # Plot results
    plt.figure(figsize=(9, 6.5))
    
    # Plot analytical yield curve
    plt.plot(analytical_strains, analytical_stresses, '--', color='#d62728', linewidth=2.0, 
             label=f'Homogeneous Softening Analytical path')
    
    # Plot simulation points and line
    plt.plot(strains, stresses, 'o-', color='#1f77b4', linewidth=2.0, markersize=8, label='GeoFLAC Simulation')
 
    # Add label details
    plt.text(min(strains) * 0.7, -ucs_peak_mpa * 0.4, 
             f"Mohr-Coulomb Parameters:\n"
             f"Peak Cohesion ($c_1$) = {cohesion1/1e6:.1f} MPa\n"
             f"Residual Cohesion ($c_2$) = {cohesion2/1e6:.1f} MPa\n"
             f"Friction angle ($\\phi$) = {friction_angle_deg:.1f}$^\\circ$\n"
             f"Effective Modulus ($E_{{eff}}$) = {E_eff/1e3:.1f} GPa\n"
             f"Peak UCS ($\\sigma_C^{{peak}}$) = {ucs_peak_mpa:.2f} MPa\n"
             f"Residual UCS ($\\sigma_C^{{residual}}$) = {ucs_residual_mpa:.2f} MPa",
             fontsize=10.5, bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))
 
    plt.title(r'Elastoplastic Compression: $\sigma_{zz}$ vs. Strain $\epsilon_{zz}$ (Mohr-Coulomb)', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel(r'Vertical Strain ($\epsilon_{zz}$)', fontsize=12, labelpad=10)
    plt.ylabel(r'Total Stress $\sigma_{zz}$ (MPa)', fontsize=12, labelpad=10)
    
    plt.xlim(min(strains) * 1.1, 0.0002)
    plt.ylim(-ucs_peak_mpa * 1.25, 5.0)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=11, loc='best')
    
    # Adjust layout and save figure
    plt.tight_layout()
    os.makedirs('images', exist_ok=True)
    plot_path = 'images/stress_strain_yield.png'
    plt.savefig(plot_path, dpi=300)
    print(f"\nPlot saved successfully to '{plot_path}'")

if __name__ == '__main__':
    main()
