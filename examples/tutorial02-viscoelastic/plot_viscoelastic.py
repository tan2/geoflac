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
    
    times_myr = []
    stresses_mpa = []
    
    print("Frame\tTime (Kyr)\tStress (xx, MPa)")
    for i in range(1, fl.nrec + 1):
        # Time in Myr
        t_myr = fl.time[i-1]
        t_kyr = t_myr * 1000.0
        
        # GeoFLAC outputs deviatoric stresses to sxx.0, and pressure (mean stress) to pres.0.
        # Both are stored in Kilobars (1 kb = 100 MPa = 10^8 Pa).
        sxx_dev = fl.read_sxx(i)
        pres = fl.read_pres(i)
        
        # Total horizontal stress = deviatoric horizontal stress + mean stress (pressure)
        sxx_total = sxx_dev + pres
        
        # Average over the entire model domain
        mean_sxx_kb = np.mean(sxx_total)
        
        # Convert stress to MPa (1 Kilobar = 100 MPa)
        mean_sxx_mpa = mean_sxx_kb * 100.0
        
        times_myr.append(t_myr)
        stresses_mpa.append(mean_sxx_mpa)
        
        print(f"{i}\t{t_kyr:.1f}\t\t{mean_sxx_mpa:.4f}")

    # Physical parameters for Maxwell Viscoelastic Model
    eta = 1.0e21          # Viscosity in Pa.s
    lambda_lame = 3.0e10  # Lamé parameter lambda in Pa
    mu_lame = 3.0e10      # Shear modulus mu in Pa
    
    # Grid parameters
    L = 20000.0           # Width of the domain in m
    Vx = -0.01            # Left boundary velocity in m/yr (1 cm/yr compressive)
    sec_year = 3.1558e7   # Seconds in a year
    
    # Strain rate in s^-1
    strain_rate_s = (Vx / L) / sec_year  # approx -1.584384e-14 s^-1
    
    # Effective 2D plane strain elastic modulus (E_eff = E / (1 - nu^2))
    # E_eff = 4 * mu * (lambda + mu) / (lambda + 2 * mu) = 80 GPa = 80,000 MPa
    E_eff = 4.0 * mu_lame * (lambda_lame + mu_lame) / (lambda_lame + 2.0 * mu_lame)
    
    # Effective relaxation time: tau_eff = eta_eff / E_eff = 4 * eta / E_eff
    # tau_eff = eta * (lambda + 2*mu) / (mu * (lambda + mu)) = 5.0e9 seconds (~158.438 years)
    tau_eff_s = eta * (lambda_lame + 2.0 * mu_lame) / (mu_lame * (lambda_lame + mu_lame))
    tau_eff_yr = tau_eff_s / sec_year
    tau_eff_kyr = tau_eff_yr / 1000.0
    
    # Analytical steady-state stress: sigma_xx_steady = 4 * eta * strain_rate
    sigma_xx_steady_pa = 4.0 * eta * strain_rate_s
    sigma_xx_steady_mpa = sigma_xx_steady_pa / 1.0e6  # approx -6.3375 MPa
    
    # High-resolution time array for analytical curve (in Kyr)
    t_analytical_kyr = np.linspace(0.0, max(times_myr) * 1000.0, 1000)
    t_analytical_s = t_analytical_kyr * 1000.0 * sec_year
    
    # Analytical stress over time: sigma_xx(t) = sigma_xx_steady * (1 - exp(-t / tau_eff))
    sxx_analytical_mpa = sigma_xx_steady_mpa * (1.0 - np.exp(-t_analytical_s / tau_eff_s))

    # Plot results
    plt.figure(figsize=(9, 6.5))
    
    # Plot high-resolution analytical solution
    plt.plot(t_analytical_kyr, sxx_analytical_mpa, '--', color='#d62728', linewidth=2.0, 
             label=r'Analytical Solution ($\sigma_{xx}(t) = 4\eta\dot{\epsilon}_{xx}(1 - e^{-t/\tau_{eff}})$)')
    
    # Plot simulation points
    times_kyr = np.array(times_myr) * 1000.0
    plt.plot(times_kyr, stresses_mpa, 'o', color='#1f77b4', markersize=8, label='GeoFLAC Simulation')
    plt.plot(times_kyr, stresses_mpa, '-', color='#1f77b4', linewidth=1.5, alpha=0.7)

    plt.text(18.0, sigma_xx_steady_mpa * 0.7, 
             f"Analytical parameters:\n"
             f"$\\eta = 1.0\\times 10^{{21}}$ Pa$\\cdot$s\n"
             f"$\\dot{{\\epsilon}}_{{xx}} = {strain_rate_s*1e14:.4f}\\times 10^{{-14}}$ s$^{{-1}}$\n"
             f"$\\tau_{{eff}} = {tau_eff_yr:.2f}$ years\n"
             f"$\\sigma_{{xx}}^{{steady}} = {sigma_xx_steady_mpa:.4f}$ MPa",
             fontsize=10.5, bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5'))

    plt.title(r'Viscoelastic Stress Relaxation: Total $\sigma_{xx}$ over Time', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Time (Kyr)', fontsize=12, labelpad=10)
    plt.ylabel(r'Total Stress $\sigma_{xx}$ (MPa)', fontsize=12, labelpad=10)
    
    plt.xlim(-2.0, 42.0)
    # Stress is compressive (negative), so we set the limits accordingly
    plt.ylim(sigma_xx_steady_mpa * 1.25, 0.5)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=11, loc='lower left')
    
    # Adjust layout and save figure
    plt.tight_layout()
    plot_path = 'stress_relaxation.png'
    plt.savefig(plot_path, dpi=300)
    print(f"\nPlot saved successfully to '{plot_path}'")

if __name__ == '__main__':
    main()
