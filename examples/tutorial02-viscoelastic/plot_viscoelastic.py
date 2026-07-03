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
        t_myr = fl.time[i-1]
        t_kyr = t_myr * 1000.0
        
        sxx_dev = fl.read_sxx(i)
        pres = fl.read_pres(i)
        sxx_total = sxx_dev + pres
        
        mean_sxx_kb = np.mean(sxx_total)
        mean_sxx_mpa = mean_sxx_kb * 100.0
        
        times_myr.append(t_myr)
        stresses_mpa.append(mean_sxx_mpa)
        
        print(f"{i}\t{t_kyr:.4f}\t\t{mean_sxx_mpa:.4f}")

    # Physical parameters for Maxwell Viscoelastic Model
    eta = 1.0e21          # Viscosity in Pa.s
    lambda_lame = 3.0e10  # Lamé parameter lambda in Pa
    mu_lame = 3.0e10      # Shear modulus mu in Pa
    
    # Grid parameters
    L = 20000.0           # Width of the domain in m
    Vx = -0.001           # Compression velocity in m/yr (0.1 cm/yr compressive)
    sec_year = 3.1558e7   # Seconds in a year
    
    # Strain rate in s^-1
    strain_rate_s = (Vx / L) / sec_year  # approx -1.584384e-15 s^-1
    
    # Effective relaxation time
    tau_eff_s = eta * (lambda_lame + 2.0 * mu_lame) / (mu_lame * (lambda_lame + mu_lame))
    tau_eff_yr = tau_eff_s / sec_year
    
    # Analytical steady-state stress under compression
    sigma_xx_steady_pa = 4.0 * eta * strain_rate_s
    sigma_xx_steady_mpa = sigma_xx_steady_pa / 1.0e6  # approx -6.3375 MPa
    
    # High-resolution time array for two-stage analytical curve (in Kyr)
    t_analytical_kyr = np.linspace(0.0, 20.0, 1000)
    sxx_analytical_mpa = []
    
    # Stress at 10.0 Kyr (end of stage 1 compression)
    t_switch_kyr = 10.0
    t_switch_s = t_switch_kyr * 1000.0 * sec_year
    sigma_at_switch = sigma_xx_steady_mpa * (1.0 - np.exp(-t_switch_s / tau_eff_s))
    
    for t_kyr in t_analytical_kyr:
        if t_kyr <= t_switch_kyr:
            # Stage 1: Viscoelastic loading/compression
            t_s = t_kyr * 1000.0 * sec_year
            s = sigma_xx_steady_mpa * (1.0 - np.exp(-t_s / tau_eff_s))
        else:
            # Stage 2: Pure Maxwell relaxation (strain rate = 0)
            t_relax_s = (t_kyr - t_switch_kyr) * 1000.0 * sec_year
            s = sigma_at_switch * np.exp(-t_relax_s / tau_eff_s)
        sxx_analytical_mpa.append(s)
        
    sxx_analytical_mpa = np.array(sxx_analytical_mpa)

    # Plot results
    plt.figure(figsize=(9, 6.5))
    
    # Plot high-resolution analytical solution
    plt.plot(t_analytical_kyr, sxx_analytical_mpa, '--', color='#d62728', linewidth=2.0, 
             label='Analytical Solution (Two-Stage)')
    
    # Plot simulation points
    times_kyr = np.array(times_myr) * 1000.0
    plt.plot(times_kyr, stresses_mpa, 'o', color='#1f77b4', markersize=8, label='GeoFLAC Simulation')
    plt.plot(times_kyr, stresses_mpa, '-', color='#1f77b4', linewidth=1.5, alpha=0.7)

    # Add vertical line indicating the restart / relaxation boundary
    plt.axvline(x=10.0, color='gray', linestyle=':', linewidth=1.5)
    plt.text(10.2, sigma_xx_steady_mpa * 0.95, 'Velocity set to 0.0 (Restart)', fontsize=10.5, color='gray', rotation=90, va='bottom')

    plt.title(r'Viscoelastic Stress Relaxation: Loading (0-10 Kyr) & Relaxation (10-20 Kyr)', fontsize=13, fontweight='bold', pad=15)
    plt.xlabel('Time (Kyr)', fontsize=12, labelpad=10)
    plt.ylabel(r'Total Stress $\sigma_{xx}$ (MPa)', fontsize=12, labelpad=10)
    
    plt.xlim(-1.0, 21.0)
    plt.ylim(sigma_xx_steady_mpa * 1.25, 0.5)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=11, loc='lower left')
    
    # Adjust layout and save figure
    plt.tight_layout()
    os.makedirs('images', exist_ok=True)
    plot_path = 'images/stress_relaxation.png'
    plt.savefig(plot_path, dpi=300)
    print(f"\nPlot saved successfully to '{plot_path}'")

if __name__ == '__main__':
    main()
