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
    
    print("Frame\tTime (Myr)\tStrain (zz)\tStress (zz, MPa)")
    for i in range(1, fl.nrec + 1):
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
        
        print(f"{i}\t{fl.time[i-1]:.4f}\t\t{mean_ezz:g}\t\t{mean_szz_mpa:.2f}")

    # Plot results
    plt.figure(figsize=(8, 6))
    plt.plot(strains, stresses, 'o-', color='#1f77b4', linewidth=2.5, markersize=8, label='Simulation')
    
    # Under 2D plane strain conditions (epsilon_yy = 0) and free lateral boundaries (sigma_xx = 0):
    # The Lame parameters are: lambda = 3.0e10 Pa, mu = 3.0e10 Pa.
    # Poisson's ratio: nu = lambda / (2 * (lambda + mu)) = 0.25
    # Young's Modulus: E = 2 * mu * (1 + nu) = 75 GPa = 75,000 MPa
    # Effective vertical modulus under plane strain with free horizontal boundaries:
    # E_eff = E / (1 - nu^2) = 4 * mu * (lambda + mu) / (lambda + 2 * mu) = 80 GPa = 80,000 MPa
    E_eff = 80000.0  # in MPa
    
    analytical_strains = np.array([min(strains), max(strains)])
    analytical_stresses = E_eff * analytical_strains
    plt.plot(analytical_strains, analytical_stresses, '--', color='#d62728', linewidth=1.5, label=f'Analytical (E_eff = {E_eff/1e3:.1f} GPa)')

    plt.title(r'Vertical Stress ($\sigma_{zz}$) vs. Vertical Strain ($\epsilon_{zz}$)', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel(r'Vertical Strain ($\epsilon_{zz}$)', fontsize=12, labelpad=10)
    plt.ylabel(r'Vertical Stress ($\sigma_{zz}$, MPa)', fontsize=12, labelpad=10)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(fontsize=11, loc='best')
    
    # Adjust layout and save figure
    plt.tight_layout()
    plot_path = 'stress_strain.png'
    plt.savefig(plot_path, dpi=300)
    print(f"\nPlot saved successfully to '{plot_path}'")

if __name__ == '__main__':
    main()
