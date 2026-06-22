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
    import matplotlib.colors as mcolors
except ImportError:
    print("Error: matplotlib is required to plot. Please install it.")
    sys.exit(1)

def main():
    # Instantiate Flac reader
    fl = flac.Flac()
    nrec = fl.nrec
    print(f"Number of available records: {nrec}")
    if nrec == 0:
        print("Error: No simulation data found. Please run the simulation first.")
        sys.exit(1)

    # Create a directory for plots if it doesn't exist
    os.makedirs('plots', exist_ok=True)
    
    # We will generate a premium plot for the final state
    final_frame = nrec
    time_myr = fl.time[final_frame - 1]  # time in Myr
    print(f"Plotting final frame {final_frame} at time {time_myr:.2f} Myr...")
    
    # Read mesh coordinates
    x, z = fl.read_mesh(final_frame)
    X_km = x
    Z_km = z
    
    # Read variables
    aps = fl.read_aps(final_frame)
    phase = fl.read_phase(final_frame)
    vx, vz = fl.read_vel(final_frame)
    temp = fl.read_temperature(final_frame)
    
    vx_cm_yr = vx
    vz_cm_yr = vz
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Subplot 1: Accumulated Plastic Strain (aps) and Temperature Isotherms
    ax1 = axes[0]
    cmap_aps = plt.colormaps.get_cmap('YlOrRd')
    cmap_aps.set_under('#fcfcfc')
    
    # Plot aps
    im1 = ax1.pcolormesh(X_km, Z_km, aps, cmap=cmap_aps, vmin=0.01, vmax=2.0, shading='flat')
    ax1.plot(X_km[:, 0], Z_km[:, 0], color='black', linewidth=1.5)
    
    # Plot temperature isotherms (contour lines)
    cs = ax1.contour(X_km, Z_km, temp, levels=[100, 200, 300, 400, 500, 600, 700], colors='#333333', linewidths=0.8, alpha=0.8)
    ax1.clabel(cs, inline=True, fmt='%d°C', fontsize=8)
    
    cbar1 = fig.colorbar(im1, ax=ax1, orientation='vertical', pad=0.02, shrink=0.8)
    cbar1.set_label('Accumulated Plastic Strain', fontsize=11, fontweight='bold')
    
    # Overlay velocity arrows (downsample for clarity)
    skip_x = 4
    skip_z = 1
    ax1.quiver(X_km[::skip_z, ::skip_x], Z_km[::skip_z, ::skip_x], 
               vx_cm_yr[::skip_z, ::skip_x], vz_cm_yr[::skip_z, ::skip_x],
               color='#1d3557', scale=20.0, width=0.0015, headwidth=4, headlength=5,
               label='Velocity field (cm/yr)')
    ax1.legend(loc='upper left', framealpha=0.9)
    
    ax1.set_ylabel('Depth (km)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Thermochron-Wedge: Shear Localization & Isotherms at {time_myr:.2f} Myr', 
                  fontsize=14, fontweight='bold', pad=12)
    ax1.grid(True, linestyle=':', alpha=0.5, color='gray')
    ax1.set_ylim(-30, 5)
    ax1.set_aspect('equal')
    
    # Subplot 2: Lithological Phases
    ax2 = axes[1]
    # We have up to 19 phases. We will use a colormap to distinguish them
    cmap_phase = plt.colormaps.get_cmap('tab20')
    
    im2 = ax2.pcolormesh(X_km, Z_km, phase, cmap=cmap_phase, shading='flat', vmin=0.5, vmax=19.5)
    ax2.plot(X_km[:, 0], Z_km[:, 0], color='black', linewidth=1.5)
    
    # Overlay grid lines to show mesh deformation
    for i in range(0, x.shape[0], 4):
        ax2.plot(x[i, :], z[i, :], color='black', alpha=0.15, linewidth=0.5)
    for j in range(0, x.shape[1], 2):
        ax2.plot(x[:, j], z[:, j], color='black', alpha=0.15, linewidth=0.5)
        
    ax2.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Depth (km)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Deformed Lithological Phases & Grid Mesh at {time_myr:.2f} Myr', 
                  fontsize=14, fontweight='bold', pad=12)
    ax2.grid(True, linestyle=':', alpha=0.5, color='gray')
    ax2.set_ylim(-30, 5)
    ax2.set_aspect('equal')
    
    # Add a general colorbar for phases to identify different geological layers
    cbar2 = fig.colorbar(im2, ax=ax2, orientation='vertical', pad=0.02, shrink=0.8)
    cbar2.set_label('Lithology/Phase ID', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('thermochron_wedge.png', dpi=300)
    plt.savefig('plots/final_state_thermochron_wedge.png', dpi=300)
    print("Saved 'thermochron_wedge.png' and 'plots/final_state_thermochron_wedge.png'")
    
    # Create an evolution plot
    if nrec >= 3:
        fig_evo, axes_evo = plt.subplots(3, 1, figsize=(12, 12), sharex=True, sharey=True)
        frames_to_plot = [1, nrec // 2 + 1, nrec]
        
        for idx, f_idx in enumerate(frames_to_plot):
            ax = axes_evo[idx]
            x_f, z_f = fl.read_mesh(f_idx)
            aps_f = fl.read_aps(f_idx)
            t_f = fl.time[f_idx - 1]
            
            im = ax.pcolormesh(x_f, z_f, aps_f, cmap=cmap_aps, vmin=0.01, vmax=2.0, shading='flat')
            ax.plot(x_f[:, 0], z_f[:, 0], color='black', linewidth=1.5)
            ax.set_title(f'Accumulated Plastic Strain at t = {t_f:.2f} Myr', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (km)', fontsize=10)
            ax.set_ylim(-30, 5)
            ax.grid(True, linestyle=':', alpha=0.5)
            ax.set_aspect('equal')
            
        axes_evo[-1].set_xlabel('Distance (km)', fontsize=11, fontweight='bold')
        fig_evo.colorbar(im, ax=axes_evo.tolist(), orientation='vertical', pad=0.02, shrink=0.6, label='Plastic Strain')
        plt.savefig('plots/evolution_thermochron_wedge.png', dpi=300)
        print("Saved 'plots/evolution_thermochron_wedge.png'")

if __name__ == '__main__':
    main()
