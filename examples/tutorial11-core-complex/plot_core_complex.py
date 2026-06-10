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
    os.makedirs('images', exist_ok=True)
    
    # We will generate a premium plot for the final state
    final_frame = nrec
    time_myr = fl.time[final_frame - 1]  # time in Myr
    print(f"Plotting final frame {final_frame} at time {time_myr:.2f} Myr...")
    
    # Read mesh coordinates
    x, z = fl.read_mesh(final_frame)
    X_km = x / 1e3
    Z_km = z / 1e3
    
    # Read variables
    aps = fl.read_aps(final_frame)
    phase = fl.read_phase(final_frame)
    temp = fl.read_temperature(final_frame)
    vx, vz = fl.read_vel(final_frame)
    
    # Convert velocity to cm/yr
    # 1 m/s = 100 cm/s * 3.1536e7 s/yr = 3.1536e9 cm/yr
    vx_cm_yr = vx * 3.1536e9
    vz_cm_yr = vz * 3.1536e9
    
    # Create the final state premium plot
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Subplot 1: Accumulated Plastic Strain (aps) and Detachment Faults
    ax1 = axes[0]
    # Use a warm/hot colormap for strain localization
    cmap_aps = plt.colormaps.get_cmap('YlOrRd')
    # Set background for zero strain
    cmap_aps.set_under('#fcfcfc')
    
    # Plot aps
    im1 = ax1.pcolormesh(X_km, Z_km, aps, cmap=cmap_aps, vmin=0.01, vmax=4.0, shading='flat')
    cbar1 = fig.colorbar(im1, ax=ax1, orientation='vertical', pad=0.02, shrink=0.8)
    cbar1.set_label('Accumulated Plastic Strain', fontsize=11, fontweight='bold')
    
    # Overlay velocity arrows (downsample for clarity)
    # nx is 51, nz is 16. Let's skip every 2 elements in x and 1 in z
    skip_x = 2
    skip_z = 1
    ax1.quiver(X_km[::skip_z, ::skip_x], Z_km[::skip_z, ::skip_x], 
               vx_cm_yr[::skip_z, ::skip_x], vz_cm_yr[::skip_z, ::skip_x],
               color='#1d3557', scale=15.0, width=0.0015, headwidth=4, headlength=5,
               label='Velocity field (cm/yr)')
    ax1.legend(loc='upper right', framealpha=0.9)
    
    ax1.set_ylabel('Depth (km)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Metamorphic Core Complex: Strain Localization & Exhumation at {time_myr:.2f} Myr', 
                  fontsize=14, fontweight='bold', pad=12)
    ax1.grid(True, linestyle=':', alpha=0.5, color='gray')
    ax1.set_ylim(-30, 2)
    
    # Subplot 2: Lithological Phases & Temperature Structure (Isotherms)
    ax2 = axes[1]
    # Custom discrete colormap for phases
    # Phase 1: Upper Crust (light beige/brown), Phase 2: Lower Crust (light green/grey)
    colors_phase = ['#e5c494', '#8da0cb']
    cmap_phase = mcolors.ListedColormap(colors_phase)
    
    im2 = ax2.pcolormesh(X_km, Z_km, phase, cmap=cmap_phase, shading='flat', vmin=0.5, vmax=2.5)
    
    # Legend for phases
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e5c494', edgecolor='gray', label='Upper Crust (Phase 1)'),
        Patch(facecolor='#8da0cb', edgecolor='gray', label='Lower Crust (Phase 2)')
    ]
    ax2.legend(handles=legend_elements, loc='upper right', framealpha=0.9)
    
    # Overlay temperature isotherms
    # Compute temperature node coordinates (temp is on nodes, shape (nx, nz))
    levels = [200, 300, 400, 500, 600, 700]
    cs = ax2.contour(X_km, Z_km, temp, levels=levels, colors='#e63946', linewidths=1.5, linestyles='--')
    ax2.clabel(cs, inline=True, fmt='%d°C', fontsize=9, colors='#e63946')
    
    ax2.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Depth (km)', fontsize=12, fontweight='bold')
    ax2.set_title(f'Crustal Structure & Thermal Geotherm (Isotherms) at {time_myr:.2f} Myr', 
                  fontsize=14, fontweight='bold', pad=12)
    ax2.grid(True, linestyle=':', alpha=0.5, color='gray')
    ax2.set_ylim(-30, 2)
    
    plt.tight_layout()
    plt.savefig('images/core_complex.png', dpi=300)
    print("Saved 'images/core_complex.png'")
    
    # Let's also create an evolution plot if we have multiple frames
    if nrec >= 3:
        fig_evo, axes_evo = plt.subplots(3, 1, figsize=(12, 12), sharex=True, sharey=True)
        # Select 3 frames: initial (1), middle (nrec // 2 + 1), final (nrec)
        frames_to_plot = [1, nrec // 2 + 1, nrec]
        
        for idx, f_idx in enumerate(frames_to_plot):
            ax = axes_evo[idx]
            x_f, z_f = fl.read_mesh(f_idx)
            aps_f = fl.read_aps(f_idx)
            t_f = fl.time[f_idx - 1]
            
            im = ax.pcolormesh(x_f/1e3, z_f/1e3, aps_f, cmap=cmap_aps, vmin=0.01, vmax=4.0, shading='flat')
            ax.set_title(f'Accumulated Plastic Strain at t = {t_f:.2f} Myr', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (km)', fontsize=10)
            ax.set_ylim(-30, 2)
            ax.grid(True, linestyle=':', alpha=0.5)
            
        axes_evo[-1].set_xlabel('Distance (km)', fontsize=11, fontweight='bold')
        fig_evo.colorbar(im, ax=axes_evo.tolist(), orientation='vertical', pad=0.02, shrink=0.6, label='Plastic Strain')
        plt.savefig('images/evolution_core_complex.png', dpi=300)
        print("Saved 'images/evolution_core_complex.png'")
        
if __name__ == '__main__':
    main()
