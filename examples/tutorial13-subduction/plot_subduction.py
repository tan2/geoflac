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
    
    # Plot final frame
    frame = nrec
    time_myr = fl.time[frame - 1]  # time in Myr
    print(f"Plotting frame {frame} at time {time_myr:.2f} Myr...")
    
    # Read mesh coordinates
    x, z = fl.read_mesh(frame)
    X_km = x / 1e3
    Z_km = z / 1e3
    
    # Read variables
    phase = fl.read_phase(frame)
    temp = fl.read_temperature(frame)
    vx, vz = fl.read_vel(frame)
    
    # Convert velocity to cm/yr
    # 1 m/s = 100 cm/s * 3.1536e7 s/yr = 3.1536e9 cm/yr
    vx_cm_yr = vx * 3.1536e9
    vz_cm_yr = vz * 3.1536e9
    
    # Try reading magma and melt fractions
    try:
        fmagma = fl.read_fmagma(frame)
        has_magma = True
    except Exception:
        fmagma = np.zeros_like(phase, dtype=float)
        has_magma = False
        print("Note: Magma fraction file (fmagma.0) not found or could not be read.")

    try:
        fmelt = fl.read_fmelt(frame)
        has_melt = True
    except Exception:
        fmelt = np.zeros_like(phase, dtype=float)
        has_melt = False
        print("Note: Melt fraction file (fmelt.0) not found or could not be read.")

    # 18 phases defined in subduction.inp
    # Map each phase ID to a specific name and hexadecimal color
    phase_info = {
        1:  ("Basalt (Anhydrous)", "#a6cee3"),
        2:  ("Cont. Crust (Upper)", "#fdbf6f"),
        3:  ("Basalt (Oceanic)", "#1f78b4"),
        4:  ("Olivine (Mantle)", "#33a02c"),
        5:  ("Schist", "#cab2d6"),
        6:  ("Cont. Crust (Lower)", "#ff7f00"),
        7:  ("Basalt (Oceanic II)", "#b2df8a"),
        8:  ("Olivine (Mantle II)", "#6a3d9a"),
        9:  ("Serpentinite (Mantle wedge)", "#ffff99"),
        10: ("Sedimentary Rock", "#b15928"),
        11: ("Sediment", "#fb9a99"),
        12: ("Weak Cont. Crust", "#fb8072"),
        13: ("Eclogite (Metamorphic Slab)", "#e31a1c"),
        14: ("Arc Crust (Volcanics)", "#e7298a"),
        15: ("Weak Mid Crust", "#bc80bd"),
        16: ("Hydrated Mantle (Melt Source)", "#8dd3c7"),
        17: ("Metamorphic Sediment", "#bebada"),
        18: ("Dry Mantle (Deep)", "#80b1d3")
    }

    # Define color array for ListedColormap
    max_phase = 18
    color_list = ["#ffffff"] * max_phase
    for pid, (_, color) in phase_info.items():
        color_list[pid - 1] = color
    cmap_phase = mcolors.ListedColormap(color_list)

    # ----------------------------------------------------
    # Plot 1: Full Subduction Zone Profile
    # ----------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(15, 6))
    
    # Plot phase distribution
    im1 = ax1.pcolormesh(X_km, Z_km, phase, cmap=cmap_phase, vmin=0.5, vmax=max_phase+0.5, shading='flat')
    
    # Add temperature isotherms (contours)
    # temp is a nodal variable
    cs1 = ax1.contour(X_km, Z_km, temp, levels=[200, 400, 600, 800, 1000, 1200], 
                       colors='#222222', linewidths=0.8, linestyles='--')
    ax1.clabel(cs1, inline=True, fmt='%d°C', fontsize=8, colors='#222222')
    
    # Overlay velocity vectors (downsampled for clarity)
    # Mesh size is 186x65. Let's skip every 5 in X and 2 in Z
    skip_x = 5
    skip_z = 2
    ax1.quiver(X_km[::skip_z, ::skip_x], Z_km[::skip_z, ::skip_x], 
               vx_cm_yr[::skip_z, ::skip_x], vz_cm_yr[::skip_z, ::skip_x],
               color='black', scale=15.0, width=0.001, alpha=0.6,
               label='Velocity field (cm/yr)')

    # Add legend for key phases
    # Include only the phases that are active and important in subduction
    from matplotlib.patches import Patch
    key_phases = [3, 4, 9, 13, 14, 16]
    legend_elements = [
        Patch(facecolor=phase_info[pid][1], edgecolor='gray', label=f"{phase_info[pid][0]} (Phase {pid})")
        for pid in key_phases
    ]
    ax1.legend(handles=legend_elements, loc='lower left', framealpha=0.9, fontsize=9)
    
    ax1.set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Depth (km)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Ocean-Ocean Subduction: Phase Distribution & Geotherm at {time_myr:.2f} Myr', 
                  fontsize=14, fontweight='bold', pad=12)
    ax1.grid(True, linestyle=':', alpha=0.3, color='gray')
    ax1.set_xlim(0, 960)
    ax1.set_ylim(-300, 10)
    
    plt.tight_layout()
    plt.savefig('subduction_full_zone.png', dpi=300)
    plt.savefig('plots/subduction_full_zone.png', dpi=300)
    plt.close()
    
    # ----------------------------------------------------
    # Plot 2: Zoomed-in Mantle Wedge and Arc Melting Zone
    # ----------------------------------------------------
    fig2, axes2 = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Set zoom limits (focus on the mantle wedge / subduction interface)
    x_zoom_min, x_zoom_max = 300, 700
    z_zoom_min, z_zoom_max = -150, 10

    # Panel 2A: Zoomed Phase & Velocity Field
    ax_a = axes2[0]
    ax_a.pcolormesh(X_km, Z_km, phase, cmap=cmap_phase, vmin=0.5, vmax=max_phase+0.5, shading='flat')
    cs_a = ax_a.contour(X_km, Z_km, temp, levels=[400, 600, 800, 1000, 1200], 
                         colors='#222222', linewidths=1.0, linestyles='--')
    ax_a.clabel(cs_a, inline=True, fmt='%d°C', fontsize=8, colors='#222222')
    
    # Downsample velocity quiver for zoom
    skip_x_z = 2
    skip_z_z = 1
    ax_a.quiver(X_km[::skip_z_z, ::skip_x_z], Z_km[::skip_z_z, ::skip_x_z], 
               vx_cm_yr[::skip_z_z, ::skip_x_z], vz_cm_yr[::skip_z_z, ::skip_x_z],
               color='black', scale=10.0, width=0.0015, headwidth=4, headlength=5)
    
    ax_a.set_ylabel('Depth (km)', fontsize=11, fontweight='bold')
    ax_a.set_title('Mantle Wedge: Phase Transformations & Dehydration Flow (Zoom)', 
                  fontsize=13, fontweight='bold', pad=10)
    ax_a.grid(True, linestyle=':', alpha=0.3, color='gray')
    ax_a.set_xlim(x_zoom_min, x_zoom_max)
    ax_a.set_ylim(z_zoom_min, z_zoom_max)
    
    # Panel 2B: Magma & Melt Fraction in the Wedge
    ax_b = axes2[1]
    
    # Plot temperature isotherms as background
    ax_b.contourf(X_km, Z_km, temp, levels=15, cmap='Oranges', alpha=0.1)
    
    # Plot fmagma or fmelt if available
    if has_magma or has_melt:
        plot_var = fmagma if has_magma else fmelt
        var_label = 'Magma Fraction' if has_magma else 'Melt Fraction'
        cmap_magma = plt.colormaps.get_cmap('Reds')
        cmap_magma.set_under(alpha=0.0) # transparent for zero magma
        
        im_magma = ax_b.pcolormesh(X_km, Z_km, plot_var, cmap=cmap_magma, vmin=0.001, vmax=0.10, shading='flat')
        cbar_magma = fig2.colorbar(im_magma, ax=ax_b, orientation='vertical', pad=0.02, shrink=0.8)
        cbar_magma.set_label(f'{var_label} (Vol %)', fontsize=10, fontweight='bold')
    else:
        # Fallback representation of melting zone based on temperature & hydrated phase
        # Just plot a colored zone where phase == 16 (hydrated mantle)
        hydrated_mask = (phase == 16).astype(float)
        cmap_fallback = mcolors.ListedColormap(['none', '#8dd3c7'])
        im_fallback = ax_b.pcolormesh(X_km, Z_km, hydrated_mask, cmap=cmap_fallback, shading='flat', vmin=0.5, vmax=1.5)
        
        # Legend indicating hydrated mantle
        legend_elements_b = [Patch(facecolor='#8dd3c7', edgecolor='gray', label='Hydrated Mantle wedge (Potential Melt Source)')]
        ax_b.legend(handles=legend_elements_b, loc='upper left')

    # Add slab slab contours for reference (basalt #3 and eclogite #13)
    slab_mask = ((phase == 3) | (phase == 13)).astype(float)
    ax_b.contour(X_km[:-1, :-1] + 0.5 * (X_km[1:, 1:] - X_km[:-1, :-1]), 
                 Z_km[:-1, :-1] + 0.5 * (Z_km[1:, 1:] - Z_km[:-1, :-1]), 
                 slab_mask, levels=[0.5], colors='#1f78b4', linewidths=1.2)
    
    # Overlay geotherm contours
    cs_b = ax_b.contour(X_km, Z_km, temp, levels=[600, 800, 1000, 1200], 
                         colors='#444444', linewidths=0.8, linestyles='--')
    ax_b.clabel(cs_b, inline=True, fmt='%d°C', fontsize=8, colors='#444444')
    
    ax_b.set_xlabel('Distance (km)', fontsize=11, fontweight='bold')
    ax_b.set_ylabel('Depth (km)', fontsize=11, fontweight='bold')
    ax_b.set_title('Magma Generation & Diking Conduit Zone', fontsize=13, fontweight='bold', pad=10)
    ax_b.grid(True, linestyle=':', alpha=0.3, color='gray')
    ax_b.set_xlim(x_zoom_min, x_zoom_max)
    ax_b.set_ylim(z_zoom_min, z_zoom_max)
    
    plt.tight_layout()
    plt.savefig('subduction_mantle_wedge.png', dpi=300)
    plt.savefig('plots/subduction_mantle_wedge.png', dpi=300)
    plt.close()
    
    print("Saved 'subduction_full_zone.png' and 'subduction_mantle_wedge.png'")
    
    # ----------------------------------------------------
    # Plot 3: Subduction Evolution Sequence
    # ----------------------------------------------------
    if nrec >= 3:
        fig_evo, axes_evo = plt.subplots(3, 1, figsize=(14, 11), sharex=True, sharey=True)
        frames_to_plot = [1, nrec // 2 + 1, nrec]
        
        for idx, f_idx in enumerate(frames_to_plot):
            ax = axes_evo[idx]
            x_f, z_f = fl.read_mesh(f_idx)
            phase_f = fl.read_phase(f_idx)
            t_f = fl.time[f_idx - 1]
            
            ax.pcolormesh(x_f/1e3, z_f/1e3, phase_f, cmap=cmap_phase, vmin=0.5, vmax=max_phase+0.5, shading='flat')
            ax.set_title(f'Subduction Evolution at t = {t_f:.2f} Myr', fontsize=12, fontweight='bold')
            ax.set_ylabel('Depth (km)', fontsize=10)
            ax.grid(True, linestyle=':', alpha=0.3)
            ax.set_ylim(-300, 10)
            ax.set_xlim(0, 960)
            
        axes_evo[-1].set_xlabel('Distance (km)', fontsize=11, fontweight='bold')
        plt.savefig('plots/subduction_evolution.png', dpi=300)
        print("Saved 'plots/subduction_evolution.png'")

if __name__ == '__main__':
    main()
