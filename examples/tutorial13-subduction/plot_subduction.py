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

    # Create a directory for images if it doesn't exist
    os.makedirs('images', exist_ok=True)
    
    # 18 phases defined in subduction.inp
    # Map each phase ID to a specific name and hexadecimal color
    phase_info = {
        1:  ("Sediment", "#fdbf6f"),
        2:  ("Continental Crust", "#fdbf6f"),
        3:  ("Oceanic Crust (Basalt)", "#1f78b4"),
        4:  ("Mantle Olivine", "#7fc97f"),
        5:  ("Schist", "#cab2d6"),
        6:  ("Lower Cont. Crust", "#ff7f00"),
        7:  ("Basaltic Crust II", "#b2df8a"),
        8:  ("Mantle Olivine II", "#a6d854"),
        9:  ("Serpentinite", "#ffff99"),
        10: ("Sedimentary Rock", "#b15928"),
        11: ("Sediment II", "#fb9a99"),
        12: ("Weak Crust", "#fb8072"),
        13: ("Eclogite (Slab)", "#e31a1c"),
        14: ("Arc Crust (Volcanics)", "#e7298a"),
        15: ("Weak Mid Crust", "#bc80bd"),
        16: ("Hydrated Mantle (Water Source)", "#8dd3c7"),
        17: ("Metamorphic Sediment", "#bebada"),
        18: ("Dry Mantle", "#33a02c")
    }

    # Define active phases and a compact ListedColormap containing only their colors
    active_phases = [1, 3, 4, 5, 9, 13, 16]
    cmap_active = mcolors.ListedColormap([phase_info[pid][1] for pid in active_phases])

    # ----------------------------------------------------
    # Plot: Subduction Evolution Sequence
    # ----------------------------------------------------
    if nrec >= 3:
        fig_evo, axes_evo = plt.subplots(3, 1, figsize=(14, 11), sharex=True, sharey=True)
        frames_to_plot = [1, nrec // 2 + 1, nrec]
        
        im = None
        for idx, f_idx in enumerate(frames_to_plot):
            ax = axes_evo[idx]
            x_f, z_f = fl.read_mesh(f_idx)
            phase_f = fl.read_phase(f_idx)
            temp_f = fl.read_temperature(f_idx)
            t_f = fl.time[f_idx - 1]
            
            # Map original phase IDs to 0..6 for compact mapping
            phase_mapped = np.zeros_like(phase_f, dtype=int)
            for p_idx, pid in enumerate(active_phases):
                phase_mapped[phase_f == pid] = p_idx
                
            im = ax.pcolormesh(x_f, z_f, phase_mapped, cmap=cmap_active, vmin=-0.5, vmax=len(active_phases)-0.5, shading='flat')
            
            # Overlay temperature contours
            cs = ax.contour(x_f, z_f, temp_f, levels=[200, 400, 600, 800, 1000, 1200], 
                             colors='#333333', linewidths=0.8, linestyles='--')
            ax.clabel(cs, inline=True, fmt='%d°C', fontsize=8, colors='#333333')
            
            ax.set_title(f'Subduction Evolution at t = {t_f:.2f} Myr', fontsize=12, fontweight='bold', pad=5)
            ax.set_ylabel('Depth (km)', fontsize=11, fontweight='bold')
            ax.grid(True, linestyle=':', alpha=0.3, color='gray')
            ax.set_ylim(-300, 10)
            ax.set_xlim(0, 960)
            ax.set_aspect('equal')
            
        axes_evo[-1].set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
        
        # Set vertical space between the stacked subplots
        fig_evo.subplots_adjust(hspace=0.20)
        
        # Add a single colorbar for all panels
        cbar = fig_evo.colorbar(im, ax=axes_evo.tolist(), orientation='vertical', pad=0.02, shrink=0.6)
        cbar.set_label('Lithological Phases', fontsize=11, fontweight='bold')
        
        # Use rock types instead of numbers for the colorbar labels
        cbar.set_ticks(range(len(active_phases)))
        cbar.set_ticklabels([phase_info[pid][0] for pid in active_phases])
        
        plot_path = 'images/subduction_evolution.png'
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved '{plot_path}' successfully.")
    else:
        print("Note: Evolution plot requires at least 3 records of output.")

if __name__ == '__main__':
    main()
