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
    nrec = fl.nrec
    print(f"Number of available records: {nrec}")
    if nrec == 0:
        print("Error: No simulation data found. Please run the simulation first.")
        sys.exit(1)

    # Create a directory for plots if it doesn't exist
    os.makedirs('images', exist_ok=True)
    
    # Use a warm/hot colormap for strain localization
    cmap_aps = plt.colormaps.get_cmap('YlOrRd')
    # Set background for zero strain
    cmap_aps.set_under('#fcfcfc')
    
    # We will select 3 frames for the evolution plot:
    # 1. Initial (1)
    # 2. Middle (nrec // 2 + 1)
    # 3. Final (nrec)
    frames_to_plot = [1, nrec // 2 + 1, nrec]
    
    fig, axes = plt.subplots(3, 1, figsize=(11, 12), sharex=True, sharey=True)
    
    # Downsampling parameters for velocity vectors (nx=86, nz=15)
    skip_x = 3
    skip_z = 2
    
    im = None
    for idx, f_idx in enumerate(frames_to_plot):
        ax = axes[idx]
        x_f, z_f = fl.read_mesh(f_idx)
        aps_f = fl.read_aps(f_idx)
        vx_f, vz_f = fl.read_vel(f_idx)
        t_f = fl.time[f_idx - 1]
        
        # Plot plastic strain
        im = ax.pcolormesh(x_f, z_f, aps_f, cmap=cmap_aps, vmin=0.01, vmax=4.0, shading='flat')
        
        # Plot topography line
        ax.plot(x_f[:, 0], z_f[:, 0], color='black', linewidth=1.5)
        
        # Plot velocity vectors (downsampled)
        q = ax.quiver(x_f[::skip_z, ::skip_x], z_f[::skip_z, ::skip_x], 
                      vx_f[::skip_z, ::skip_x], vz_f[::skip_z, ::skip_x],
                      color='#1d3557', scale=15.0, width=0.0015, headwidth=4, headlength=5)
        
        # Add velocity scale key at the bottom right corner of each subplot
        ax.quiverkey(q, X=0.90, Y=0.08, U=1.0, label='1 cm/yr', labelpos='E', 
                     coordinates='axes', fontproperties={'weight': 'bold', 'size': 9})
        
        ax.set_title(f'Accumulated Plastic Strain & Velocity Field at t = {t_f:.2f} Myr', 
                     fontsize=12, fontweight='bold', pad=8)
        ax.set_ylabel('Depth (km)', fontsize=11, fontweight='bold')
        ax.set_ylim(-13, 5)
        ax.grid(True, linestyle=':', alpha=0.5, color='gray')
        ax.set_aspect('equal')
        
    axes[-1].set_xlabel('Distance (km)', fontsize=12, fontweight='bold')
    
    # Reduce vertical space between stacked subplots
    fig.subplots_adjust(hspace=0.02)
    
    # Add common colorbar
    cbar = fig.colorbar(im, ax=axes.tolist(), orientation='vertical', pad=0.02, shrink=0.6)
    cbar.set_label('Accumulated Plastic Strain', fontsize=11, fontweight='bold')
    
    # Save the evolution figure
    plot_path = 'images/evolution_core_complex.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved '{plot_path}' successfully.")

if __name__ == '__main__':
    main()
