#!/usr/bin/python3
'''
Read flac density and coordinate data and return gravity and topography at regular spacing.
'''

import sys
import argparse
import numpy as np
import flac

gn = 6.6726e-11
ocean_density = 1030.
sediment_density = 2200. # from flac_gravity3.py inside compute_gravity

def find_trench_index(zz):
    imax = zz.argmax()
    i = zz[:imax].argmin()
    return i

def find_peaks(zz):
    nx = zz.shape[0]
    itrench = find_trench_index(zz)
    peaks = []
    for i in range(itrench+1, nx-1):
        left_dz = zz[i] - zz[i-1]
        right_dz = zz[i+1] - zz[i]
        if left_dz * right_dz <= 0 and (left_dz > 0 or right_dz < 0):
            peaks.append(i)
    peaks.append(nx-1)
    return peaks

def compute_gravity(frame, fill_basin=False, erode_forearc=False, ref_grav="mean"):
    fl = flac.Flac()
    x, z = fl.read_mesh(frame)  # in km
    x *= 1e3
    z *= 1e3

    xx = x[:,0]
    zz = z[:,0].copy() # copy because we might modify zz when filling basin
    xmin = x[0,0]
    xmax = x[-1,0]

    cx, cz = flac.elem_coord(x, z)

    xdiag1 = x[0:-1, 0:-1] - x[1:, 1:]
    zdiag1 = z[0:-1, 0:-1] - z[1:, 1:]
    xdiag2 = x[1:, 0:-1] - x[0:-1, 1:]
    zdiag2 = z[1:, 0:-1] - z[0:-1, 1:]
    area = 0.5 * np.abs(xdiag1*zdiag2 - xdiag2*zdiag1)

    rho = fl.read_density(frame)  # in kg/m^3

    if erode_forearc:
        # anything above sea level is removed
        rho[cz > 0] = 0

    mass = rho * area

    px = np.linspace(xmin, xmax, num=5*fl.nx)
    if erode_forearc:
        pz_height = max(0.0, np.max(zz)) + 4e3
    else:
        pz_height = np.max(zz) + 4e3
    pz = np.ones(px.shape) * pz_height

    topo = np.interp(px, x[:,0], z[:,0])

    grav = np.empty(px.shape)
    for i in range(len(grav)):
        dx = px[i] - cx
        dz = pz[i] - cz
        dist2 = (dx**2 + dz**2)
        grav[i] = 2 * gn * np.sum(mass * dz / dist2)

    topomod = None
    if fill_basin:
        peaks = find_peaks(zz)
        itrench = find_trench_index(zz)
        basin_depth = -2000
        sed_thickness = np.zeros(fl.nx)

        for ii in range(len(peaks)-1):
            fill_height = min((basin_depth, zz[peaks[ii]], zz[peaks[ii+1]]))
            for i in range(peaks[ii], peaks[ii+1]):
                if zz[i] < fill_height:
                    sed_thickness[i] = fill_height - zz[i]
                    zz[i] = fill_height

        for i in range(itrench, fl.nx-1):
            sedz = 0.5 * (sed_thickness[i] + sed_thickness[i+1])
            if sedz > 0:
                midz = 0.5 * (zz[i] + zz[i+1])
                midx = 0.5 * (xx[i] + xx[i+1])
                m = (xx[i+1] - xx[i]) * sedz * sediment_density
                dx = px - midx
                dz = pz - midz
                dist2 = (dx**2 + dz**2)
                grav += 2 * gn * m * dz / dist2

        # modify topo for output
        topomod = topo.copy()
        if erode_forearc:
            topomod[topo > 0] = 0
        
        for ii in range(len(peaks)-1):
            fill_height = min((basin_depth, topo[peaks[ii]], topo[peaks[ii+1]]))
            for i in range(peaks[ii], peaks[ii+1]):
                if topo[i] < fill_height:
                    topomod[i] = fill_height

    # ocean
    for i in range(fl.nx-1):
        midz = 0.5 * (zz[i] + zz[i+1])
        if midz < 0:
            midx = 0.5 * (xx[i] + xx[i+1])
            m = (xx[i+1] - xx[i]) * -midz * ocean_density
            dx = px - midx
            dz = pz - midz
            dist2 = (dx**2 + dz**2)
            grav += 2 * gn * m * dz / dist2

    # left boundary
    for i in range(fl.nz-1):
        sigma = rho[0,i] * (z[0,i] - z[0,i+1])
        dx = px - xmin
        dz = pz - 0.5 * (z[0,i] + z[0,i+1])
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if zz[0] < 0:
        sigma = ocean_density * -zz[0]
        dx = px - xmin
        dz = pz - 0.5 * zz[0]
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    # right boundary
    for i in range(fl.nz-1):
        sigma = rho[-1,i] * (z[-1,i] - z[-1,i+1])
        dx = xmax - px
        dz = pz - 0.5 * (z[-1,i] + z[-1,i+1])
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if zz[-1] < 0:
        sigma = ocean_density * -zz[-1]
        dx = xmax - px
        dz = pz - 0.5 * zz[-1]
        angle = np.arctan2(dx, dz)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if ref_grav == "right":
        grav -= grav[-1]
    else:
        grav -= np.mean(grav)

    return px, topo, topomod, grav

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute gravity and topography from FLAC output.')
    parser.add_argument('frame', type=int, help='Output frame number')
    parser.add_argument('-s', '--sediment', action='store_true', help='Fill basins with sediment')
    parser.add_argument('-e', '--erode', action='store_true', help='Erode forearc height to sea level')
    parser.add_argument('-r', '--reference', choices=['mean', 'right'], default='mean',
                        help='Reference gravity setting: subtract mean (mean) or subtract far-right (right)')
    
    args = parser.parse_args()

    px, topo, topomod, gravity = compute_gravity(args.frame, fill_basin=args.sediment, erode_forearc=args.erode, ref_grav=args.reference)
    
    if args.sediment:
        flac.printing(px, topo, topomod, gravity)
    else:
        flac.printing(px, topo, gravity)
