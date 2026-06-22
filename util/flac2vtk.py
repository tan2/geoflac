#!/usr/bin/python3

'''Convert the binary output of flac to VTK (vts) files.
'''
from __future__ import print_function
import sys, os
import zlib, base64, glob
import numpy as np
import flac

ndims = 2
nstress = 4

# whether write VTK files in compressed binary (base64 encoded) or uncompressed plain ascii
output_in_binary = True

max_thermochron_temp = 620.0

def parse_chron_file(filepath):
    systems = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    current_system = None
    rates = []
    temps = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith(';'):
            continue
        if line.startswith('0 (#'):
            if current_system is not None:
                systems.append((np.array(rates), np.array(temps)))
            current_system = line
            rates = []
            temps = []
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rate = float(parts[0])
                temp = float(parts[1])
                rates.append(rate)
                temps.append(temp)
            except ValueError:
                pass
    if current_system is not None:
        systems.append((np.array(rates), np.array(temps)))
    return systems

def main(path, start=1, end=-1, thermochron=False):

    # changing directory
    os.chdir(path)

    fl = flac.Flac()

    nx, nz = fl.nx, fl.nz
    nex = nx - 1
    nez = nz - 1
    if end == -1:
        end = fl.nrec

    if start == -1:
        vtslist = sorted(glob.glob('flac.*.vts'))
        lastframe = int(vtslist[-1][5:-4]) if vtslist else 0
        start = lastframe + 1

    if thermochron:
        try:
            from scipy.interpolate import griddata
            import pandas as pd
        except ImportError:
            print("Error: scipy and pandas are required for thermochronology interpolation. Please install them.")
            sys.exit(1)

        script_dir = os.path.dirname(os.path.abspath(__file__))
        chron_file = os.path.join(script_dir, 'thermo_chron.dat')
        if not os.path.exists(chron_file):
            print(f"Error: thermochron reference file not found at {chron_file}")
            sys.exit(1)

        if not os.path.exists('coolingrate.0'):
            print("Error: coolingrate.0 not found in the current directory.")
            print("Please make sure you have run the model with the updated outflac.f90 to generate it.")
            sys.exit(1)

        print(f"Reading thermochron reference from {chron_file}...")
        systems = parse_chron_file(chron_file)
        nchron = len(systems)
        thermochron_system_names = ['ZFT', 'ZHe', 'AFT', 'AHE', 'orthoclase', 'biotite', 'muscovite', 'hornblende'][:nchron]
        print(f"Loaded {nchron} thermochronometer systems: {thermochron_system_names}")

        # Read max_markers from _markers.0 to pre-allocate correctly
        if not os.path.exists('_markers.0'):
            print("Error: _markers.0 not found in the current directory.")
            sys.exit(1)
        tmp_markers = np.loadtxt('_markers.0')
        if tmp_markers.ndim == 1:
            tmp_markers = tmp_markers.reshape((1, -1))
        max_markers = int(np.max(tmp_markers[:, 2]))

        # Initialize tracking arrays on markers of size (max_markers + 1)
        unreset_time = -100.0
        chron_time = [np.full(max_markers + 1, unreset_time, dtype=np.float32) for _ in range(nchron)]
        chron_if = [np.ones(max_markers + 1, dtype=np.int32) for _ in range(nchron)]

    loop_start = 1 if thermochron else start
    for i in range(loop_start, end+1):
        if thermochron:
            t_frame = fl.time[i - 1]
            # Read temperature and cooling rate at nodes
            x, z = fl.read_mesh(i)
            T_nodes = fl.read_temperature(i)
            CR_nodes = fl.read_coolingrate(i)

            # Read markers
            x_markers, z_markers, age_markers, phase_markers, ID_markers, a1_markers, a2_markers, ntriag_markers = fl.read_markers(i)
            nmarkers = len(ID_markers)

            if nmarkers > 0:

                # 1. Compute temperature and cooling rate on each marker
                t_temps = flac.marker_interpolate_node(ntriag_markers, a1_markers, a2_markers, nz, T_nodes)
                t_rates = flac.marker_interpolate_node(ntriag_markers, a1_markers, a2_markers, nz, CR_nodes)

                # Filter out markers with temperature > max_thermochron_temp
                valid = (t_temps <= max_thermochron_temp)
                x_markers = x_markers[valid]
                z_markers = z_markers[valid]
                ID_markers = ID_markers[valid]
                ntriag_markers = ntriag_markers[valid]
                t_temps = t_temps[valid]
                t_rates = t_rates[valid]
                nmarkers = len(ID_markers)

                # 2. Update closure states and closure times
                for c in range(nchron):
                    rates_ref, temps_ref = systems[c]
                    cooling_rate_val = np.maximum(1e-10, -t_rates)
                    t_closure = np.interp(cooling_rate_val, rates_ref, temps_ref)

                    marker_ids = ID_markers
                    if_closed = chron_if[c][marker_ids]
                    c_time = chron_time[c][marker_ids]

                    open_mask = (t_temps > t_closure)
                    close_mask = (~open_mask) & (if_closed == 0)

                    if_closed[open_mask] = 0
                    c_time[open_mask] = t_frame

                    if_closed[close_mask] = 1
                    c_time[close_mask] = t_frame

                    chron_if[c][marker_ids] = if_closed
                    chron_time[c][marker_ids] = c_time

            if i >= start:
                marker_ages = [np.zeros(nmarkers, dtype=np.float32) for _ in range(nchron)]
                if nmarkers > 0:
                    for c in range(nchron):
                        c_time = chron_time[c][ID_markers]
                        ages = (t_frame - c_time).astype(np.float32)
                        unreset_mask = (c_time == unreset_time)
                        ages[unreset_mask] = np.nan
                        marker_ages[c] = ages

                # 4. Interpolate ages from markers to nodes
                nodal_ages = [np.zeros((nx, nz)) for _ in range(nchron)]
                # Compute element index for each marker
                ntriag_clean = np.maximum(1, ntriag_markers)
                k_m = (ntriag_clean - 1) % 2 + 1
                j_el = np.clip(((ntriag_clean - k_m) // 2) % (nz - 1), 0, nz - 2)
                i_el = np.clip(((ntriag_clean - k_m) // 2) // (nz - 1), 0, nx - 2)

                if nmarkers > 0:
                    points = np.column_stack((x_markers, z_markers))
                    grid_x, grid_z = x, z
                    for c in range(nchron):
                        ages = marker_ages[c]
                        valid = ~np.isnan(ages)
                        if np.any(valid):
                            points_valid = points[valid]
                            ages_valid = ages[valid]
                            nodal_age = griddata(points_valid, ages_valid, (grid_x, grid_z), method='linear')
                            nan_mask = np.isnan(nodal_age)
                            if np.any(nan_mask):
                                df = pd.DataFrame(nodal_age)
                                df = df.interpolate(method='linear', limit_direction='both', axis=0)
                                df = df.interpolate(method='linear', limit_direction='both', axis=1)
                                df = df.ffill(axis=0).bfill(axis=0).ffill(axis=1).bfill(axis=1)
                                nodal_age = df.to_numpy()
                            
                            # Mask out nodes if all markers in their elements are NaN
                            node_is_nan = np.ones((nx, nz), dtype=bool)
                            valid_i = i_el[valid]
                            valid_j = j_el[valid]
                            node_is_nan[valid_i, valid_j] = False
                            node_is_nan[valid_i + 1, valid_j] = False
                            node_is_nan[valid_i, valid_j + 1] = False
                            node_is_nan[valid_i + 1, valid_j + 1] = False
                            
                            nodal_age[node_is_nan] = np.nan
                            nodal_ages[c] = nodal_age
                        else:
                            nodal_ages[c] = np.full((nx, nz), np.nan)

        if i < start:
            continue

        print('Writing record #%d, model time=%.3e' % (i, fl.time[i-1]), end='\r')
        sys.stdout.flush()
        fvts = open('flac.%06d.vts' % i, 'w')
        vts_header(fvts, nex, nez, fl.time[i-1], fl.steps[i-1])

        # node-based field
        fvts.write('  <PointData>\n')

        vx, vz = fl.read_vel(i)
        # VTK requires vector field (velocity, coordinate) has 3 components.
        # Allocating a 3-vector tmp array for VTK data output.
        tmp = np.zeros((fl.nx, fl.nz, 3), dtype=vx.dtype)
        tmp[:,:,0] = vx
        tmp[:,:,1] = vz
        # vts requires x-axis increment fastest, swap axes order
        vts_dataarray(fvts, tmp.swapaxes(0,1), 'Velocity', 3)

        if thermochron:
            a = T_nodes
        else:
            a = fl.read_temperature(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Temperature')

        if thermochron:
            for name, data in zip(thermochron_system_names, nodal_ages):
                vts_dataarray(fvts, data.swapaxes(0, 1), f'age_{name}')

        x0, z0 = fl.read_original_mesh(i)
        vts_dataarray(fvts, x0.swapaxes(0,1), 'x0')
        vts_dataarray(fvts, z0.swapaxes(0,1), 'z0')

        fvts.write('  </PointData>\n')

        # element-based field
        fvts.write('  <CellData>\n')

        # logrithm of strain rate 2nd invariant
        a = fl.read_srII(i)
        srat = a
        vts_dataarray(fvts, a.swapaxes(0,1), 'Strain rate')

        srxx, srzz, srxz = fl.read_strain_rate(i)
        vts_dataarray(fvts, srxx.swapaxes(0,1), 'Sr xx')
        vts_dataarray(fvts, srzz.swapaxes(0,1), 'Sr zz')
        vts_dataarray(fvts, srxz.swapaxes(0,1), 'Sr xz')

        sr1 = compute_s1(srxx, srzz, srxz)
        vts_dataarray(fvts, sr1.swapaxes(0,1), 'Strain rate 1-axis', 3)

        a = fl.read_eII(i)
        eii = a
        vts_dataarray(fvts, a.swapaxes(0,1), 'eII')

        exx, ezz, exz = fl.read_strain(i)
        e1 = compute_s1(exx, ezz, exz)
        vts_dataarray(fvts, e1.swapaxes(0,1), 'Strain 1-axis', 3)

        a = fl.read_density(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Density')

        a = fl.read_area(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Area')

        a = fl.read_aps(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Plastic strain')

        a = fl.read_sII(i)
        sii = a
        vts_dataarray(fvts, a.swapaxes(0,1), 'Stress')

        sxx = fl.read_sxx(i)
        vts_dataarray(fvts, sxx.swapaxes(0,1), 'Sxx')

        syy = fl.read_szz(i)
        vts_dataarray(fvts, syy.swapaxes(0,1), 'Syy')

        szz = fl.read_szz(i)
        vts_dataarray(fvts, szz.swapaxes(0,1), 'Szz')

        sxz = fl.read_sxz(i)
        vts_dataarray(fvts, sxz.swapaxes(0,1), 'Sxz')

        pressure = fl.read_pres(i)
        vts_dataarray(fvts, pressure.swapaxes(0,1), 'Pressure')

        # compression axis of stress
        a = compute_s1(sxx, szz, sxz)
        vts_dataarray(fvts, a.swapaxes(0,1), 's1', 3)

        a = fl.read_fmelt(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Melt fraction')

        a = fl.read_fmagma(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Magma fraction')

        a = fl.read_diss(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Dissipation')

        # logrithm of effective viscosity
        eff_visc = np.log10(a + 1e-45) + 8 - srat
        vts_dataarray(fvts, eff_visc.swapaxes(0,1), 'Eff. Visc')

        a = fl.read_visc(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Viscosity')

        a = fl.read_phase(i)
        vts_dataarray(fvts, a.swapaxes(0,1), 'Phase')

        # Work done by stress
        a = sii * 1e8 * eii
        vts_dataarray(fvts, a.swapaxes(0,1), 'Work')

        fvts.write('  </CellData>\n')

        # coordinate
        if not thermochron:
            x, z = fl.read_mesh(i)
        tmp[:,:,0] = x
        tmp[:,:,1] = z
        fvts.write('  <Points>\n')
        vts_dataarray(fvts, tmp.swapaxes(0,1), '', 3)
        fvts.write('  </Points>\n')

        vts_footer(fvts)
        fvts.close()

    print()
    return


def compute_s1(sxx, szz, sxz):
    mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
    theta = 0.5 * np.arctan2(2*sxz,  sxx-szz)

    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    nx, nz = sxx.shape
    tmp = np.zeros((nx, nz, 3), dtype=sxx.dtype)
    tmp[:,:,0] = mag * np.sin(theta)
    tmp[:,:,1] = mag * np.cos(theta)
    return tmp


def vts_dataarray(f, data, data_name=None, data_comps=None):
    if data.dtype in (int, np.int32, np.int_):
        dtype = 'Int32'
        data = data.astype(np.int32)
    elif data.dtype in (float, np.single, np.double,
                        np.float32, np.float64, np.float128):
        dtype = 'Float32'
        data = data.astype(np.float32)
    else:
        raise Error('Unknown data type: ' + name)

    name = ''
    if data_name:
        name = 'Name="{0}"'.format(data_name)

    ncomp = ''
    if data_comps:
        ncomp = 'NumberOfComponents="{0}"'.format(data_comps)

    if output_in_binary:
        fmt = 'binary'
    else:
        fmt = 'ascii'
    header = '<DataArray type="{0}" {1} {2} format="{3}">\n'.format(
        dtype, name, ncomp, fmt)
    f.write(header)

    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tobytes()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header).decode('ascii'))
        f.write(base64.standard_b64encode(b).decode('ascii'))
    else:
        data.tofile(f, sep=' ')
    f.write('\n</DataArray>\n')
    return


def vts_header(f, nex, nez, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<StructuredGrid WholeExtent="0 {0} 0 {1} 0 0">
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {3}
  </DataArray>
</FieldData>
<Piece Extent="0 {0} 0 {1} 0 0">
'''.format(nex, nez, time, step))
    return


def vts_footer(f):
    f.write(
'''</Piece>
</StructuredGrid>
</VTKFile>
''')
    return


if __name__ == '__main__':
    thermochron = False
    args = sys.argv[1:]
    if '-t' in args:
        thermochron = True
        args.remove('-t')

    if len(args) < 1:
        print('''usage: flac2vtk.py [-t] path [frame_min [frame_max]]

Processing flac data output to VTK format.

If -t is specified, compute thermochronology closure ages of markers and write them to the vts files.
If frame_min is -1, start from the latest vts file.
If frame_max is -1 or not given, processing to latest frames.
If both frame_min and frame_max are not given, processing all frames''')
        sys.exit(1)

    path = args[0]
    start = 1
    end = -1
    if len(args) >= 2:
        start = int(args[1])
        if len(args) >= 3:
            end = int(args[2])

    main(path, start, end, thermochron=thermochron)
