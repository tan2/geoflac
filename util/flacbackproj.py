#!/usr/bin/python3

'''Back-project marker output of flac to VTK (vtp) files.
Usage: flacbackproj.py <path>
'''


import sys, os, glob
import numpy as np
import flac
from flac2vtk import vts_dataarray
from flacmarker2vtk import vtp_header, vtp_footer

path = sys.argv[1]
os.chdir(path)

try:
    fl= flac.Flac()
    m0 = fl.read_markers(0)
except:
    fl = flac.FlacFromVTK()
    m0 = fl.read_markers(0)

x_orig = m0[0]
z_orig = m0[1]
id_orig = m0[4]

for frame in fl.frames:
    i = frame

    m_new = fl.read_markers(frame-1)
    x_new = m_new[0]
    z_new = m_new[1]
    id_new = m_new[4]

    # markers existed in m0
    s = id_new <= id_orig.max()

    si = id_new[s] - 1
    z0 = z_orig[si]
    x0 = x_orig[si]

    nmarkers = si.size

    print('Writing record #%d, model time=%.3e, %d markers' % (i, fl.time[i-1], nmarkers))
    fvtp = open('flacbp.%06d.vtp' % i, 'w')
    vtp_header(fvtp, nmarkers, fl.time[i-1], fl.steps[i-1])

    # point-based data
    fvtp.write('  <PointData>\n')
    vts_dataarray(fvtp, x0, 'x0', 1)
    vts_dataarray(fvtp, z0, 'z0', 1)
    fvtp.write('  </PointData>\n')

    # point coordinates

    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    tmp = np.zeros((nmarkers, 3), dtype=x_new.dtype)
    tmp[:,0] = x_new[s]
    tmp[:,1] = z_new[s]
    fvtp.write('  <Points>\n')
    vts_dataarray(fvtp, tmp, '', 3)
    fvtp.write('  </Points>\n')

    vtp_footer(fvtp)
    fvtp.close()

