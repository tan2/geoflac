
import numpy as np
import flac

fl = flac.Flac()

frame = 1

x, z = fl.read_mesh(frame)
T = fl.read_temperature(frame)
ph = fl.read_phase(frame)

xm, zm, age, phase, idm, a1, a2, ntriag = fl.read_markers(frame)

# test: xx and xm should be very close
xx = flac.marker_interpolate_node(ntriag, a1, a2, fl.nz, x)
print( np.max(np.abs(xx - xm) / xm) )  # relative error, should be small

# nodal field
Tm = flac.marker_interpolate_node(ntriag, a1, a2, fl.nz, T)

# elem field
pm = flac.marker_interpolate_elem(ntriag, fl.nz, ph)

