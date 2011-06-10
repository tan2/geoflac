#!/usr/bin/env python
'''Interpolate flac fields to a high resolution regular grid suitable for plotting

Usage: flac_interpolate.py frame field

frame (starts from 1)
field can be either 'phase', 'temperature', 'aps', 'density', 'eII', 'sII',
                    'sxx', 'szz', 'sxz', 'srII', 'pres', 'diss', 'visc'
'''
import sys, os
import numpy as np
import flac


# domain bounds (in km)
xmin = 300
xmax = 700
ymin = -150
ymax = 10

# grid resolution (in km)
dx = .8
dy = .2


def excluding(x, y, f, xmin, xmax, ymin, ymax):
    '''Excluding nodes of (x, y, f) that are outside the domain bounds'''
    # bool * bool is element-wise logical AND
    ind = (xmin <= x) * (x <= xmax) * (ymin <= y) * (y <= ymax)
    return x[ind], y[ind], f[ind]


def clip_topo(x, y, f, x0, y0):
    '''Setting f=NaN for nodes above (x0, y0).
    x, y, f x0, y0 are 2d arrays.'''
    xx0 = x0[:,0]
    yy0 = y0[:,0]
    xx = x[0,:]
    yy = np.interp(xx, xx0, yy0)

    for i in range(len(xx)):
        ind = y[:,i] > yy[i]
        f[ind,i] = np.nan
    return np.ma.masked_array(f, mask=np.isnan(f))


def interpolate(frame, field):
    # read coordinate of nodes
    fl = flac.Flac()
    xx, yy = fl.read_mesh(frame)
    x, y = flac.make_uniform_grid(xmin, xmax, ymin, ymax, dx, dy)

    if field == 'phase':
        ## phase
        # read marker location, age and phase
        mx, my, mage, mphase = fl.read_markers(frame)
        mx, my, mphase = excluding(mx, my, mphase, xmin-dx, xmax+dx, ymin-dy, ymax+dy)
        ph = flac.nearest_neighbor_interpolation2d(mx, my, mphase, x, y)
        #ph = flac.neighborhood_interpolation2d(mx, my, mphase, x, y, 7*dx, 7*dy)
        f = ph.astype(np.float32)
    elif field in ('temperature', 'aps', 'density', 'eII', 'sII',
                   'sxx', 'szz', 'sxz', 'srII', 'pres', 'diss', 'visc'):
        # read field
        f = getattr(fl, 'read_'+field)(frame)
        cx, cy = flac.elem_coord(xx, yy)
        cx, cy, f = excluding(cx, cy,f, xmin-dx, xmax+dx, ymin-dy, ymax+dy)
        f = flac.gaussian_interpolation2d(cx, cy, f, x, y)
    else:
        raise RuntimeError('unknown field %s' % field)

    f = clip_topo(x, y, f, xx, yy)
    return cx, cy, f


if __name__ == '__main__':

    if len(sys.argv) < 3:
        print __doc__
        exit(1)

    frame = int(sys.argv[1])
    field = sys.argv[2].lower()
    cx, cy, f = interpolate(frame, field)
    flac.printing(cx, cy, f)
