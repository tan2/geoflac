#!/usr/bin/env python
'''
Read flac density and coordinate data and return gravity and topography at regular spacing.

Similar to flac_gravity.py, but interpolate density to finer grid before computing gravity.
The resulting gravity field is more oscillating than that of flac_gravity.py
'''

import sys
import numpy as np
import flac
import flac_interpolate as fi

gn = 6.6726e-11
ocean_density = 1030.


def compute_gravity(frame, xresl, yresl):
    '''Interploate density to a uniform grid before computing gravity.
    xresl and yresl are the spacing of the new uniform grid in km'''

    ## read data
    fl = flac.Flac()
    xorig, yorig = fl.read_mesh(frame)  # in km

    # new refined uniform grid
    xmin = xorig[0,0]
    xmax = xorig[-1,0]
    fi.xmin = xmin
    fi.xmax = xmax
    d = 5 # in km
    fi.ymin = np.floor(yorig.min()/d + 1) * d
    fi.ymax = np.floor(yorig.max()/d + 1) * d
    fi.dx = xresl
    fi.dy = yresl

    # interpolate density to new grid
    #'''
    x, y, rho = fi.interpolate(frame, 'density')

    # km to meter
    x *= 1e3
    y *= 1e3

    # fill water and air
    ind = (y < 0) * rho.mask
    rho.mask[ind] = False
    rho.data[ind] = ocean_density
    rho = rho.filled(0)
    #'''

    xresl *= 1e3
    yresl *= 1e3
    xmax *= 1e3
    xmin *= 1e3

    area = xresl * yresl
    mass = rho * area

    ## calculate gravity at these points
    # px in uniform spacing
    px = x[0,:]
    # py is 4km above the highest topography to avoid high frequency oscillation
    py_height = y.max() + 4e3
    py = np.ones(px.shape) * py_height

    # original topography defined in px grid
    topo = np.interp(px, x[:,0], y[:,0])

    ## contribution of material inside the model
    grav = np.empty(px.shape)
    for i in range(len(grav)):
        dx = px[i] - x
        dy = py[i] - y
        dist2 = (dx**2 + dy**2)
        # downward component of gravity of an infinite long line source of density anomaly
        # see Turcotte & Schubert, 1st Ed., Eq 5-104
        grav[i] = 2 * gn * np.sum(mass * dy / dist2)

    # contribution of material outside left boundary
    # assuming the leftmost element extend to negative infinity
    # see Turcotte & Schubert, 1st Ed., Eq 5-106
    nz = y.shape[-1]
    for i in range(nz):
        sigma = rho[0,i] * yresl
        dx = px - xmin
        dy = py - y[0,i]
        angle = np.arctan2(dx, dy)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    # ditto for right boundary
    for i in range(nz):
        sigma = rho[-1,i] * yresl
        dx = xmax - px
        dy = py - y[-1,i]
        angle = np.arctan2(dx, dy)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    ##
    grav -= np.mean(grav)

    return px, topo, grav


if __name__ == '__main__':
    frame = int(sys.argv[1])
    xresl = float(sys.argv[2])
    yresl = float(sys.argv[3])

    x, topo, gravity = compute_gravity(frame, xresl, yresl)

