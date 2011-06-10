#!/usr/bin/env python

import sys
import numpy as np
import flac

gn = 6.6726e-11
ocean_density = 1030.
free_air_factor = 0.3086e-5  # in m/s^2


def compute_gravity(frame):

    ## read data
    fl = flac.Flac()
    x, y = fl.read_mesh(frame)  # in km
    x *= 1e3
    y *= 1e3
    rho = fl.read_density(frame)  # in kg/m^3

    # surface coordinates
    xx = x[:,0]
    yy = y[:,0]
    xmin = x[0,0]
    xmax = x[-1,0]

    # center of elements
    cx, cy = flac.elem_coord(x, y)

    # area of elements = 0.5 * | AC x BD | = 0.5 * |xdiag1*ydiag2 - xdiag2*ydiag1|
    #  A -- D
    #  |    |
    #  B -- C
    xdiag1 = x[0:-1, 0:-1] - x[1:, 1:]
    ydiag1 = y[0:-1, 0:-1] - y[1:, 1:]
    xdiag2 = x[1:, 0:-1] - x[0:-1, 1:]
    ydiag2 = y[1:, 0:-1] - y[0:-1, 1:]
    area = 0.5 * np.abs(xdiag1*ydiag2 - xdiag2*ydiag1)
    mass = rho * area


    ## calculate gravity at these points
    # px in uniform spacing
    px = np.linspace(xmin, xmax, num=5*fl.nx)
    # py is 2km above the highest topography to avoid high frequency oscillation
    py_height = np.max(yy) + 2e3
    py = np.ones(px.shape) * py_height

    # original topography defined in px grid
    #origy = np.interp(px, x[:,0], y[:,0])
    #py = origy.copy()
    #py[py<0] = 0
    #py += 2e3

    ## contribution of material inside the model
    grav = np.empty(px.shape)
    for i in range(len(grav)):
        dx = px[i] - cx
        dy = py[i] - cy
        dist2 = (dx**2 + dy**2)
        # downward component of gravity of an infinite long line source of density anomaly
        # see Turcotte & Schubert, 1st Ed., Eq 5-104
        grav[i] = 2 * gn * np.sum(mass * dy / dist2)

    ## contribution of ocean
    for i in range(fl.nx-1):
        midy = 0.5 * (yy[i] + yy[i+1])
        if midy < 0:
            midx = 0.5 * (xx[i] + xx[i+1])
            m = (xx[i+1] - xx[i]) * -midy * ocean_density
            dx = px - midx
            dy = py - midy
            dist2 = (dx**2 + dy**2)
            grav += 2 * gn * m * dy / dist2

    # contribution of material outside left boundary
    # assuming the leftmost element extend to negative infinity
    # see Turcotte & Schubert, 1st Ed., Eq 5-106
    for i in range(fl.nz-1):
        sigma = rho[0,i] * (y[0,i] - y[0,i+1])
        dx = px - xmin
        dy = py - 0.5 * (y[0,i] + y[0,i+1])
        angle = np.arctan2(dx, dy)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if yy[0] < 0:
        sigma = ocean_density * -yy[0]
        dx = px - xmin
        dy = py - 0.5 * yy[0]
        angle = np.arctan2(dx, dy)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)


    # ditto for right boundary
    for i in range(fl.nz-1):
        sigma = rho[-1,i] * (y[-1,i] - y[-1,i+1])
        dx = xmax - px
        dy = py - 0.5 * (y[-1,i] + y[-1,i+1])
        angle = np.arctan2(dx, dy)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    if yy[-1] < 0:
        sigma = ocean_density * -yy[-1]
        dx = xmax - px
        dy = py - 0.5 * yy[-1]
        angle = np.arctan2(dx, dy)
        grav += 2 * gn * sigma * (0.5*np.pi - angle)

    ## free-air correction
    grav += free_air_factor * py

    ##
    grav -= np.mean(grav)

    return px, grav


if __name__ == '__main__':
    frame = int(sys.argv[1])
    x, gravity = compute_gravity(frame)


