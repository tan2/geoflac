#!/usr/bin/python

import sys, os
import numpy as np
from math import sqrt
from scipy.special import erf
from matplotlib import pyplot as plt


def half_space_cooling_T(z, Tsurf, Tmantle,  age_in_myrs):
    diffusivity = 1e-6
    myrs2sec = 86400 * 365.2425e6

    T = Tsurf + (Tmantle - Tsurf) * erf(z /
            sqrt(4 * diffusivity * age_in_myrs * myrs2sec) )
    return T


def get_visc(edot, T, n, A, E):
    '''edot: second invariant of strain rate
    T: temperature in Celsius
    n, A, E: viscosity parameters

    return viscosity in Pascal.s
    '''
    R = 8.31448  # gas constant
    pow = 1.0/n - 1
    pow1 = -1.0/n
    visc = 0.25 * (edot**pow) * (0.75*A)**pow1 * np.exp(E / (n * R * (T + 273))) * 1e6
    return visc


def visc_profile(z, T, edot, layerz, nAEs):
    '''Viscosity profile of multi-layers
    z: numpy array of depth (in meters)
    T: array of temperature (in Celsius)
    edot: strain rate (in 1/second)
    layerz: (0, z1, z2, ...) the depth interface of the layers
    nAEs: ( ..., (n, A, E), ...) visc parameters of each layers
    '''

    if layerz[0] != 0:
        print("Error: layerz[0] is not 0", layerz)
    nlayers = len(layerz)
    layerz = tuple(layerz) + (z[-1],)  # deepest depth

    viscp = np.zeros_like(z)
    for i in range(nlayers):
        n, A, E = nAEs[i][:]
        vs = get_visc(edot, T, n, A, E)

        # find depth range of each layer
        z0, z1 = layerz[i], layerz[i+1]
        n0 = (z >= z0).argmax()
        n1 = (z >= z1).argmax()
        #print(i, z0, fz1, n0, n1)

        viscp[n0:n1] = vs[n0:n1]
    return viscp



if __name__ == "__main__":

    layerz = (0, 15e3, 30e3)   # 1st elem must be 0
    nAEs = ( (3.05, 1.13e+2, 2.00e+5),
             (2.00, 1.00e-5, 1.67e+5),
             (3.00, 7.00e+4, 5.10e+5) )
    edot = 1e-14  # high strain rate
    #edot = 1e-15  # low strain rate

    # upper bound of z
    deepz = layerz[-1] * 2

    z = np.linspace(0, deepz, num=101)
    T = half_space_cooling_T(z, 10, 1330, 150)

    visc = visc_profile(z, T, edot, layerz, nAEs)

