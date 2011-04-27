#!/usr/bin/env python

import sys
try:
    import numpy as np
except ImportError:
    print 'Error: Failed importing "numpy" module.'
    print 'Please install the module and add it to PYTHONPATH environment variable.'
    sys.exit(1)

try:
    # python 2.7 or later
    from collections import Counter
except ImportError:
    try:
        # backport from http://code.activestate.com/recipes/576611
        from counter import Counter
    except ImportError:
        pass


# the precision of saved file
doubleprecision = False
sizeofint = 4
if doubleprecision:
    sizeoffloat = 8
    default_dtype = np.double
else:
    sizeoffloat = 4
    default_dtype = np.single



class Flac(object):
    '''Read Flac data file
    '''

    def __init__(self, swap_endian=False):
        self.swap_endian = swap_endian
        self.reload()
        return


    def reload(self):
        # read record number
        tmp = np.fromfile('_contents.0', sep=' ')
        tmp.shape = (-1, 3)
        self.time = tmp[:,2]
        self.nrec = len(self.time)

        # read grid size
        nex, nez = np.fromfile('nxnz.0', sep=' ', dtype=int)
        self.nx, self.nz = nex+1, nez+1
        self.nnodes = self.nx * self.nz
        self.nelements = nex * nez

        return


    def read_mesh(self, frame):
        columns = 2
        f = open('mesh.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        x, z = self._read_data(f, columns)
        self._reshape_nodal_fields(x, z)
        return x, z


    def read_vel(self, frame):
        columns = 2
        f = open('vel.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        vx, vz = self._read_data(f, columns)
        self._reshape_nodal_fields(vx, vz)
        return vx, vz


    def read_temperature(self, frame):
        columns = 1
        f = open('temperature.0')
        offset = (frame-1) * columns * self.nnodes * sizeoffloat
        f.seek(offset)
        T = self._read_data(f, columns)
        self._reshape_nodal_fields(T)
        return T


    def read_aps(self, frame):
        columns = 1
        f = open('aps.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        aps = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(aps)
        return aps


    def read_density(self, frame):
        columns = 1
        f = open('density.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        density = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(density)
        return density



    def read_eII(self, frame):
        columns = 1
        f = open('eII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        eII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(eII)
        return eII


    def read_sII(self, frame):
        columns = 1
        f = open('sII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sII)
        return sII


    def read_sxx(self, frame):
        columns = 1
        f = open('sxx.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sxx = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sxx)
        return sxx


    def read_sxz(self, frame):
        columns = 1
        f = open('sxz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        sxz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(sxz)
        return sxz


    def read_szz(self, frame):
        columns = 1
        f = open('szz.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        szz = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(szz)
        return szz


    def read_srII(self, frame):
        columns = 1
        f = open('srII.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        srII = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(srII)
        return srII


    def read_pres(self, frame):
        columns = 1
        f = open('pres.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        pres = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(pres)
        return pres


    def read_diss(self, frame):
        columns = 1
        f = open('diss.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        diss = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(diss)
        return diss


    def read_visc(self, frame):
        columns = 1
        f = open('visc.0')
        offset = (frame-1) * columns * self.nelements * sizeoffloat
        f.seek(offset)
        visc = self._read_data(f, columns, count=self.nelements)
        self._reshape_elemental_fields(visc)
        return visc


    def read_phase(self, frame):
        columns = 1
        f = open('phase.0')
        offset = (frame-1) * columns * self.nelements * sizeofint
        f.seek(offset)
        phase = self._read_data(f, columns, count=self.nelements, dtype=np.int32)
        self._reshape_elemental_fields(phase)
        return phase


    def read_markers(self, frame):
        # read tracer size
        tmp = np.fromfile('_markers.0', sep=' ')
        tmp.shape = (-1, 4)
        n = int(tmp[frame-1,2])

        suffix = '.%06d.0' % frame
        dead = self._read_data('markdead' + suffix, count=n, dtype=np.int32)

        tmp = self._read_data('markx' + suffix, count=n)
        x = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data('marky' + suffix, count=n)
        y = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data('markage' + suffix, count=n)
        age = self._remove_dead_markers(tmp, dead)

        tmp = self._read_data('markphase' + suffix, count=n, dtype=np.int32)
        phase = self._remove_dead_markers(tmp, dead)

        return x, y, age, phase


    def read_tracers(self):
        # read tracer size
        tmp = np.fromfile('_tracers.0', sep=' ')
        tmp.shape = (-1, 4)
        ntracerrec = tmp.shape[0]
        ntracers = int(tmp[0,1])
        n = ntracerrec * ntracers

        time = self._read_data('outtracktime.0', count=n)
        x = self._read_data('outtrackxx.0', count=n)
        x.shape = (ntracerrec, ntracers)

        y = self._read_data('outtrackyy.0', count=n)
        y.shape = (ntracerrec, ntracers)

        T = self._read_data('outtracktemp.0', count=n)
        T.shape = (ntracerrec, ntracers)

        p = self._read_data('outtrackpres.0', count=n)
        p.shape = (ntracerrec, ntracers)

        e = self._read_data('outtrackstrain.0', count=n)
        e.shape = (ntracerrec, ntracers)

        phase = self._read_data('outtrackphase.0', count=n)
        phase.shape = (ntracerrec, ntracers)

        return x, y, T, p, e, phase


    def _read_data(self, fileobj, columns=1,
                   count=None, dtype=None):
        '''Read data from a file-like object 'fileobj'.

        The 'dtype' specifies the storage type, default to single precision
        float.
        '''

        # number of nodes
        if count is None:
            count = self.nnodes

        # total number of items
        n = columns * count

        if dtype is None:
            dtype = default_dtype

        result = np.fromfile(fileobj, dtype, n)
        if self.swap_endian:
            result.byteswap(True)

        if columns == 1:
            return result
        else:
            result.shape = (columns, -1)
            return tuple(result[i,:] for i in range(columns))


    def _reshape_nodal_fields(self, *argv):
        for x in argv:
            x.shape = (self.nx, self.nz)
        return


    def _reshape_elemental_fields(self, *argv):
        for x in argv:
            x.shape = (self.nx-1, self.nz-1)
        return


    def _remove_dead_markers(self, a, dead):
        b = a[dead==1]
        return b


def elem_coord(x, y):
    '''Turning nodal coordinates to element coordinates'''
    cx = (x[:-1, :-1] + x[:-1, 1:] + x[1:, :-1] + x[1:, 1:]) / 4
    cy = (y[:-1, :-1] + y[:-1, 1:] + y[1:, :-1] + y[1:, 1:]) / 4
    return cx, cy


def make_uniform_grid(xmin, xmax, ymin, ymax, dx, dy):
    # grid size
    nx = (xmax - xmin) / dx + 1
    ny = (ymax - ymin) / dy + 1

    # generate uniform grid
    xx = np.linspace(xmin, xmax, nx)
    yy = np.linspace(ymin, ymax, ny)

    x, y = np.meshgrid(xx, yy)
    return x, y


def nearest_neighbor_interpolation2d(x0, y0, f0, x, y):
    '''Interpolating field f0, which is defined on (x0, y0)
    to a new grid (x, y) using nearest neighbor method'''

    if x0.shape != y0.shape:
        raise Exception('x0 and y0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != y.shape:
        raise Exception('x and y arrays have different shape')

    dx = x[0,1] - x[0,0]
    dy = y[1,0] - y[0,0]

    ny, nx = x.shape
    f = np.zeros(x.shape)
    for i in range(ny):
        for j in range(nx):
            dist2 = ((x[i,j] - x0) / dx)**2 + ((y[i,j] - y0) / dy)**2
            ind = np.argmin(dist2)
            f[i,j] = f0[ind]

    return f


def neighborhood_interpolation2d(x0, y0, f0, x, y, dx, dy):
    '''Interpolating field f0, which is defined on (x0, y0)
    to a new grid (x, y) using neighborhood'''

    if f0.dtype not in (bool, np.bool,
                        int, np.int, np.int8, np.int16, np.int32, np.int64):
        raise Exception('f0 must have integer value')

    if x0.shape != y0.shape:
        raise Exception('x0 and y0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != y.shape:
        raise Exception('x and y arrays have different shape')

    ny, nx = x.shape
    f = np.zeros(x.shape, dtype=f0.dtype)
    for i in range(ny):
        for j in range(nx):
            #dist2 = ((x[i,j] - x0) / dx)**2 + ((y[i,j] - y0) / dy)**2
            #ind = np.argmin(dist2)
            xx = x[i,j]
            yy = y[i,j]
            ind = (x0 >= xx-dx) * (x0 <= xx+dx) * (y0 >= yy-dy) * (y0 <= yy+dy)
            g = f0[ind]
            if len(g) == 0:
                #print (i, j), (xx, yy), ind
                #raise Exception('No point inside domain bounds. Increase (dx, dy)!')
                f[i,j] = sys.maxint
            else:
                z = Counter(g).most_common(1)
                f[i,j] = z[0][0]

    return f


def gaussian_interpolation2d(x0, y0, f0, x, y):
    '''Interpolating field f0, which is defined on (x0, y0)
    to a new grid (x, y) using nearest neighbor method'''

    if x0.shape != y0.shape:
        raise Exception('x0 and y0 arrays have different shape')

    if x0.shape != f0.shape:
        raise Exception('x0 and f0 arrays have different shape')

    if x.shape != y.shape:
        raise Exception('x and y arrays have different shape')

    # using 1d index for x0, y0, f0
    x0 = x0.flat
    y0 = y0.flat
    f0 = f0.flat

    dx = 1.5 * (x[0,1] - x[0,0])
    dy = 1.5 * (y[1,0] - y[0,0])

    f = np.zeros(x.shape)
    g = np.zeros(x.shape)
    for i in range(len(x0)):
        weight = np.exp(-((x - x0[i]) / dx)**2 - ((y - y0[i]) / dy)**2)
        f += weight * f0[i]
        g += weight

    return f / g


def printing(*args, **kwd):
    '''printing(arg1, arg2, ..., stream=None)

    stream: None -- output to stdout
            filename -- output to file
            file-like object -- output to the object
    '''

    stream = kwd.get('stream')
    if stream == None:
        stream = sys.stdout
    elif isinstance(stream, (str, unicode)):
        filename = stream
        stream = open(filename, 'w')

    # flatten numpy arrays
    try:
        args = tuple(x.flat for x in args)
    except AttributeError:
        pass

    narg = len(args)
    fmt = '%.15e' + '\t%.15e'*(narg-1)

    for items in zip(*args):
        print >> stream, fmt % tuple(items)
    return



## This is an example on how to use this module.
## Before running this module, 'cd' to a directory containing the flac data.
if __name__ == '__main__':

    fl = Flac()

    # read the last record of the mesh and temperature
    x, y = fl.read_mesh(fl.nrec-1)
    T = fl.read_temperature(fl.nrec-1)

    # print (x, y, T) to screen
    printing(x, y, T)

    print '# time =', fl.time[fl.nrec-1], 'Myrs'
