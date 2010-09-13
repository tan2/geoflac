#!/usr/bin/env python

import sys

# the precision of saved file
doubleprecision = False
if doubleprecision:
    sizeoffloat = 8
    default_typecode = 'd'
else:
    sizeoffloat = 4
    default_typecode = 'f'



class FlacBase(object):
    '''Base class for reading Flac files. _read_data() needs to be
    defined by derived class.'''

    def __init__(self, swap_endian=False):
        self.swap_endian = swap_endian

        # read record number
        headers = open('_contents.0').readlines()
        header = headers[-1].split()
        self.nrec = int(header[0])

        # read grid size
        header = open('nxnz.0').readline().split()
        self.nx = int(header[0]) + 1
        self.nz = int(header[1]) + 1

        self.nodes = self.nx * self.nz
        self.elements = (self.nx - 1) * (self.nz - 1)

        # read tracer size
        headers = open('_tracers.0').readlines()
        header = headers[-1].split()
        self.ntracerrec = int(header[0])
        self.ntracers = int(header[1])
        return


    def read_mesh(self, frame):
        columns = 2
        f = open('mesh.0')
        offset = frame * columns * self.nodes * sizeoffloat
        f.seek(offset)
        self.x, self.z = self._read_data(f, columns)
        return


    def read_vel(self, frame):
        columns = 2
        f = open('vel.0')
        offset = frame * columns * self.nodes * sizeoffloat
        f.seek(offset)
        self.vx, self.vz = self._read_data(f, columns)
        return


    def read_temperature(self, frame):
        columns = 1
        f = open('temperature.0')
        offset = frame * columns * self.nodes * sizeoffloat
        f.seek(offset)
        self.T = self._read_data(f, columns)
        return


    def read_time(self):
        f = open('time.0')
        self.time = self._read_data(f, 1, num=self.nrec)

        return


    def ij2node(self, i, j):
        return (j-1) + (i-1) * (self.nz)



class _Flac(FlacBase):
    '''Read Flac data file using stdlib "array" module. The data will be read in as 1D arrays.'''

    def _read_data(self, fileobj, columns,
                    num=None, typecode=None):
        '''Read data from a file-like object 'fileobj'.

        The 'typecode' specifies the storage type, default to single precision
        float. Other typecode is possible, but not implemented.
        '''

        # number of nodes
        if num is None:
            num = self.nodes

        # total number of items
        n = columns * num

        if typecode is None:
            typecode = default_typecode

        a = array.array(typecode)
        a.fromfile(fileobj, n)
        if self.swap_endian:
            a.byteswap()

        #print a
        if columns == 1:
            return a
        else:
            # split the data
            result = range(columns)
            for i in range(columns):
                result[i] = a[i*num:(i+1)*num]

            # convert list into tuple for unpacking
            return tuple(result)



class _Flac_numpy(FlacBase):
    '''Read Flac data file using "numpy" module.
    This class is preferred over class _Flac.
    '''


    def _read_data(self, fileobj, columns,
                    num=None, typecode=None):
        '''Read data from a file-like object 'fileobj'.

        The 'dtype' specifies the storage type, default to single precision
        float.
        '''

        # number of nodes
        if num is None:
            num = self.nodes

        # total number of items
        n = columns * num

        if typecode is None:
            typecode = default_typecode

        result = numpy.fromfile(fileobj, typecode, n)
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


    def read_mesh(self, frame):
        FlacBase.read_mesh(self, frame)
        self._reshape_nodal_fields(self.x, self.z)
        return


    def read_vel(self, frame):
        FlacBase.read_vel(self, frame)
        self._reshape_nodal_fields(self.vx, self.vz)
        return


    def read_temperature(self, frame):
        FlacBase.read_temperature(self, frame)
        self._reshape_nodal_fields(self.T)
        return


    def read_tracers(self):
        self.tracer_time = numpy.zeros(self.ntracerrec, dtype='f')
        self.tracer_x = numpy.zeros((self.ntracerrec, self.ntracers), dtype='f')
        self.tracer_y = numpy.zeros((self.ntracerrec, self.ntracers), dtype='f')
        self.tracer_temp = numpy.zeros((self.ntracerrec, self.ntracers), dtype='f')
        self.tracer_pres = numpy.zeros((self.ntracerrec, self.ntracers), dtype='f')
        self.tracer_strain = numpy.zeros((self.ntracerrec, self.ntracers), dtype='f')
        self.tracer_phase = numpy.zeros((self.ntracerrec, self.ntracers), dtype='f')

        self.tracer_id = self._read_data('outtrackID.0', 1, num=self.ntracers, typecode='f')

        f1 = open('outtracktime.0')
        f2 = open('outtrackxx.0')
        f3 = open('outtrackyy.0')
        f4 = open('outtracktemp.0')
        f5 = open('outtrackpres.0')
        f6 = open('outtrackstrain.0')
        f7 = open('outtrackphase.0')

        for i in range(self.ntracerrec):
            time = self._read_data(f1, 1, num=self.ntracers, typecode='f')
            x = self._read_data(f2, 1, num=self.ntracers, typecode='f')
            y = self._read_data(f3, 1, num=self.ntracers, typecode='f')
            temperature = self._read_data(f4, 1, num=self.ntracers, typecode='f')
            pressure = self._read_data(f5, 1, num=self.ntracers, typecode='f')
            strain = self._read_data(f6, 1, num=self.ntracers, typecode='f')
            phase = self._read_data(f7, 1, num=self.ntracers, typecode='f')

            self.tracer_x[i,:] = x
            self.tracer_y[i,:] = y
            self.tracer_temp[i,:] = temperature
            self.tracer_pres[i,:] = pressure
            self.tracer_strain[i,:] = strain
            self.tracer_phase[i,:] = phase

            # all tracers have the same time
            self.tracer_time[i] = time[0]

        return




## setup alias
try:
    import numpy
    Flac = _Flac_numpy
except ImportError:
    import array
    Flac = _Flac



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

    # read the last record from the mesh file
    fl.read_mesh(fl.nrec-1)

    # print x, z to screen
    printing(fl.x, fl.z)

    fl.read_time()
    print fl.time
