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

        # read header files
        header = open('_contents.0').readline().split()
        self.nrec = int(header[0])

        header = open('nxnz.0').readline().split()
        self.nx = int(header[0]) + 1
        self.nz = int(header[1]) + 1

        self.nodes = self.nx * self.nz
        self.elements = (self.nx - 1) * (self.nz - 1)

        return


    def read_mesh(self, frame):
        columns = 2
        f = open('mesh.0')
        offset = frame * columns * self.nodes * sizeoffloat
        f.seek(offset)
        self.z, self.x = self._read_data(f, columns)
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
    This class is preferred over class Flac.
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
