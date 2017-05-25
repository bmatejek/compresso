cimport cython
cimport numpy as np
import numpy as np
import ctypes
import os.path
import h5py

# import c++ functions
cdef extern from "../c++/compresso.h" namespace "compresso":
    unsigned char *Compress(unsigned long *data, long res[3], long steps[3])
    unsigned long *Decompress(unsigned long *compressed_data)


class compresso(object):
    @staticmethod
    def name():
        return 'Compresso'

    @staticmethod
    def Compress(data, steps):
        # get the shape 
        resolution = data.shape
        
        # convert data to c++ array
        nentries = data.size
        cdef np.ndarray[unsigned long, ndim=3, mode='c'] cpp_data
        cpp_data = np.ascontiguousarray(data, dtype=ctypes.c_uint64)
        cdef long *cpp_resolution = [resolution[0], resolution[1], resolution[2]]
        cdef long *cpp_steps = [steps[0], steps[1], steps[2]]
        cdef unsigned char *cpp_compressed_data = Compress(&(cpp_data[0,0,0]), cpp_resolution, cpp_steps)

        length = 0
        for i in range(8):
            length += ord(cpp_compressed_data[i]) * (2 ** (8 * i))

        cdef unsigned char[:] tmp_compressed_data = <unsigned char[:length]> cpp_compressed_data
        compressed_data = (np.asarray(tmp_compressed_data))[8:]

        return compressed_data

    @staticmethod
    def Decompress(data):
        return 0