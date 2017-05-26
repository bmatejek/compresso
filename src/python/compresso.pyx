cimport cython
cimport numpy as np
import numpy as np
import ctypes

# import c++ functions
cdef extern from "../c++/compresso.h" namespace "compresso":
    unsigned char *Compress(unsigned long *data, long res[3], long steps[3])
    unsigned long *Decompress(unsigned char *compressed_data)


class compresso(object):
    @staticmethod
    def name():
        return 'Compresso'

    @staticmethod
    def compress(data, steps):
        # get the shape 
        resolution = data.shape
        
        # convert data to c++ array
        nentries = data.size
        cdef np.ndarray[unsigned long, ndim=3, mode='c'] cpp_data
        cpp_data = np.ascontiguousarray(data, dtype=ctypes.c_uint64)
        cdef long *cpp_resolution = [resolution[0], resolution[1], resolution[2]]
        cdef long *cpp_steps = [steps[0], steps[1], steps[2]]
        cdef unsigned char *cpp_compressed_data = Compress(&(cpp_data[0,0,0]), cpp_resolution, cpp_steps)

        # convert the first 8 bytes to one unsigned long
        length = 0
        for i in range(8):
            length += (cpp_compressed_data[i]) * (2 ** (8 * i))

        # convert the c++ array to a numpy array
        cdef unsigned char[:] tmp_compressed_data = <unsigned char[:length]> cpp_compressed_data
        compressed_data = (np.asarray(tmp_compressed_data))[8:]
        
        return compressed_data

    @staticmethod
    def decompress(data):
        # get the resolution
        (zres, yres, xres) = (0, 0, 0)
        for i in range(8):
            zres += (data[i]) * (2 ** (8 * i))
            yres += (data[i + 8]) * (2 ** (8 * i))
            xres += (data[i + 16]) * (2 ** (8 * i))

        # convert the data to c++ array
        cdef np.ndarray[unsigned char, ndim=1, mode='c'] cpp_data
        cpp_data = np.ascontiguousarray(data, dtype=ctypes.c_uint8)
        cdef unsigned long *cpp_decompressed_data = Decompress(&(cpp_data[0]))
        cdef unsigned long[:] tmp_decompressed_data = <unsigned long[:zres*yres*xres]> cpp_decompressed_data
        decompressed_data = np.reshape(np.asarray(tmp_decompressed_data), (zres, yres, xres))

        return decompressed_data