cimport cython
cimport numpy as np
import numpy as np
import ctypes
import os.path
import h5py

# import c++ functions
cdef extern from "../c++/compresso.h" namespace "compresso":
    unsigned long *Compress(unsigned long *data, long res[3], long steps[3])
    unsigned long *Decompress(unsigned long *compressed_data)


def ReadH5File(filename, dataset):
    with h5py.File(filename, 'r') as hf:
        data = np.array(hf[dataset])

    return data


class compresso(object):
    @staticmethod
    def name():
        return 'Compresso'

    @staticmethod
    def Compress(filename, steps):
        # get the filename extension 
        _, extension = os.path.splitext(filename)

        if extension == '.h5': data = ReadH5File(filename, 'main')

        data = data[0:10,:,:]
        
        # get the shape 
        resolution = data.shape
        
        # convert data to c++ array
        nentries = data.size
        cdef np.ndarray[unsigned long, ndim=3, mode='c'] cpp_data
        cpp_data = np.ascontiguousarray(data, dtype=ctypes.c_uint64)
        cdef long *cpp_resolution = [resolution[0], resolution[1], resolution[2]]
        cdef long *cpp_steps = [steps[0], steps[1], steps[2]]
        cdef unsigned long *cpp_compressed_data = Compress(&(cpp_data[0,0,0]), cpp_resolution, cpp_steps)
        zres, yres, xres = data.shape
        import struct
        from PIL import Image
        # read the boundaries
        with open('boundaries', 'rb') as fp:
            boundaries = np.zeros((zres, yres, xres), dtype=np.uint8)

            for iz in range(zres):
                for iy in range(yres):
                    for ix in range(xres):
                        boundaries[iz,iy,ix], = struct.unpack('B', fp.read(1))

            image = Image.fromarray(255 * boundaries[2,:,:])
            image.save('boundaries.png')
        with open('components', 'rb') as fp:
            components = np.zeros((zres, yres, xres), dtype=np.uint64)

            for iz in range(zres):
                for iy in range(yres):
                    for ix in range(xres):
                        components[iz,iy,ix], = struct.unpack('Q', fp.read(8))

            components_rgb = np.zeros((yres, xres, 3), dtype=np.uint8)
            for iy in range(yres):
                for ix in range(xres):
                    red = ((np.uint64(107 * components[0,iy,ix])) % 700) % 255
                    green = ((np.uint64(509 * components[0,iy,ix])) % 900) % 255
                    blue = ((np.uint64(200 * components[0,iy,ix])) % 777) % 255

                    components_rgb[iy,ix,0] = red
                    components_rgb[iy,ix,1] = green
                    components_rgb[iy,ix,2] = blue

            image = Image.fromarray(components_rgb)
            image.save('components.png')
        return 0

    @staticmethod
    def Decompress(data):
        return 0