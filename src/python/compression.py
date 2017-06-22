from compresso import Compresso
import h5py
import numpy as np
import bz2

with h5py.File('ac3_rhoana.h5', 'r') as hf:
    data = np.array(hf['main'])

data = data[0:1,0:512,0:512]

data = data.astype(np.uint64)
compressed_data = Compresso.compress(data, data.shape, (1, 8, 8))
decompressed_data = Compresso.decompress(compressed_data)
assert (np.array_equal(data, decompressed_data))

import struct
with open('output.b2c', 'wb') as fd:
   for byte in compressed_data:
       fd.write(struct.pack('B', byte))
