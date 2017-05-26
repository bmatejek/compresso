#include <unordered_map>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include "compresso.h"


// dimension constants

static const short RN_X = 2;
static const short RN_Y = 1;
static const short RN_Z = 0;

static long row_size = -1;
static long sheet_size = -1;
static long grid_size = -1;



///////////////////////////////////
//// INTERNAL HELPER FUNCTIONS ////
///////////////////////////////////

static long
IndicesToIndex(long ix, long iy, long iz)
{
    return iz * sheet_size + iy * row_size + ix;
}



/////////////////////////////////////////
//// UNION-FIND CLASS FOR COMPONENTS ////
/////////////////////////////////////////

class UnionFindElement {
public:
    UnionFindElement(unsigned long label):
    label(label),
    parent(this),
    rank(0)
    {}

public:
    unsigned long label;
    UnionFindElement *parent;
    int rank;
};



UnionFindElement *
Find(UnionFindElement *x)
{
    if (x->parent != x) {
        x->parent = Find(x->parent);
    }
    return x->parent;
}



void
Union(UnionFindElement *x, UnionFindElement *y)
{
    UnionFindElement *xroot = Find(x);
    UnionFindElement *yroot = Find(y);

    if (xroot == yroot) return;

    // merge the two roots
    if (xroot->rank < yroot->rank) {
        xroot->parent = yroot;
    }
    else if (xroot->rank > yroot->rank) {
        yroot->parent = xroot;
    }
    else {
        yroot->parent = xroot;
        xroot->rank = xroot->rank + 1;
    }
}



///////////////////////////////
//// COMPRESSION ALGORITHM ////
///////////////////////////////

static bool *
ExtractBoundaries(unsigned long *data, long res[3])
{
    // create the boundaries array
    bool *boundaries = new bool[grid_size];
    if (!boundaries) { fprintf(stderr, "Failed to allocate memory for boundaries...\n"); return NULL; }

    // determine which pixels differ from east or south neighbors
    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                boundaries[iv] = false;

                // check the east neighbor
                if (ix < res[RN_X] - 1) {
                    if (data[iv] != data[IndicesToIndex(ix + 1, iy, iz)]) boundaries[iv] = true;
                }
                if (iy < res[RN_Y] - 1) {
                    if (data[iv] != data[IndicesToIndex(ix, iy + 1, iz)]) boundaries[iv] = true;
                }
            }
        }
    }

    // return the boundaries array
    return boundaries;
}



static unsigned long *
ConnectedComponents(bool *boundaries, long res[3])
{
    // create the connected components grid
    unsigned long *components = new unsigned long[grid_size];
    if (!components) { fprintf(stderr, "Failed to allocate memory for connected components...\n"); return NULL; }

    // initialize to zero
    for (long iv = 0; iv < grid_size; ++iv)
        components[iv] = 0;

    // run connected components for every slice
    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        // create a vector of union find elements
        std::vector<UnionFindElement *> union_find = std::vector<UnionFindElement *>();

        // current label in connected component
        int curlab = 1;
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                // continue if boundary
                if (boundaries[iv]) continue;

                // only consider the pixel to the north and west
                long north = IndicesToIndex(ix - 1, iy, iz);
                long west = IndicesToIndex(ix, iy - 1, iz);

                unsigned long neighbor_labels[2] = { 0, 0 };

                // get the labels for the relevant neighbor
                if (ix > 0) neighbor_labels[0] = components[north];
                if (iy > 0) neighbor_labels[1] = components[west];

                // if the neighbors are boundary, create new label
                if (!neighbor_labels[0] && !neighbor_labels[1]) {
                    components[iv] = curlab;

                    // add to union find structure
                    union_find.push_back(new UnionFindElement(0));

                    // update the next label
                    curlab++;
                }
                // the two pixels have equal non-trivial values
                else if (neighbor_labels[0] == neighbor_labels[1])
                    components[iv] = neighbor_labels[0];
                else {
                    if (!neighbor_labels[0]) components[iv] = neighbor_labels[1];
                    else if (!neighbor_labels[1]) components[iv] = neighbor_labels[0];
                    else {
                        // take the minimum value
                        components[iv] = std::min(neighbor_labels[0], neighbor_labels[1]);

                        // set the equivalence relationship
                        Union(union_find[neighbor_labels[0] - 1], union_find[neighbor_labels[1] - 1]);
                    }
                }
            }
        }

        // reset the current label to 1
        curlab = 1;

        // create the connected components in order
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                if (boundaries[iv]) continue;

                // get the parent for this component
                UnionFindElement *comp = Find(union_find[components[iv] - 1]);
                if (!comp->label) {
                    comp->label = curlab;
                    curlab++;
                }

                components[iv] = comp->label;
            }
        }

        // free memory
        for (unsigned long iv = 0; iv < union_find.size(); ++iv)
            delete union_find[iv];
    }

    // return the connected components array
    return components;
}



static std::vector<unsigned long> *
IDMapping(unsigned long *components, unsigned long *data, long res[3])
{
    // create a vector of the ids
    std::vector<unsigned long> *ids = new std::vector<unsigned long>();

    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        // create a set for this individual slice
        std::unordered_set<unsigned long> hash_map = std::unordered_set<unsigned long>();

        // iterate over the entire slice
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                // get the component label
                unsigned long component_id = components[iv];

                // if this component does not belong yet, add it
                if (!hash_map.count(component_id)) {
                    hash_map.insert(component_id);

                    // add the segment id
                    unsigned long segment_id = data[iv] + 1;
                    ids->push_back(segment_id);
                }
            }
        }
    }

    // return the mapping
    return ids;
}



static unsigned long *
EncodeBoundaries(bool *boundaries, long res[3], long steps[3])
{
    // determine the number of bloxks in each direction
    long nblocks[3];
    for (int dim = 0; dim <= 2; ++dim) {
        nblocks[dim] = (long) (ceil((double)res[dim] / steps[dim]) + 0.5);
    }
    long nwindows = nblocks[RN_Z] * nblocks[RN_Y] * nblocks[RN_X];

    // create an empty array for the encodings
    unsigned long *boundary_data = new unsigned long[nwindows];
    if (!boundary_data) { fprintf(stderr, "Failed to allocate memory for boundary windows...\n"); return NULL; }
    for (long iv = 0; iv < nwindows; ++iv)
        boundary_data[iv] = 0;

    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                // no encoding for non-boundaries
                if (!boundaries[iv]) continue;

                // find the block from the index
                long zblock = iz / steps[RN_Z];
                long yblock = iy / steps[RN_Y];
                long xblock = ix / steps[RN_X];

                // find the offset within the block
                long zoffset = iz % steps[RN_Z];
                long yoffset = iy % steps[RN_Y];
                long xoffset = ix % steps[RN_X];

                long block = zblock * (nblocks[RN_Y] * nblocks[RN_X]) + yblock * nblocks[RN_X] + xblock;
                long offset = zoffset * (steps[RN_Y] * steps[RN_X]) + yoffset * steps[RN_X] + xoffset;

                boundary_data[block] += (1LU << offset);
            }
        }
    }

    // return the encodings
    return boundary_data;
}



static std::vector<unsigned long> *
ValueMapping(unsigned long *boundary_data, long nwindows)
{
    // get a list of unique values
    std::vector<unsigned long> *values = new std::vector<unsigned long>();
    std::unordered_set<unsigned long> hash_map = std::unordered_set<unsigned long>();

    // go through all of the boundary data to create array of values
    for (long iv = 0; iv < nwindows; ++iv) {
        if (!hash_map.count(boundary_data[iv])) {
            hash_map.insert(boundary_data[iv]);
            values->push_back(boundary_data[iv]);
        }
    }

    // sort the values
    sort(values->begin(), values->end());

    // create mapping from values to indices
    std::unordered_map<unsigned long, unsigned long> mapping = std::unordered_map<unsigned long, unsigned long>();
    for (unsigned long iv = 0; iv < values->size(); ++iv) {
        mapping[(*values)[iv]] = iv;
    }

    // update boundary data
    for (long iv = 0; iv < nwindows; ++iv) {
        boundary_data[iv] = mapping[boundary_data[iv]];
    }

    // return the values
    return values;
}



static std::vector<unsigned long> *
EncodeIndeterminateLocations(bool *boundaries, unsigned long *data, long res[3])
{
    std::vector<unsigned long> *locations = new std::vector<unsigned long>();

    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                if (!boundaries[iv]) continue;
                else if (iy > 0 && !boundaries[IndicesToIndex(ix, iy - 1, iz)]) continue;
                else if (ix > 0 && !boundaries[IndicesToIndex(ix - 1, iy, iz)]) continue;
                else {
                    long north = IndicesToIndex(ix - 1, iy, iz);
                    long south = IndicesToIndex(ix + 1, iy, iz);
                    long east = IndicesToIndex(ix, iy - 1, iz);
                    long west = IndicesToIndex(ix, iy + 1, iz);
                    long up = IndicesToIndex(ix, iy, iz + 1);
                    long down = IndicesToIndex(ix, iy, iz - 1);

                    // see if any of the immediate neighbors are candidates
                    if (ix > 0 && !boundaries[north] && data[north] == data[iv]) locations->push_back(0);
                    else if (ix < res[RN_X] - 1 && !boundaries[south] && data[south] == data[iv]) locations->push_back(1);
                    else if (iy > 0 && !boundaries[east] && data[east] == data[iv]) locations->push_back(2);
                    else if (iy < res[RN_Y] - 1 && !boundaries[west] && data[west] == data[iv]) locations->push_back(3);
                    else if (iz > 0 && !boundaries[down] && data[down] == data[iv]) locations->push_back(4);
                    else if (iz < res[RN_Z] - 1 && !boundaries[up] && data[up] == data[iv]) locations->push_back(5);
                    else locations->push_back(data[IndicesToIndex(ix, iy, iz)] + 6);
                }
            }
        }
    }

    return locations;
}



static unsigned long
WriteUint32(unsigned char *data, unsigned long offset, unsigned long value)
{
    data[offset + 3] = (value >> 24);
    data[offset + 2] = (value >> 16) % 256;
    data[offset + 1] = (value >> 8) % 256;
    data[offset] = (value) % 256;

    // return the new offset
    return offset + 4;
}



static unsigned long
WriteUint64(unsigned char *data, unsigned long offset, unsigned long value)
{
    data[offset + 7] = (value >> 56);
    data[offset + 6] = (value >> 48) % 256;
    data[offset + 5] = (value >> 40) % 256;
    data[offset + 4] = (value >> 32) % 256;
    data[offset + 3] = (value >> 24) % 256;
    data[offset + 2] = (value >> 16) % 256;
    data[offset + 1] = (value >> 8) % 256;
    data[offset] = value % 256;

    // return the new offset
    return offset + 8;
}



static unsigned int
ReadUint32(unsigned char *data, unsigned long &offset)
{
    unsigned int value = 0;
    for (int iv = 0; iv < 4; ++iv) 
        value += data[offset + iv] << (8 * iv);

    offset += 4;

    return value;
}



static unsigned long
ReadUint64(unsigned char *data, unsigned long &offset)
{
    unsigned long value = 0;
    for (int iv = 0; iv < 8; ++iv)
        value += ((unsigned long)data[offset + iv]) << (8 * iv);

    offset += 8;

    return value;
}




unsigned char *
compresso::Compress(unsigned long *data, long res[3], long steps[3])
{
    row_size = res[RN_X];
    sheet_size = res[RN_X] * res[RN_Y];
    grid_size = res[RN_X] * res[RN_Y] * res[RN_Z];

    // determine the number of bloxks in each direction
    long nblocks[3];
    for (int dim = 0; dim <= 2; ++dim) {
        nblocks[dim] = (long) (ceil((double)res[dim] / steps[dim]) + 0.5);
    }
    long nwindows = nblocks[RN_Z] * nblocks[RN_Y] * nblocks[RN_X];

    // get the boundary voxels
    bool *boundaries = ExtractBoundaries(data, res);
    if (!boundaries) return NULL;

    // get the connected components
    unsigned long *components = ConnectedComponents(boundaries, res);
    if (!components) return NULL;

    // get the ids
    std::vector<unsigned long> *ids = IDMapping(components, data, res);
    if (!ids) return NULL;

    // encode the boundary data
    unsigned long *boundary_data = EncodeBoundaries(boundaries, res, steps);
    if (!boundary_data) return NULL;

    // get the values from the boundary data
    std::vector<unsigned long> *values = ValueMapping(boundary_data, nwindows);
    if (!values) return NULL;

    // get the locations
    std::vector<unsigned long> *locations = EncodeIndeterminateLocations(boundaries, data, res);
    if (!locations) return NULL;

    // create an array of bytes
    unsigned short header_size = 11;
    unsigned long ids_size = ids->size();
    unsigned long values_size = values->size();
    unsigned long locations_size = locations->size();
    /* TODO fix this */
    unsigned long byte_size = 4;

    // count the number of output uint64 bytes
    unsigned long nuint64_bytes = (header_size + ids_size + values_size + locations_size) * 8;
    unsigned long nuint32_bytes = nwindows * 4;

    // get the total number of bytes
    unsigned long nbytes = nuint64_bytes + nuint32_bytes;
    
    // create the compressed data output
    unsigned char *compressed_data = new unsigned char[nbytes];
    if (!compressed_data) return NULL;

    // counter for the byte offset
    unsigned long offset = 0;

    // write the header
    offset = WriteUint64(compressed_data, offset, nbytes);
    offset = WriteUint64(compressed_data, offset, res[RN_Z]);
    offset = WriteUint64(compressed_data, offset, res[RN_Y]);
    offset = WriteUint64(compressed_data, offset, res[RN_X]);
    offset = WriteUint64(compressed_data, offset, ids_size);
    offset = WriteUint64(compressed_data, offset, values_size);
    offset = WriteUint64(compressed_data, offset, byte_size);
    offset = WriteUint64(compressed_data, offset, locations_size);
    offset = WriteUint64(compressed_data, offset, steps[RN_Z]);
    offset = WriteUint64(compressed_data, offset, steps[RN_Y]);
    offset = WriteUint64(compressed_data, offset, steps[RN_X]);

    unsigned int *contracted_data = new unsigned int[nwindows];
    if (!contracted_data) return NULL;
    for (long iv = 0; iv < nwindows; ++iv)
        contracted_data[iv] = (unsigned int)boundary_data[iv];    

    // write the id values
    for (unsigned long iv = 0; iv < ids_size; ++iv)
        offset = WriteUint64(compressed_data, offset, (*ids)[iv]);
    for (unsigned long iv = 0; iv < values_size; ++iv)
        offset = WriteUint64(compressed_data, offset, (*values)[iv]);
    for (unsigned long iv = 0; iv < locations_size; ++iv)
        offset = WriteUint64(compressed_data, offset, (*locations)[iv]);
    for (long iv = 0; iv < nwindows; ++iv)
        offset = WriteUint32(compressed_data, offset, contracted_data[iv]);

    delete[] contracted_data;
    delete[] boundaries;
    delete[] components;
    delete values;
    delete locations;
    delete ids;
    delete[] boundary_data;

    return compressed_data;
}



/////////////////////////////////
//// DECOMPRESSION ALGORITHM ////
/////////////////////////////////

static bool *
DecodeBoundaries(unsigned long *boundary_data, std::vector<unsigned long> *values, long res[3], long steps[3])
{
    // determine the number of bloxks in each direction
    long nblocks[3];
    for (int dim = 0; dim <= 2; ++dim) {
        nblocks[dim] = (long) (ceil((double)res[dim] / steps[dim]) + 0.5);
    }

    bool *boundaries = new bool[grid_size];
    for (long iv = 0; iv < grid_size; ++iv)
        boundaries[iv] = false;

    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                // find the block from the index
                long zblock = iz / steps[RN_Z];
                long yblock = iy / steps[RN_Y];
                long xblock = ix / steps[RN_X];

                // find the offset within the block
                long zoffset = iz % steps[RN_Z];
                long yoffset = iy % steps[RN_Y];
                long xoffset = ix % steps[RN_X];

                long block = zblock * (nblocks[RN_Y] * nblocks[RN_X]) + yblock * nblocks[RN_X] + xblock;
                long offset = zoffset * (steps[RN_Y] * steps[RN_X]) + yoffset * steps[RN_X] + xoffset;

                unsigned long value = (*values)[boundary_data[block]];
                if ((value >> offset) % 2) boundaries[iv] = true;
            }
        }
    }

    return boundaries;
}



static unsigned long *
IDReverseMapping(unsigned long *components, std::vector<unsigned long> *ids, long res[3])
{
    unsigned long *decompressed_data = new unsigned long[grid_size];
    for (long iv = 0; iv < grid_size; ++iv)
        decompressed_data[iv] = 0;

    int ids_index = 0;
    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        // create mapping (not memory efficient but FAST!!)
        // number of components is guaranteed to be less than ids->size()
        unsigned long *mapping = new unsigned long[ids->size()];
        for (unsigned long iv = 0; iv < ids->size(); ++iv)
            mapping[iv] = 0;

        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);

                if (!mapping[components[iv]]) {
                    mapping[components[iv]] = (*ids)[ids_index];
                    ids_index++;
                }

                decompressed_data[iv] = mapping[components[iv]] - 1;
            }
        }
    }

    return decompressed_data;
}



static void
DecodeIndeterminateLocations(bool *boundaries, unsigned long *decompressed_data, std::vector<unsigned long> *locations, long res[3])
{
    long index = 0;

    // go through all voxels
    for (long iz = 0; iz < res[RN_Z]; ++iz) {
        for (long iy = 0; iy < res[RN_Y]; ++iy) {
            for (long ix = 0; ix < res[RN_X]; ++ix) {
                long iv = IndicesToIndex(ix, iy, iz);
                
                // get the north and west neighbors
                long north = IndicesToIndex(ix - 1, iy, iz);
                long west = IndicesToIndex(ix, iy - 1, iz);

                if (!boundaries[iv]) continue;
                else if (ix > 0 && !boundaries[north]) {
                    decompressed_data[iv] = decompressed_data[north];
                }
                else if (iy > 0 && !boundaries[west]) {
                    decompressed_data[iv] = decompressed_data[west];
                }
                else {
                    unsigned long offset = (*locations)[index];
                    if (offset == 0) decompressed_data[iv] = decompressed_data[IndicesToIndex(ix - 1, iy, iz)];
                    else if (offset == 1) decompressed_data[iv] = decompressed_data[IndicesToIndex(ix + 1, iy, iz)];
                    else if (offset == 2) decompressed_data[iv] = decompressed_data[IndicesToIndex(ix, iy - 1, iz)];
                    else if (offset == 3) decompressed_data[iv] = decompressed_data[IndicesToIndex(ix, iy + 1, iz)];
                    else if (offset == 4) decompressed_data[iv] = decompressed_data[IndicesToIndex(ix, iy, iz - 1)];
                    else if (offset == 5) decompressed_data[iv] = decompressed_data[IndicesToIndex(ix, iy, iz + 1)];
                    else decompressed_data[iv] = offset - 6;
                    index++;
                }
            }
        }
    }
}



unsigned long *
compresso::Decompress(unsigned char *compressed_data)
{
    // extract the header
    unsigned long offset = 0;
    long res[3];
    long steps[3];
    res[RN_Z] = ReadUint64(compressed_data, offset);
    res[RN_Y] = ReadUint64(compressed_data, offset);
    res[RN_X] = ReadUint64(compressed_data, offset);
    unsigned long ids_size = ReadUint64(compressed_data, offset);
    unsigned long values_size = ReadUint64(compressed_data, offset);
    unsigned long byte_size = ReadUint64(compressed_data, offset);
    unsigned long locations_size = ReadUint64(compressed_data, offset);
    steps[RN_Z] = ReadUint64(compressed_data, offset);
    steps[RN_Y] = ReadUint64(compressed_data, offset);
    steps[RN_X] = ReadUint64(compressed_data, offset);

    // set global variables
    row_size = res[RN_X];
    sheet_size = res[RN_Y] * res[RN_X];
    grid_size = res[RN_Z] * res[RN_Y] * res[RN_X];

    // determine the number of blocks in the z, y, and x dimensions
    long nblocks[3];
    for (int dim = 0; dim <= 2; ++dim) {
        nblocks[dim] = (long) (ceil((double)res[dim] / steps[dim]) + 0.5);
    }
    long nwindows = nblocks[RN_Z] * nblocks[RN_Y] * nblocks[RN_X];

    // allocate memory for all arrays
    std::vector<unsigned long> *ids = new std::vector<unsigned long>();
    if (!ids) return NULL;
    std::vector<unsigned long> *values = new std::vector<unsigned long>();
    if (!values) return NULL;
    std::vector<unsigned long> *locations = new std::vector<unsigned long>();
    if (!locations) return NULL;
    unsigned long *boundary_data = new unsigned long[nwindows];    

    for (unsigned long iv = 0; iv < ids_size; ++iv)
        ids->push_back(ReadUint64(compressed_data, offset));
    for (unsigned long iv = 0; iv < values_size; ++iv)
        values->push_back(ReadUint64(compressed_data, offset));
    for (unsigned long iv = 0; iv < locations_size; ++iv)
        locations->push_back(ReadUint64(compressed_data, offset));
    /* TODO allow different data types */
    for (long iv = 0; iv < nwindows; ++iv)
        boundary_data[iv] = ReadUint32(compressed_data, offset);

    // get the boundaries from the data
    bool *boundaries = DecodeBoundaries(boundary_data, values, res, steps);
    if (!boundaries) return NULL;

    // get the connected components
    unsigned long *components = ConnectedComponents(boundaries, res);
    if (!components) return NULL;

    // decompress the data
    unsigned long *decompressed_data = IDReverseMapping(components, ids, res);
    if (!decompressed_data) return NULL;

    // decode the final indeterminate locations
    DecodeIndeterminateLocations(boundaries, decompressed_data, locations, res);

    // free memory
    delete[] boundaries;
    delete[] components;
    delete[] boundary_data;
    delete ids;
    delete values;
    delete locations;

    return decompressed_data;

}