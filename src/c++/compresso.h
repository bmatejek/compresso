namespace compresso {
    unsigned char *Compress(unsigned long *data, long resolution[3], long steps[3]);

    unsigned long *Decompress(unsigned char *compressed_data);
};