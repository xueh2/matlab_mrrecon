function im = load_GT_UT_ref_h5(filename, tag)
im = h5read(filename, tag);