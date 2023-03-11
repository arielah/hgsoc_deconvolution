# The data from the Vazquez-Garcia et al paper (GSE180661) has over a million
# cells and subsequently their HDF5 file is 40Gb in a compressed form and far
# too big to load into R at once on my machine, let alone run BayesPrism on it.
# After several dead ends, my workaround is to convert the HDF5 file to an mtx
# file, which is more compatible with the downstream packages I'm using, and
# which can be broken up into chunks. So this file reads in the HDF5 file and
# writes the data to a series of mtx files, 100,000 cells at a time.

import numpy as np
import h5py
import yaml
import os

with open('../../config.yml', 'r') as cfg:
    params = yaml.safe_load(cfg)


# Number of cells per chunk. R definitely can't handle 1.3 million at a time,
# it could probably do 200K at a time but I decided not to push it.
chunk_size = 100000

# Load data
infile = os.path.join(params["local_data_path"],
                      "Vazquez-Garcia",
                      "GSE180661_matrix.h5")

f = h5py.File(infile, 'r')

# The actual read count values
data = f['X']['data']  # Actual read count values

# The row (gene) number for each data point
indices = f['X']['indices']  # The row (gene) number for each data point

# The point in indices where the new column starts, e.g. "the next 5132
# values are in cell X". The indices values should wrap around there.
indptr = f['X']['indptr']

# The barcode labels for each cell/column
cells = np.array(f['obs'][:])

# The gene names for each row
genes = np.array(f['var'][:])

outdir = os.path.join(params["local_data_path"], "Vazquez-Garcia")

# A list of all mtx matrices in the chunk, because it's easy to
# append to and then concatenate at the end
mtx_list = []

# Iterate through the cells
for i in range(len(indptr)-1):

    # Format data in mtx format: row number, column number, value
    index_left = indptr[i]
    index_right = indptr[i+1]
    mini_index = indices[index_left:index_right]+1
    mini_data = data[index_left:index_right]
    mini_column = np.full(len(mini_data), i+1)
    mini_mtx = np.column_stack((mini_index, mini_column, mini_data))
    mtx_list.append(mini_mtx)
        
    # Write file every 100000 cells or at the end of the file
    if (i % chunk_size == 0 and i > 0) or i == len(indptr)-2:

        # Create directory for this chunk called "part_x"
        if i % chunk_size == 0:
            partno = int(i / chunk_size)
        elif i == len(indptr) - 2:
            partno = int(i / chunk_size) + 1
        partdir = "part_" + str(partno)
        chunkpath = os.path.join(outdir, partdir)
        os.makedirs(chunkpath, exist_ok=True)
        chunkfile = os.path.join(chunkpath, "matrix.mtx")
        
        # Concatenate all mtx content
        mtx_chunk = np.concatenate(mtx_list)

        with open(chunkfile, 'a') as o:
            # Write header so downstream package recognizes it as an mtx file
            # Coordinate shows it's listing row # and column #, not a 2D array
            # Integer means to expect all values to be integers (read counts)
            # General means there's not a symmetry constraint for the matrix
            o.write("%%MatrixMarket matrix coordinate integer general\n")

            # First line of mtx file is the total number of rows, cols, and
            # values to expect. It assumes a zero index, so we write the total
            # size of the final matrix but only populate it with the values
            # for some of the cells, so it has lots of empty cells. This is
            # resolved downstream.
            first_line = np.array([[len(genes), len(cells), len(mtx_chunk)]])
            np.savetxt(o, first_line, fmt='%i')
           
            # Save actual mtx content
            np.savetxt(o, mtx_chunk, fmt='%i')

        # Save a copy of the cell and barcode labels to the chunk directory
        featurefile = os.path.join(chunkpath, "features.tsv")
        np.savetxt(featurefile, genes, fmt="%s", delimiter="\t")
        cellfile = os.path.join(chunkpath, "barcodes.tsv")
        np.savetxt(cellfile, cells, fmt="%s")

        # Clear list to start the next file chunk
        mtx_list = []

print("Done")
