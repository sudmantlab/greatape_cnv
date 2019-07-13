# hdf5 example
# Storing the CNV data of an individual in the following format:
# Each chromosome is its own group, e.g. '/chr1/'
# Within each chromosome group, the start and end positions are stored in a two-dimensional array: '/chr1/Window' -> [[0, 10], [10, 20], ...]
# Also within each chromosome group, the number of reads in each region are stored in a one-dimensional array: '/chr1/Depth' -> [50, 40, ...]

import h5py
import numpy as np
import pysam


def get_depth(bam_file, window_size):
	
	samfile = pysam.AlignmentFile(bam_file, "rb")
	f = h5py.File("mytestfile.hdf5", "w")

	for chr in range(1, 24):
		chrom = 'chr' + str(chr)
		grp = f.create_group(chrom)
		window = []
		count = []
		index = 0
		while(True):
			try:
				ct = samfile.count(chrom, index, index + window_size) # is this the correct count?
				window.append([index, index + window_size])
				count.append(ct)
				index += window_size
			except ValueError:
				break

		window_dset = grp.create_dataset('Window', data=np.array(window))
		count_dset = grp.create_dataset('Depth', data=np.array(count))


	samfile.close()
	f.close()

# Try making an hdf5 file in Snakemake with each bam file, and then opening the hdf5 file and writing in it