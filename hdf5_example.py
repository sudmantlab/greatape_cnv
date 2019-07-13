# hdf5 example
# Storing the CNV data of an individual in the following format:
# Each chromosome is its own group, e.g. '/chr1/'
# Within each chromosome group, the start and end positions are stored in a two-dimensional array: '/chr1/Window' -> [[0, 10], [10, 20], ...]
# Also within each chromosome group, the number of reads in each region are stored in a one-dimensional array: '/chr1/Depth' -> [50, 40, ...]


import h5py
import numpy as np


def get_depth():

	f = h5py.File("mytestfile.hdf5", "w")

	window_size = 10

	for chr in range(1, 3):
		chrom = 'chr' + str(chr)
		grp = f.create_group(chrom)
		print('GROUP NAME: ', grp.name)
		window = []
		count = []
		index = 0

		while index <= 20:
			ct = index
			window.append([index, index + window_size])
			count.append(ct)
			index += window_size

		window_dset = grp.create_dataset('Window', data = np.array(window))
		count_dset = grp.create_dataset('Depth', data = np.array(count))

		print('CHECK KEYS: ', list(f.keys()))
		print('WINDOW SHAPE: ', window_dset.shape)
		print('WHAT IS WINDOW? ', window_dset)


		print('COUNT SHAPE: ', count_dset.shape)
		print('WHAT IS COUNT? ', count_dset)


	f.close()



get_depth()