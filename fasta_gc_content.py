import h5py
import numpy as np
import pysam
import pysamstats
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
from pyfasta import Fasta

"""
Calculates the percentage of GC content per sliding window within a reference fasta file
GC Percentage per window = (num_G + num_C)/(window_size - num_N)
Inputs:
	-fasta_file: input fasta file, e.g. 'panTro6.fa'
	-hdf5_file: name of the hdf5 output file, e.g. 'panTro6_gc.hdf5'
 	-window_size: size of the window used; normally 1000 bp
	-window_slide: the step size of the sliding window; normally 200 bp
Output:
An hdf5 file where each chromosome is its own group, e.g. '/chr1/'
Within each group, there are three datasets:
	-'Start': start index of the window
	-'End': end index of the window; should be start + window size except for possibly the last window
	-'GC Percent': GC percentage of the window, which equals (num_G + num_C)/(window_size - num_N)
"""



def ref_gc_hdf5(fasta_file, hdf5_file, window_size, window_slide):

	f_hdf5 = h5py.File(hdf5_file, "w")
	f_fasta = Fasta(fasta_file)

	for chr_name in f_fasta.keys():
		
		grp = f_hdf5.create_group(chr_name)

		chr_length = len(f_fasta[chr_name])

		start = np.array([], dtype='uint32')
		end = np.array([], dtype='uint32')
		count = np.array([], dtype='float32')

		index = 0
		while ((index + window_size) < chr_length):
			seq = f_fasta[chr_name][index:(index+window_size)]
			num_gc = 0
			num_N = 0
			for base in seq:
				if (base == 'G' or base == 'C'):
					num_gc += 1
				elif (base == 'N'):
					num_N += 1
			if (window_size == num_N):
				gc_percent = 0.0
			else:
				gc_percent = float(num_gc)/(len(seq) - num_N)
			count = np.append(count, gc_percent)
			start = np.append(start, index)
			end = np.append(end, index + window_size)
			index += window_slide

		start_dset = grp.create_dataset('Start', data=start)
		end_dset = grp.create_dataset('End', data=end)
		count_dset = grp.create_dataset('GC Percent', data=count)



	f_fasta.close()
	f_hdf5.close()

#ref_gc_hdf5('panTro6.fa', 'panTro6_gc.hdf5', 1000, 200)