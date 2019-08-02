#bam -> hdf5
# Storing the average depth data of a sample in the following format:
# Each chromosome is its own group, e.g. '/chr1/'
#There are also groups for contigs that do not map to the reference genome ('ChrUn_...')
# Within each chromosome group, there are three arrays: 
#'Start' = the start index of each window, 0-indexed
#'End' = the end index of each window, end exclusive (ex. 0-200, includes 200)
#'Depth' = the total depth over the window divided by the size of the window

import h5py
import numpy as np
import pysam
import pysamstats
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def create(bam_file, hdf5_file, window_size, window_slide):

	bamfile = pysam.AlignmentFile(bam_file, "rb")
	f = h5py.File(hdf5_file, "w")
	bam_ref = bamfile.references
	bam_lengths = bamfile.lengths
	chrom_lengths = {i[0]:i[1] for i in zip(bam_ref, bam_lengths)}

	for item in chrom_lengths.items():

		chr_name = item[0]
		chr_length = item[1]
		grp = f.create_group(item[0])

		start = np.array([], dtype='uint32')
		end = np.array([], dtype='uint32')
		count = np.array([], dtype='float32')

		read_vals = pysamstats.load_coverage(bamfile, chrom=chr_name, start=0, end=chr_length)

		index = 0
		while ((index + window_size) < chr_length):
			window_reads = read_vals.reads_all[index:(index+window_size)]
			tot_reads = np.sum(window_reads)
			avg_depth = tot_reads/window_size
			count = np.append(count, avg_depth)
			start = np.append(start, index)
			end = np.append(end, index + window_size)
			index += window_slide

		start_dset = grp.create_dataset('Start', data=start)
		end_dset = grp.create_dataset('End', data=end)
		count_dset = grp.create_dataset('Depth', data=count)



	bamfile.close()
	f.close()


my_parser = argparse.ArgumentParser()
my_parser.add_argument('--bam_file', action='store', type=str, required=True)
my_parser.add_argument('--hdf5_file', action='store', type=str, required=True)
my_parser.add_argument('--window_size', action='store', type=int, required=True)
my_parser.add_argument('--window_slide', action='store', type=int, required=True)

args = vars(my_parser.parse_args())

create(args['bam_file'], args['hdf5_file'], args['window_size'], args['window_slide'])



#create('SRR726352.bam', 'SRR726352_chr1.hdf5', 1000, 200)


