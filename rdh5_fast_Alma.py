import pysam
import pysamstats
import tables
import numpy as np
import argparse
import time
import math
import h5py
import os
from scipy.ndimage.interpolation import shift
import seaborn as sns
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from .hdf5_plotting import *

"""
python -u rdh5_fast.py create --fn_bam SRR748186.bam --fn_out SRR748186_fast.h5 --stats
Note: start and end indices of each window are inclusive
"""


# todo: write better documentation for all of this
def create_h5(args):
	fn_bam = args.fn_bam
	window = args.window
	step = args.step
	fn_out = args.fn_out
	chrs = args.chrs
	bam_index = args.bam_index
	stats = args.stats
	ref = args.ref
	# copynum = args.copynum

	assert window > step, "Window must be greater than step"
	assert window > 0, "Window cannot equal zero"
	assert step > 0, "Step cannot equal zero"

	# if fn_out wasn't specified, name it the same as fn_bam
	if fn_out is None:
		fn_out = fn_bam[:len(fn_bam) - 4] + ".h5"
	else:
		assert fn_out[len(fn_out) - 3:] == ".h5", "Include \'.h5\' extension at the end of fn_out"

	# for generating stats only when HDF5 file exists
	if stats and os.path.exists(fn_out):
		gen_stats(fn_bam, window, step, fn_out, bam_index)

	# for generating copy number only when HDF5 file exists
	# if copynum and os.path.exists(fn_out):
	# 	gen_copynum(fn_out)

	# finding index file
	if bam_index is None:
		bam_file = pysam.AlignmentFile(fn_bam, "rb") # will open index file automatically if in the same folder
	else:
		bam_file = pysam.AlignmentFile(fn_bam, "rb", index_filename=bam_index) # will open index file automatically if in the same folder

	refs = bam_file.references
	contig_lengths = {x[0]:x[1] for x in zip(refs, bam_file.lengths)}

	# checking whether chromosomes were specified
	if len(chrs) != 0:
		temp = {}
		for c in chrs:
			temp[c] = contig_lengths[c]
		contig_lengths = temp

	out_file = tables.open_file(fn_out, mode='w')
	depth_group = out_file.create_group(out_file.root, "depth")

	for r in contig_lengths.items():
		print("--------------------------------------------------")
		chrom_name = r[0]
		chrom_length = r[1]
		print("Chromosome name:", chrom_name)
		print("chromosome length:", chrom_length)
		chrom_group = out_file.create_group(depth_group, chrom_name)

		# Generate array of read-depth per base using pysamstats
		# cvg_array.dytpe = dtype((numpy.record, [('chrom', 'S27'), ('pos', '<i4'), ('reads_all', '<i4'), ('reads_pp', '<i4')]))
		# A numpy.recarray with each column with following format: (name, dtype)
		# 'S27' = <27-byte string, '<i4' = little-endian 32-bit int
		# chrom = reference/chromosome name
		# pos = base position in reference/chromosome
		# reads_all = Number of reads aligned at the position. N.b., this is really the total, i.e., includes reads where the mate is unmapped or otherwise not properly paired.
		# reads_pp = Number of reads flagged as properly paired by the aligner
		t = time.time()
		cvg_array = pysamstats.load_coverage(bam_file, chrom=chrom_name, start=0, end=chrom_length, pad=True)
		cvg = cvg_array.reads_all
		print("Time to load coverage array for", chrom_name, ":", (time.time()-t) / 60, "minutes")

		size = math.ceil(chrom_length / step - window / step + 1) # number of windows in chrom
		if chrom_length < window:
			size = 1
		print("Number of windows in chrom {}: {}".format(chrom_name, size))

		sums = np.zeros(size, dtype=np.uint32) # an array to hold window sums
		read_depth = np.zeros(size, dtype=np.float16) # array that holds final average read depths
		start = np.zeros(size, dtype=np.uint32) # array that holds start index of each window
		end = np.zeros(size, dtype=np.uint32) # array that holds final index of each window

		cvg_cum_sum = np.cumsum(cvg) # calculate cumulative sum of coverage array
		cvg_cum_sum_shift = shift(cvg_cum_sum, window, cval=0) # shift cumsum array to the right by WINDOW spaces, filling in a zero in each space
		cvg_final = cvg_cum_sum - cvg_cum_sum_shift # subtract the two arrays to get the window sums for step of 1

		# extract window sums at each step of STEP
		i = window - 1
		j = 0
		while i < cvg_final.size:
			sums[j] = cvg_final[i]
			start[j] = i - window + 1
			end[j] = i
			i += step
			j += 1

		# divide the window sums by WINDOW to get the mean
		np.divide(sums, window, out=read_depth)

		# tail case
		if i - step < cvg_final.size - 1:
			read_depth[size - 1] = np.mean(cvg[i - window + 1 : chrom_length])
			start[size - 1] = i - window + 1
			end[size - 1] = chrom_length - 1

		# save arrays in hdf5
		out_file.create_carray(chrom_group, name="read_depth", atom=tables.Float16Atom(), obj=read_depth)
		out_file.create_carray(chrom_group, name="start", atom=tables.UInt32Atom(), obj=start)
		out_file.create_carray(chrom_group, name="end", atom=tables.UInt32Atom(), obj=end)

		print("--------------------------------------------------")


	bam_file.close()
	out_file.close()

	#Generate stats, such as mean, median, and max:
	gen_stats(fn_bam, window, step, fn_out, bam_index)

	#Generate read depth percentiles:
	ref_gc = ref[:len(ref) - 3] + '_gc.hdf5'
	#First generate GC percentile per window for reference genome if DNE:
	if not os.path.exists(fn_out):
		ref_gc_hdf5(ref, ref_gc, window, step)

	gc_percentiles(fn_bam, ref_gc)

	#Generate GC-corrected read depth and copy number:
	gc_correction(sample_hdf5, sample_gc_hdf5, ref_hdf5)

	# generate stats
	if stats:
		gen_stats(fn_bam, window, step, fn_out, bam_index)

	# generate copynum
	# if copynum:
	# 	gen_copynum(fn_out)

def gen_stats(bam, window, step, fn_out, idx):
	f = h5py.File(fn_out, mode='r+') # read/write mode iff file exists
	stats_group = f.create_group("stats")
	depth_group = f["depth"]
	stats = {"window": window, "step": step, "BAM_path": bam}

	if idx is None and os.path.isfile(bam + ".bai"):
		stats["BAM_index_path"] = bam + ".bai"
	else:
		stats["BAM_index_path"] = idx

	# calculate mean, median, and max for each chromosome's read_depth array
	chrs = list(depth_group)
	chrs_len = len(chrs)
	means = np.zeros(chrs_len, dtype=np.float16)
	medians = np.zeros(chrs_len, dtype=np.float16)
	maxes = np.zeros(chrs_len, dtype=np.float16)
	for i in range(chrs_len):
		rd_arr = depth_group[chrs[i]]["read_depth"][:]
		means[i] = np.mean(rd_arr)
		medians[i] = np.median(rd_arr)
		maxes[i] = np.max(rd_arr)

	# calculate mean, median, and max for entire genome
	stats["mean"] = np.mean(means)
	stats["median"] = np.median(medians)
	stats["max"] = np.max(maxes)

	stats_group.create_dataset('stats', data=str(stats))

	f.close()

# different implementation of get_stats where stats are saved as groups instead of in a json
# todo: working progress
def gen_stats2(bam, window, step, fn_out, idx):
	f = h5py.File(fn_out, mode='r+') # read/write mode iff file exists
	stats_group = f.create_group("stats")
	depth_group = f.get("depth")
	stats = {"window": window, "step": step, "BAM_path": bam}

	if idx is None and os.path.isfile(bam + ".bai"):
		stats["BAM_index_path"] = bam + ".bai"
	else:
		stats["BAM_index_path"] = idx

	# calculate mean, median, and max for each chromosome's read_depth array
	chrs = list(depth_group.__iter__())
	chrs_len = len(chrs)
	means = np.zeros(chrs_len, dtype=np.float16)
	medians = np.zeros(chrs_len, dtype=np.float16)
	maxes = np.zeros(chrs_len, dtype=np.int32)
	for i in range(chrs_len):
		rd_arr = depth_group[chrs[i]]["read_depth"][:]
		means[i] = np.mean(rd_arr)
		medians[i] = np.median(rd_arr)
		maxes[i] = np.max(rd_arr)

	# calculate mean, median, and max for entire genome
	stats["mean"] = np.mean(means)
	stats["median"] = np.median(medians)
	stats["max"] = np.max(maxes)

	for i in stats.items():
		stats_group.create_dataset(np.array(i[0]), data=np.array(i[1])) # fixme

	f.close()

def get_stats(hdf5):
	f = h5py.File(hdf5, 'r')
	stats_dict = f.get("stats").get("stats")
	if stats_dict is None:
		print("Stats have not been generated. Rerun with create and --gen_stats flags.")
		exit(1)
	stats_dict = eval(stats_dict[...][()]) # https://docs.scipy.org/doc/numpy-1.14.1/reference/arrays.scalars.html
	return stats_dict

# todo: working progress
def get_stats2(hdf5):
	f = h5py.File(hdf5, 'r')
	stats_list = list(f['stats'])
	stats_dict = {}
	for s in stats_list:
		stats_dict[s] = stats_list[s][...][()]
	return stats_dict

# # todo: test this
# # todo: save copynums in HDF5 file?
# def gen_copynum(hdf5):
# 	stats_dict = get_stats(hdf5)
# 	mean = stats_dict["mean"] # todo: mean is GC-corrected? read_depth arrays are GC-corrected?
# 	scalar = 2 / mean # assuming mean cvg is equivalent to copy number of 2, copynum = cvg * 2 / mean
# 	file = tables.open_file(hdf5)
# 	chrs = list(file.get_node(where=file.root.depth))
# 	for node in chrs:
# 		rd_arr = node.read_depth[:]
# 		cn_arr = np.multiply(rd_arr, scalar)
# 		file.create_carray(where=node, name="copynum", atom=tables.Float16Atom(), data=cn_arr)


def summary_h5(args):
	fn_h5 = args.fn_h5
	stats = args.stats
	if stats:
		print(get_stats(fn_h5))



"""
The function gc_percentiles:
Calculates the average read depth of windows for every GC pertentile
Input: 
The sample hdf5 file and the reference genome hdf5 file containing the GC percentage per window (calculated with the function ref_gc_hdf5)
Output: 
Adds a group 'GC' to the sample hdf5 file containing four datasets:
The dataset 'GC Percentiles':
	-contains the GC integer percentiles from 0 to 100 (101 values)
	-for each percentile p, the percentile bin is the rounded GC percent to the nearest integer
The dataset 'Average Depth':
	-contains the average of the read depths of all windows across the entire genome per GC percentile
The dataset 'Num Windows':
	-contains the number of windows across the entire genome per GC percentile
The dataset 'Poly Fit':
	-contains 101 values of a degree-3 polynomial fit to all the read depths vs GC percents across the entire genome
"""


def gc_percentiles(sample_hdf5, ref_hdf5):

	sample_file = h5py.File(sample_hdf5, 'r')
	ref_file = h5py.File(ref_hdf5, 'r')

	gc_content = {val:[] for val in np.arange(0,101)}
	num_windows = {val:[0] for val in np.arange(0,101)}

	sample_chrs = sample_file['depth']

	total_gc_percents = np.array([], dtype='float16')
	total_read_depths = np.array([], dtype='uint32')

	for chr_name in list(sample_chrs.keys()):

		read_depths = sample_chrs[chr_name]['read_depth']
		gc_percent = ref_file[chr_name]['GC_Percent']
		
		total_read_depths = np.append(total_read_depths, read_depths[...])
		total_gc_percents = np.append(total_gc_percents, gc_percent[...])

		if len(read_depths) != len(gc_percent):
			print('Read depth length: ', len(read_depths))
			print('GC Percent length: ', len(gc_percent))
			break

		zipped = zip(gc_percent[...], read_depths[...])
		[gc_content[(int(np.around(100*i[0])))].append(i[1]) for i in zipped]

		[num_windows[(int(np.around(100*i)))].append(1) for i in gc_percent[...]]
		
	#If there are GC percentiles with 0 windows, we set these GC percentiles to 0.	
	for val in gc_content:
		if len(gc_content[val]) == 0:
			gc_content[val] = 0


	#POLYNOMIAL FIT:
	z = np.polyfit(total_gc_percents, total_read_depths, 3)
	f = np.poly1d(z)

	# calculate new x's and y's
	x_new = np.linspace(0, 100, 101)
	y_new = f(x_new)


	avg = [np.average(np.array(gc_content[perc])) for perc in gc_content]
	tot = [np.sum(np.array(num_windows[perc])) for perc in num_windows]

	"""
	The group 'GC' is created in the sample hdf5 file in order to store the following datasets:
	The dataset 'GC Percentiles' contains the integer GC percentiles from 0 to 100
	The dataset 'Average Depth' contains the average read depths across the entire genome of windows per GC percentile
	The dataset 'Num Windows' contains the number of windows across the entire genome per GC percentile
	The dataset 'Poly Fit' contains 101 values of a degree-3 polynomial fit to all the read depths vs GC percents across the entire genome
	"""

	#f_new = h5py.File('SRR726352_gc_faster_3.hdf5', 'w')
	grp = sample_file.create_group('GC')
	percentages = grp.create_dataset('GC Percentiles', data=np.arange(0,101))
	avg_gc = grp.create_dataset('Average Depth', data=np.array(avg))
	num_wind = grp.create_dataset('Num Windows', data=np.array(tot))
	polyfit = grp.create_dataset('Poly Fit', data=np.array(y_new))
	
	print('Average Read Depth: ', avg)
	print('Num Windows: ', tot)

	#f_new.close()
	sample_file.close()
	ref_file.close()








"""
Corrects the read depth of each window across the genome according to its GC Content
Assumption: the read depth should be the same across all GC percentiles
However, read depth varies in a predictable way according to GC percentage
Fix: we multiply the read depth by (mean read depth)/(average read depth in that GC percentile) to adjust for GC content
Inputs:
	-sample_hdf5: hdf5 file of the sample
	-ref_hdf5: hdf5 of the reference genome containing the GC percentage for each window across the genome ('panTro6_gc.hdf5')
Output:
	-Adds a dataset 'GC-Corrected Depths' to each chromosome group in the sample_hdf5 file
	-This dataset contains the GC-Corrected depth for each window
	-Calculation: GC-Corrected Depth = (original depth * mean depth across genome)/(average read depth of all windows at that GC percentile)
"""

def gc_correction(sample_hdf5, sample_gc_hdf5, ref_hdf5):

	sample_file = h5py.File(sample_hdf5, 'r+')
	sample_gc_file = h5py.File(sample_gc_hdf5, 'r+')
	ref_file = h5py.File(ref_hdf5, 'r')
	results_file = h5py.File('SRR726352_cnv_2.hdf5', 'w')

	#mean = sample_file['Stats']['Mean']

	gc_percentiles = sample_file['GC']['GC Percentiles']
	gc_avg_depth_raw = sample_file['GC']['Average Depth']


	#CHOOSE: SAVGOL FILTER OR POLYFIT:
	#POLYFIT:
	smoothed_gc = sample_file['GC']['Poly Fit'][...]
	#SAVGOL FILTER:
	smoothed_gc = savgol_filter(gc_avg_depth, 51, 3)

	chrs = sample_hdf5['depth']

	initial_mean = sample_hdf5['stats']['mean']

	new_means = np.zeros(len(list(chrs.keys())), dtype=np.float16)

	index = 0
	for chr in list(chrs.keys()):

		zipped = (ref_file[chr]['GC_Percent'][...], chrs[chr]['depth'][...])
		gc_corrected = [(i[1]*initial_mean)/float(smoothed_gc[int(np.around(100*i[0]))]) for i in zipped]

		#Saving means of each chromosome to calculate new corrected mean for CNV
		new_means[i] = np.mean(gc_corrected)
		i += 1

		gc_dset = chrs[chr].create_dataset('GC-Corrected Depths', data=np.array(gc_corrected))


	new_gc_mean = np.mean(new_means)
	cn_correction = 2.0/new_gc_mean

	for chr in list(chrs.keys()):

		cnv_correction = [cn_correction*val for val in chrs[chr]['GC-Corrected Depths'][...]]
		cn_dset = chrs[chr].create_dataset('Copy Number', data=np.array(cnv_correction))

	sample_file.close()
	ref_file.close()


"""

	gc_percentiles = sample_gc_file['GC']['GC Percentiles']
	gc_avg_depth_raw = sample_gc_file['GC']['Average Depth']
	gc_avg_depth = savgol_filter(gc_avg_depth_raw, 51, 3)	

	num_chrs = len(list(sample_file.keys()))
	means = np.array([], dtype=np.float16)
	for chr in list(sample_file.keys()):
		if len(sample_file[chr]['Depth'][...]) == 0:
			print('Chromosome with depth file = 0: ', chr)
			mean_chr = 0
		else:
			mean_chr = np.mean(sample_file[chr]["Depth"][...])
		means = np.append(means, mean_chr)
		#print(means)

	print('len of means: ', len(means))
	print('First 20 values: ', means[0:50])
	print('means: ', means)
	sm = np.sum(means)
	sm_test = np.sum(means[0:50])
	print('sum test: ', sm_test)
	total_mean = sm/num_chrs
	print('sum of means: ', sm)
	print('Total mean: ', str(total_mean))
	total_mean = np.mean(np.array(means))

	print('Total mean: ', str(total_mean))

	cn_correction = 2.0/total_mean

	for chr_name in list(sample_file.keys()):
		print(chr_name)
		grp = results_file.create_group(chr_name)

		read_depths = sample_file[chr_name]['Depth']
		gc_percents = ref_file[chr_name]['GC_Percent']

		gc_corr = np.array([], dtype='float32')

		copy_number = np.array([], dtype='float32')
		#gc_corr = sample_file.create_dataset('GC-Corrected Depth') #LATER

		if (len(read_depths) != len(gc_percents)):
			print('LENGTHS NOT THE SAME!')

		for i in np.arange(len(read_depths)):

			gc_val = np.around(100*gc_percents[i])
			gc_val = gc_val.astype(int)
			#print('GC VAL: ', gc_val)
			#gc_val = int((100*gc_percents[i])//1)
			depth_val = read_depths[i]

			corrected_depth = (depth_val*total_mean)/float(gc_avg_depth[gc_val])
			gc_corr = np.append(gc_corr, corrected_depth)

			cn_val = cn_correction*corrected_depth
			copy_number = np.append(copy_number, cn_val)


		start_dset = grp.create_dataset('GC-Corrected Depths', data=gc_corr)
		cn_dset = grp.create_dataset('Copy Number', data=copy_number)

	sample_file.close()
	ref_file.close()
	"""



"""
ref_gc_hdf5 calculates the GC percentile per sliding window across the entire genome
These windows should be the same as those calculated by create_h5
Output:
hdf5 file that contains a group for each chromosome
Each chromosome contains three datasets: Start, End, and GC_Percent (float value)
"""

def ref_gc_hdf5(fasta_file, hdf5_file, window_size, window_slide):

	f_hdf5 = h5py.File(hdf5_file, "w")
	f_fasta = pysam.Fasta(fasta_file)

	for chr_name in f_fasta.keys():
		
		grp = f_hdf5.create_group(chr_name)
		print('CHROM NAME: ', chr_name)
		chr_length = len(f_fasta[chr_name])

		num_windows = math.ceil(chr_length / window_slide - window_size / window_slide + 1)
		if chr_length < num_windows:
			num_windows = 1

		#Initialize numpy arrays to save time - around an order of magnitude faster to initalize first (or just use a list and convert at the end):
		start = np.zeros(num_windows, dtype='uint32')
		end = np.zeros(num_windows, dtype='uint32')
		count = np.zeros(num_windows, dtype='float32')

		index = 0
		for val in np.arange(num_windows):
			if val == num_windows - 1:
				end_index = chr_length
			else:
				end_index = index + window_size
			
			seq = f_fasta[chr_name][index:end_index]
			num_gc = seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')
			num_N = seq.count('N') + seq.count('n')


			# for base in seq:
			# 	if (base == 'G' or base == 'C' or base == 'g' or base == 'c'):
			# 		num_gc += 1
			# 	elif (base == 'N' or base == 'n'):
			# 		num_N += 1
			if window_size == num_N:
				gc_percent = 0.0
			else:
				gc_percent = float(num_gc)/(len(seq) - num_N)
			count[val] = gc_percent
			start[val] = index
			end[val] = end_index
			index += window_slide

		start_dset = grp.create_dataset('Start', data=start)
		end_dset = grp.create_dataset('End', data=end)
		count_dset = grp.create_dataset('GC_Percent', data=count)



	f_hdf5.close()









if __name__ == "__main__":

	# gen_stats = gen_stats2
	# get_stats = get_stats2

	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()

	parser_create = subparsers.add_parser("create")
	parser_create.add_argument("--fn_bam", '-f', required=True, help="Path to BAM file")
	parser_create.add_argument("--window", '-w', default=1000, help="Window size", type=int)
	parser_create.add_argument("--step", '-s', default=500, help="How much window should slide by", type=int)
	parser_create.add_argument("--ref", '-r', help="Which .fa reference genome file to use", type=str)
	parser_create.add_argument("--fn_out", '-o', default=None, help="Pathname to output file (include .h5 extension)")
	parser_create.add_argument("--chrs", '-c', nargs='+', default=[], help="Chromosome(s) to generate read depths for")
	parser_create.add_argument("--bam_index", '-i', help="Path to BAM index file (if not the same name and not in the same directory as the BAM file)", default=None)
	parser_create.add_argument("--stats", action='store_true', help="Generate stats (e.g. mean) for HDF5 file")
	# parser_create.add_argument("--copynum", action='store_true', help="Calculate copy number per window using average read depth and store them as arrays")
	parser_create.set_defaults(func=create_h5)

	parser_summary = subparsers.add_parser("summary")
	parser_summary.add_argument("--fn_h5", required=True, help="Pathname to HDF5 file")
	parser_summary.add_argument("--stats", action='store_true', help="Print out all stats") # todo: add list of stats available
	parser_summary.set_defaults(func=summary_h5)

	o = parser.parse_args()
	o.func(o)