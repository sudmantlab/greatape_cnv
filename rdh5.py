import pysam
import pysamstats
import tables
import numpy as np
import argparse
import time
import math
import h5py
import os

"""
python rdh5.py create --fn_bam SRR748186_short.bam --chrs chr22 --stats
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

	assert window > step, "Window must be greater than step"
	assert window > 0, "Window cannot equal zero"
	assert step > 0, "Step cannot equal zero"

	# if fn_out wasn't specified, name it the same as fn_bam
	if fn_out is None:
		fn_out = fn_bam[:len(fn_bam) - 4] + ".h5"
	else:
		assert fn_out[len(fn_out) - 3:] != ".h5", "Include \'.h5\' extension at the end of fn_out"

	# for generating stats only when HDF5 file exists
	if stats and os.path.exists(fn_out):
		gen_stats(fn_bam, window, step, fn_out, bam_index)

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

		n = window // step # number of steps in each window
		read_depth = np.zeros(size, dtype=np.float16)
		start = np.zeros(size, dtype=np.uint32)
		end = np.zeros(size, dtype=np.uint32)

		if window % step == 0:
			# an array that records averages for each block of size 'step' in 'chrom'
			step_arr = np.zeros(chrom_length // step, dtype=np.uint32)
			s = 0
			for i in range(0, chrom_length // step):
				step_arr[i] = np.sum(cvg[s:s + step])
				if i < size:  # 'step_arr' is longer than 'start' and 'end'
					start[i] = s
					end[i] = (s + window if s + window < chrom_length else chrom_length)
				s += step

			# tail case
			last = size
			if s < chrom_length:
				read_depth[size - 1] = np.mean(cvg[s:])
				last -= 1
			# setting up for loop
			arr_sum = np.sum(step_arr[0:n])
			read_depth[0] = arr_sum / window
			s, e = 0, n
			for i in range(1, last):
				arr_sum = arr_sum - step_arr[s] + step_arr[e]
				read_depth[i] = arr_sum / window
				s += 1
				e += 1

		else:
		# straightforward way of calculating average read depth for each window
			s, e = 0, window
			for i in range(size):
				start[i] = s
				end[i] = e
				read_depth[i] = np.mean(cvg[s:e])
				s += step
				e += step
				# tail case
				if e > chrom_length:
					e = chrom_length



		# save arrays in hdf5
		out_file.create_carray(chrom_group, name="read_depth", atom=tables.Float16Atom(), obj=read_depth)
		out_file.create_carray(chrom_group, name="start", atom=tables.UInt32Atom(), obj=start)
		out_file.create_carray(chrom_group, name="end", atom=tables.UInt32Atom(), obj=end)

		print("--------------------------------------------------")

	bam_file.close()
	out_file.close()

	# generate stats
	if stats:
		gen_stats(fn_bam, window, step, fn_out, bam_index)

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
	chrs = list(depth_group.__iter__())
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

def summary_h5(args):
	fn_h5 = args.fn_h5
	stats = args.stats
	if stats:
		print(get_stats(fn_h5))


if __name__ == "__main__":

	# gen_stats = gen_stats2
	# get_stats = get_stats2

	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()

	parser_create = subparsers.add_parser("create")
	parser_create.add_argument("--fn_bam", '-f', required=True, help="Path to BAM file")
	parser_create.add_argument("--window", '-w', default=1000, help="Window size", type=int)
	parser_create.add_argument("--step", '-s', default=500, help="How much window should slide by", type=int)
	parser_create.add_argument("--fn_out", '-o', default=None, help="Pathname to output file (include .h5 extension)")
	parser_create.add_argument("--chrs", '-c', nargs='+', default=[], help="Chromosome(s) to generate read depths for")
	parser_create.add_argument("--bam_index", '-i', help="Path to BAM index file (if not the same name and not in the same directory as the BAM file)", default=None)
	parser_create.add_argument("--stats", action='store_true', help="Generate stats (e.g. mean) for HDF5 file")
	parser_create.set_defaults(func=create_h5)

	parser_summary = subparsers.add_parser("summary")
	parser_summary.add_argument("--fn_h5", required=True, help="Pathname to HDF5 file")
	parser_summary.add_argument("--stats", action='store_true', help="Print out all stats") # todo: add list of stats available
	parser_summary.set_defaults(func=summary_h5)

	o = parser.parse_args()
	o.func(o)
