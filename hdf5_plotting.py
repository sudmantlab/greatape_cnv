# hdf5 example
# Storing the CNV data of an individual in the following format:
# Each chromosome is its own group, e.g. '/chr1/'
# Within each chromosome group, the start and end positions are stored in a two-dimensional array: '/chr1/Window' -> [[0, 10], [10, 20], ...]
# Also within each chromosome group, the number of reads in each region are stored in a one-dimensional array: '/chr1/Depth' -> [50, 40, ...]

import h5py
import numpy as np
import pysam
import pysamstats
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def depth_histogram(sample, chrom, start, end, window_size):
	hdf5_file = sample + '.hdf5'
	f = h5py.File(hdf5_file, "r")
	depths = f[chrom]['Depth'][...]
	ax = sns.distplot(depths, bins=np.arange(0,12.2,0.2), kde=False)
	plt.xlabel("Average Read Depth")
	plt.ylabel("Number of Windows")
	plt.title("Average Read Depth per " + str(window_size) + " bp for " + chrom)
	ax.set(xlim=(0,12))
	fig = ax.get_figure()
	pic_name = sample + '_' + chrom + '_hist.png'
	fig.savefig(pic_name)

def depth_kde(sample, chrom, start, end, window_size):
	hdf5_file = sample + '.hdf5'
	f = h5py.File(hdf5_file, "r")
	depths = f[chrom]['Depth'][...]
	ax = sns.distplot(depths, bins=np.arange(0,12.2,0.2), hist=False)
	plt.xlabel("Average Read Depth")
	plt.ylabel("Fraction of Windows")
	plt.title("KDE: Average Read Depth per " + str(window_size) + " bp for " + chrom)
	ax.set(xlim=(0,12))
	fig = ax.get_figure()
	pic_name = sample + '_' + chrom + '_kde.png'
	fig.savefig(pic_name)


def depth_scatterplot(sample, chrom, start, end, window_size):
	hdf5_file = sample + '.hdf5'
	f = h5py.File(hdf5_file, "r")
	depths = f[chrom]['Depth'][...]
	bp = f[chrom]['Start'][...]
	plt.figure(figsize=(15,8))
	ax = plt.scatter(x=bp, y=depths, s=30)
	#ax = sns.scatterplot(x=bp, y=depths, s=10)
	plt.xlabel("Base Pair Coordinates", fontsize=14)
	plt.ylabel("Average Read Depth per " + str(window_size) + "bp", fontsize=14)
	plt.title("Average Read Depth across " + chrom, fontsize=20)
	fig = ax.get_figure()
	pic_name = sample + '_' + chrom + '_scat.png'
	fig.savefig(pic_name)



my_parser = argparse.ArgumentParser()
my_parser.add_argument('--sample', action='store', type=str, required=True, help='sample name e.g SRR726352')
my_parser.add_argument('--plot_type', action='store', type=str, required=True, help='hist, scat, or kde')
my_parser.add_argument('--window_size', action='store', type=int, required=True, help='standard: 1000')
my_parser.add_argument('--chrom', action='store', type=str, required=True, help='e.g. chr1')
my_parser.add_argument('--start', action='store', type=int, required=False)
my_parser.add_argument('--end', action='store', type=int, required=False)

args = vars(my_parser.parse_args())

print('SAMPLE: ', args['sample'])
print('PLOT TYPE: ', args['plot_type'])
print('WINDOW SIZE: ', args['window_size'])
print('CHROM: ', args['chrom'])
print('START: ', args['start'])
print('END: ', args['end'])

if args['plot_type'] == 'hist':
	depth_histogram(args['sample'], args['chrom'], args['start'], args['end'], args['window_size'])
elif args['plot_type'] == 'kde':
	depth_kde(args['sample'], args['chrom'], args['start'], args['end'], args['window_size'])
elif args['plot_type'] == 'scat':
	depth_scatterplot(args['sample'], args['chrom'], args['start'], args['end'], args['window_size'])


#create('SRR726352.bam', 'SRR726352_chr1.hdf5', 1000, 200)
#depth_histogram('SRR726352_chr1.hdf5', 'chr1', 1000)
#depth_scatterplot('SRR726352_chr1.hdf5', 'chr1')
