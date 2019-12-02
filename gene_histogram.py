import os
# import pysam
# import pybedtools
import argparse
import tables
import numpy as np
import math
import matplotlib.pyplot as plt

window = None
step = None

def parse_bed(bed):
    parsed = []
    bedfile = open(bed)

    line = bedfile.readline()
    while line != "":
        parsed.append(line.split())
        line = bedfile.readline()

    bedfile.close()
    return parsed

def parse_h5_files(text):
    global window
    global step
    textfile = open(text)
    path = textfile.readline().split('\n')[0]
    files = []

    while path != "":
        # todo
        file = tables.open_file(path)
        w = file.get_node(where=file.root.stats, name='window')[...][()]
        s = file.get_node(where=file.root.stats, name='step')[...][()]
        if window is None or step is None:
            window = w
            step = s
        if w != window or s != step:
            raise Exception("Window and step sizes must be same among all HDF5 files")
        files.append(file)
        path = textfile.readline().split('\n')[0]

    return files

def get_chr_subarray_exclusive(file, chr, start, end):
    """START and END are exclusive, i.e. only copynum for
    windows that are within START and END are returned"""
    start_window = math.ceil(np.divide(start, step)) # todo: is this right?
    end_window = math.ceil(np.divide(end - window, step)) # todo: is this right?
    chrom_group = file.get_node(where=file.root.depth, name=chr)
    start_idx = file.get_node(where=file.root.depth, name=chr).start[start_window]
    end_idx = file.get_node(where=file.root.depth, name=chr).end[end_window - 1]
    return start_idx, end_idx, chrom_group.copy_number[start_window:end_window] # end index of splice is exclusive

def get_chr_subarray_inclusive(file, chr, start, end):
    """START and END are inclusive, i.e. the smallest range
    of windows that include START and END are returned"""
    start_window = start // step # todo: is this right?
    end_window = math.ceil(np.divide(end - window, step)) + 1 # todo: is this right?
    chrom_group = file.get_node(where=file.root.depth, name=chr)
    start_idx = file.get_node(where=file.root.depth, name=chr).start[start_window]
    end_idx = file.get_node(where=file.root.depth, name=chr).end[end_window - 1]
    return start_idx, end_idx, chrom_group.copy_number[start_window:end_window] # end index of splice is exclusive

get_chr_subarray = get_chr_subarray_inclusive

def make_histogram(name, chr, start, end, h5_files, bins, fn_out):
    copy_nums = []
    start_idx, end_idx = None, None
    for f in h5_files:
        if start_idx is None or end_idx is None:
            start_idx, end_idx, cn_arr = get_chr_subarray(f, chr, start, end)
        else:
            s, e, cn_arr = get_chr_subarray(f, chr, start, end)
            if s != start_idx or e != end_idx:
                print("Start and end chromosome indices are not equal")

        copy_nums.append(np.average(cn_arr))

    hist, bin_edges = np.histogram(copy_nums, bins=bins)
    plt.hist(copy_nums, bins)
    if name is None:
        plt.title("{chr}:{start}-{end}".format(chr=chr, start=start_idx, end=end_idx))
    else:
        plt.title("{name} {chr}:{start}-{end}".format(name=name, chr=chr, start=start_idx, end=end_idx))
    plt.xlabel("Copy number")
    plt.ylabel("Number of individuals")
    plt.savefig(fn_out)
    # todo

def close_files(files):
    for f in files:
        f.close()

def subarray_func(exclusive):
    global get_chr_subarray
    if exclusive:
        get_chr_subarray = get_chr_subarray_exclusive

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--chr", required=True, help="Chromosome name")
    parser.add_argument("--start", required=True, help="Chromosome start index", type=np.int64)
    parser.add_argument("--end", required=True, help="Chromsosome end index", type=np.int64)
    parser.add_argument("--name", required=False, default=None, help="Gene name")
    parser.add_argument("--exclusive", action='store_true', help="Generate stats (e.g. mean) for HDF5 file")
    parser.add_argument("--h5_list", '-l', required=True, help="Path to TXT file which is a list of all HDF5 file paths")
    parser.add_argument("--bins", default=10, required=False, help="Number of histogram bins (optional)", type=int)
    parser.add_argument("--fn_out", default="plot.png", required=False, help="Filename of output histogram plot")

    args = parser.parse_args()

    subarray_func(args.exclusive)
    files = parse_h5_files(args.h5_list)
    make_histogram(args.name, args.chr, args.start, args.end, files, args.bins, args.fn_out)
    close_files(files)
