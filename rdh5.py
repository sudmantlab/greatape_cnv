import pysam
import pysamstats
import tables
import numpy as np
import argparse
import time
import math


# todo: write better documentation for all of this
def build_h5(bam, window, step, out, chrs):
    bam_file = pysam.AlignmentFile(bam, "rb") # will open index file automatically if in the same folder
    refs = bam_file.references
    contig_lengths = {x[0]:x[1] for x in zip(refs, bam_file.lengths)}

    if len(chrs) != 0:
        temp = {}
        for c in chrs:
            temp[c] = contig_lengths[c]
        contig_lengths = temp

    out_file = tables.open_file(out, mode='w')

    for r in contig_lengths.items():
        print("--------------------------------------------------")
        chrom_name = r[0]
        chrom_length = r[1]
        print("Chromosome name:", chrom_name)
        print("chromosome length:", chrom_length)
        chrom_group = out_file.create_group(out_file.root, chrom_name)

        # Generate array of read-depth per base using pysamstats
        # cvg_array.dytpe = dtype((numpy.record, [('chrom', 'S27'), ('pos', '<i4'), ('reads_all', '<i4'), ('reads_pp', '<i4')]))
        # A numpy.recarray with each column with following format: (name, dtype)
        # 'S27' = 27-byte string, '<i4' = little-endian 32-bit int
        # chrom = reference/chromosome name
        # pos = base position in reference/chromosome
        # reads_all = Number of reads aligned at the position. N.b., this is really the total, i.e., includes reads where the mate is unmapped or otherwise not properly paired.
        # reads_pp = Number of reads flagged as properly paired by the aligner
        # pysamstats.load_coverage does not include positions where read depth is zero
        t = time.time()
        cvg_array = pysamstats.load_coverage(bam_file, chrom=chrom_name, start=0, end=chrom_length)
        cvg_pos = cvg_array.pos
        cvg_reads = cvg_array.reads_all
        cvg = np.zeros(chrom_length, dtype=np.int32)
        for i in range(len(cvg_pos)):
            cvg[cvg_pos[i]] = cvg_reads[i]
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
            step_arr = np.zeros(chrom_length // step, dtype=np.int32)
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_bam", required=True, help="Path to BAM file")
    parser.add_argument("--fn_out", required=True, help="Pathname to output file")
    parser.add_argument("--window", default=1000, help="Window size")
    parser.add_argument("--step", default=500, help="How much window should slide by")
    parser.add_argument("--chr", nargs='+', default=[], help="Chromosome(s) to generate read depths for")
    o = parser.parse_args()
    if int(o.window) < int(o.step):
        print("Window must be greater than step")
        exit(1)
    if int(o.window) == 0:
        print("Window cannot equal zero")
        exit(1)
    if int(o.step) == 0:
        print("Step cannot equal zero")
        exit(1)
    if o.fn_out[len(o.fn_out) - 3:] != ".h5":
        print("Include \'.h5\' extension at the end of output file name")
        exit(1)
    build_h5(o.fn_bam, int(o.window), int(o.step), o.fn_out, o.chr)

