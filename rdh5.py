import pysam
import pysamstats
import tables
import numpy as np
import argparse
import time


# todo: write better documentation for all of this
def build_h5(bam, window, step, out):
    bam_file = pysam.AlignmentFile(bam, "rb") # will open index file automatically if in the same folder
    refs = bam_file.references
    contig_lengths = {x[0]:x[1] for x in zip(refs, bam_file.lengths)}

    out_file = tables.open_file(out, mode='w')

    for r in contig_lengths.items():
        chrom_name = r[0]
        chrom_length = r[1]
        print("Chromosome name:", chrom_name)
        print("chromosome length:", chrom_length)
        chrom_group = out_file.create_group(out_file.root, chrom_name)
        if chrom_length == 0:
            continue

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


        size = chrom_length // step - window // step + 1 + 1 # number of windows in chrom plus one extra
        print("Number of windows in chrom {}: {}".format(chrom_name, size))
        n = window // step # number of steps in each window
        read_depth = np.zeros(size, dtype=np.float16)
        start = np.zeros(size, dtype=np.int32) # unsigned int of 4 bytes (32 bits)
        end = np.zeros(size, dtype=np.int32)

        # todo: test whether this gives a significant speedup
        # https://www.geeksforgeeks.org/window-sliding-technique/
        if window % step == 0:
            # an array that records averages for each block of size 'step' in 'chrom'
            step_arr = np.zeros(chrom_length // step + 1, dtype=np.float16)
            s, e = 0, step
            for i in range(0, chrom_length // step):
                step_arr[i] = np.mean(cvg[s:e]) / n # dividing by 'n' for faster mean calculation in next for loop
                if i < size:  # step_arr is longer than start and end
                    start[i] = s
                    end[i] = (s + window if s + window < chrom_length else chrom_length)
                s, e = e, e + step
            # tail case
            step_arr[chrom_length // step] = (np.mean(cvg[e:]))
            # start[size - 1] = s
            # end[size - 1] = chrom_length

            # initial mean calculation (tail case but it's at the beginning)
            mean = np.mean(step_arr[0:n])
            read_depth[0] = mean
            for i in range(1, read_depth.size):
                mean = mean - step_arr[i - 1] + step_arr[n]
                n += 1
                read_depth[i] = mean

        else:
            # straightforward way of calculating average read depth for each window
            s, e = 0, window
            for i in range(size):
                start[i] = s
                end[i] = e
                read_depth[i] = np.mean(cvg_array[s:e])
                s += step
                e += step
                # tail case
                if e > chrom_length:
                    e = chrom_length



        # save arrays in hdf5
        out_file.create_carray(chrom_group, name="read_depth", atom=tables.Float16Atom(), obj=read_depth)
        out_file.create_carray(chrom_group, name="start", atom=tables.Int32Atom(), obj=start)
        out_file.create_carray(chrom_group, name="end", atom=tables.Int32Atom(), obj=end)


    bam_file.close()
    out_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_bam", required=True, help="Path to BAM file")
    parser.add_argument("--fn_out", required=True, help="Pathname to output file")
    parser.add_argument("--window", default=1000, help="Window size")
    parser.add_argument("--step", default=500, help="How much window should slide by")
    o = parser.parse_args()
    build_h5(o.fn_bam, o.window, o.step, o.fn_out)


