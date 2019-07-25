import pysam
import pysamstats
import tables
import numpy as np
import argparse
import time


def build_h5(bam, window, step, out):
    bam_file = pysam.AlignmentFile(bam, "rb") # will open index file automatically if in the same folder
    refs = bam_file.references
    contig_lengths = {x[0]:x[1] for x in zip(refs, bam_file.lengths)}

    out_file = tables.open_file(out, mode='w')

    for r in contig_lengths.items():
        chrom_name = r[0]
        chrom_length = r[1]
        out_file.create_group(out_file.root, chrom_name)

        t = time.time()
        cvg_array = pysamstats.load_coverage(bam_file, chrom=chrom_name, start=0, end=chrom_length)
        print("time to load contig array for", chrom_name, ":", (time.time()-t) / 60, "minutes")

        size = chrom_length // step - window // step + 1
        n = window // step
        read_depth = np.zeros(size, dtype=np.float16)
        start = np.zeros(size, dtype=np.dtype('u4')) # unsigned int of 4 bytes (32 bits)
        end = np.zeros(size, dtype=np.dtype('u4'))

        # todo: write documentation for all of this
        if window % step == 0:
            step_arr = np.zeros(chrom_length // step + 1, dtype=np.float16)
            s, e = 0, step
            for i in range(0, chrom_length // step):
                step_arr[i] = np.mean(cvg_array[s:e]) / n # dividing by n for faster avg calculation in next for loop
                start[i] = s
                end[i] = s + window
                s, e = e, e + step
            # tail case
            step_arr[chrom_length // step] = (np.mean(cvg_array[e:]))
            start[chrom_length // step] = s
            end[chrom_length // step] = chrom_length

            mean = np.mean(step_arr[0:n])
            read_depth[0] = mean
            for i in range(1, read_depth.size):
                mean = mean - step_arr[i - 1] + step_arr[n + 1]
                n += 1
                read_depth[i] = mean

        else:
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


        # todo: create array for window % step != 0

        # save array in hdf5; todo: carray or earray?
        out_file.create_earray(out_file.root.chrom_name, name="read_depth", atom=tables.Float16Atom(), obj=read_depth)
        out_file.create_earray(out_file.root.chrom_name, name="start", atom=tables.Float16Atom(), obj=start)
        out_file.create_earray(out_file.root.chrom_name, name="end", atom=tables.Float16Atom(), obj=end)


    bam_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_bam", required=True, help="Path to BAM file")
    parser.add_argument("--fn_out", required=True, help="Pathname to output file")
    parser.add_argument("--window", help="Window size", default=1000)
    parser.add_argument("--step", help="Step size", default=500)
    o = parser.parse_args()
    build_h5(o.fn_bam, o.window, o.step, o.fn_out)


