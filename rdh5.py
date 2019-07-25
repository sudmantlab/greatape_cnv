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
        print("time to load contig array:", (time.time()-t) / 60, "minutes")

        s = 0
        e = window

        read_depth = np.array([], dtype=np.float16)
        start = np.array([], dtype=np.dtype('u4')) # unsigned int of 4 bytes (32 bits)
        end = np.array([], dtype=np.dtype('u4'))
        while e < chrom_length:
            # todo: create array for window % step != 0
            # todo: create array for window % step == 0

            # save array in hdf5
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


