import pysam
import pysamstats
import tables
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time


def build_h5(bam):
    bam_file = pysam.AlignmentFile(bam, "rb") # will open index file automatically if in the same folder
    contig_lengths = {x[0]:x[1] for x in zip(bam_file.references, bam_file.lengths)}
    refs = bam_file.references
    # for r in contig_lengths:
    t = time.time()
    carray = pysamstats.load_coverage(bam_file, chrom=refs[0], start=0, end=contig_lengths[refs[0]])
    print("time to load contig: {}".format(time.time()-t))
    plt.plot(carray.pos, carray.reads_all)
    plt.savefig('plot.png', bbox_inches='tight')
    bam_file.close()


# samfile = pysam.AlignmentFile("ex1.bam", "rb" )
# for pileupcolumn in samfile.pileup("chr1", 100, 120): # (contig, start, stop)
#     print ("\ncoverage at base %s = %s" %
#            (pileupcolumn.pos, pileupcolumn.n))
#     for pileupread in pileupcolumn.pileups:
#         if not pileupread.is_del and not pileupread.is_refskip:
#             # query position is None if is_del or is_refskip is set.
#             print ('\tbase in read %s = %s' %
#                   (pileupread.alignment.query_name,
#                    pileupread.alignment.query_sequence[pileupread.query_position]))
#
#
# mybam = pysam.AlignmentFile('/path/to/your/bamfile.bam')
#
#
# samfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="Path to BAM file")
    o = parser.parse_args()
    build_h5(o.bam)


