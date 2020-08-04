# greatape_cnv
Calculating the copy number variation of the great apes.

To view the track hubs listed below, copy-paste the hub.txt link into the link below and click 'Add Hub'. They may already show up as a track hub when you click the link.
http://www.genome.ucsc.edu/cgi-bin/hgHubConnect?hubId=20163&hgHub_do_disconnect=on&hgHubConnect.remakeTrackHub=on&hgsid=866861907_Za3YyTP6nvfDH9soIbl6H2RT2wRs

GREAT APE GENOMES:
We are characterizing structural variation among 71 chimpanzees and bonobos and comparing structural diversity to human populations.

Track hub of copy number:
https://github.com/sudmantlab/greatape_cnv/raw/master/GreatApe_Smoothed/hub.txt

Track hub of HMM-smoothed copy number:
https://github.com/sudmantlab/greatape_cnv/raw/master/chimp_cn_hmm/hub.txt

hdf5 files with copy number values:
/global/scratch2/almahalgren/chimps/repmasked_smoothed_outputs/

hdf5 files with HMM-smoothed copy number values:
/global/scratch2/almahalgren/chimps/smoothed_outputs/

Path to copy number within hdf5 files:
file['depth'][chrom]['LOWESS_smoo_med_cn']

Path to start and end arrays within hdf5 files:
file['depth'][chrom]['start'], file['depth'][chrom]['end']


ANCIENT DNA:
We are examining the structural variation in 242 individuals of European and Central Asian ancestry from 2-10,000 years ago. 

Track hub of copy number:
https://github.com/sudmantlab/greatape_cnv/raw/master/Ancient_RepMasked_Decoy/hub.txt

hdf5 files with copy number values on the Savio cluster:
/global/scratch2/almahalgren/ancient_dna/hg38_data/repmasked_smoothed_outputs/

Path to copy number within hdf5 files:
file['depth'][chrom]['LOWESS_smoo_med_cn']

Path to start and end arrays within hdf5 files:
file['depth'][chrom]['start'], file['depth'][chrom]['end']





