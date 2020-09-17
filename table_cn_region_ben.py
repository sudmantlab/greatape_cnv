import numpy as np
import h5py
import argparse
import os
import json

"""
EXAMPLE COMMAND LINE ARGUMENTS:
python table_cn_region_ben.py --chimps --indices chr1:108308030-108344867 chr1:108409400-108428733 -o gstm_chimps.tsv
"""


def cn_table(args):
    humans = args.humans
    chimps = args.chimps
    indices = args.indices
    out_name = args.out_name
    window = 1000
    chrom = str(indices[0].split(':')[0])
    print('CHROM: ', chrom)

    og_indices = []
    out_file = open(out_name,'w')
    out_file.write('INDIVIDUAL\tREGION')
    for ind in indices:
        out_file.write('\t%s' % (ind))
    out_file.write('\n')

    rounded_indices = []
    for ind in indices:
        start_ind = (int(ind.split(':')[1].split('-')[0])//window)
        end_ind = (int(ind.split(':')[1].split('-')[1])//window)
        rounded_indices.append([start_ind,end_ind])
    print('OG INDICES: ', indices)
    print('ROUNDED INDICES: ', rounded_indices)

    if humans:
        h5_path = '/global/scratch2/almahalgren/humans/repmasked_smoothed_outputs/'
        config_file = '/global/scratch2/almahalgren/humans/configs/human_heatmap_config_3rem_nospace.json'
    elif chimps:
        h5_path = '/global/scratch2/almahalgren/chimps/repmasked_smoothed_outputs/'
        config_file = '/global/scratch2/almahalgren/chimps/configs/chimp_heatmap_config.json'

    with open(config_file,'r') as json_f:
        json_data = json.load(json_f)

    samp_dict = {}
    h5s = os.listdir(h5_path)
    for reg in json_data['REGIONS'].keys():
        for s in h5s:    
            if s[-3:] == '.h5' and s[:-3] in json_data['REGIONS'][reg]:
                print(s[:-3])
                samp_path = os.path.join(h5_path,s)
                h5_file = h5py.File(samp_path,'r')
                out_file.write('%s\t%s' % (s[:-3],reg))
                for ind in rounded_indices:
                    avg = np.around(np.mean(h5_file['depth'][chrom]['LOWESS_smoo_med_cn'][ind[0]:ind[1]+1][...]),decimals=2)
                    out_file.write('\t%s' % (avg))
                out_file.write('\n')
                h5_file.close()
    
    out_file.close()
    json_f.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--humans',action='store_true')
    parser.add_argument('--chimps',action='store_true')
    parser.add_argument('--indices',required=True,nargs='+',help='Start and end indices of all SVs (assuming you are copying and pasting from the genome browser?) ex. chr1:5000-6000 chr1:7000-8000')
    parser.add_argument('--out_name','-o',required=True,help='Name of out .tsv file')
    o = parser.parse_args()
    cn_table(o)






