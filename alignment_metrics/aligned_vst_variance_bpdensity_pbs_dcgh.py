import h5py
import sys
import time
import numpy as np
import argparse
import os
import json

"""
Returns PBS, Vst, variance, BP Density, or dCGH across the genome
"""

def vst_formula(v1,n1,v2,n2,v_tot):
    if v_tot == 0:
        return 0
    else:
        return float(1 - (((v1*n1) + (v2*n2))/(v_tot*(n1+n2))))

def pbs_formula(v1,n1,v2,n2,v3,n3,v12,v13,v23):
    vst12 = vst_formula(v1, n1, v2, n2, v12)
    vst13 = vst_formula(v1, n1, v3, n3, v13)
    vst23 = vst_formula(v2, n2, v3, n3, v23)

    pbs1 = 0.5*(vst12 + vst13 - vst23)
    pbs2 = 0.5*(vst12 + vst23 - vst13)

    return max([pbs1,pbs2])

def calc_vst(start_ind,end_ind,chrom,first_region,second_region,first_dict_cn, second_dict_cn,num_first,num_second):
    first_cns = [first_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(first_dict_cn.keys())]
    first_cns = np.mean(np.array(first_cns),axis=1)
    first_var = np.var(first_cns)
    if first_region != '' and second_region != '':
        second_cns = [second_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(second_dict_cn.keys())]
        second_cns = np.mean(np.array(second_cns),axis=1)
        second_var = np.var(second_cns)
        combo_cns = np.concatenate((first_cns,second_cns))
        combo_var = np.var(combo_cns)
        if combo_var == 0:
            out_vst = 0
        else:
            out_vst = float(1 - (((first_var*num_first) + (second_var*num_second))/(combo_var*(num_first+num_second))))
    return out_vst

def calc_variance(start_ind,end_ind,chrom,first_region,second_region,first_dict_cn, second_dict_cn):
    first_cns = [first_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(first_dict_cn.keys())]
    first_cns = np.mean(np.array(first_cns),axis=1)
    first_var = np.var(first_cns)
    if first_region != '' and second_region != '':
        second_cns = [second_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(second_dict_cn.keys())]
        second_cns = np.mean(np.array(second_cns),axis=1)
        second_var = np.var(second_cns)
        combo_cns = np.concatenate((first_cns,second_cns))
        combo_var = np.var(combo_cns)
        return combo_var
    else:
        return first_var

def calc_bp_density(start_ind,end_ind,chrom,first_region,second_region,first_dict_cn, second_dict_cn,example):
    if first_region != '' and second_region != '':
        return np.mean([np.mean([np.sum(first_dict_cn[samp][chrom][start_ind:end_ind]) for samp in list(first_dict_cn.keys())]),np.mean([np.sum(second_dict_cn[samp][chrom][start_ind:end_ind]) for samp in list(second_dict_cn.keys())])])
    else:
        np.mean([np.sum(first_dict_cn[samp][chrom][start_ind:end_ind]) for samp in list(first_dict_cn.keys())])
    """
    if end_ind > len(first_dict_cn[example][chrom]):
        end_ind = len(first_dict_cn[example][chrom])
    mask_arr = [first_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(first_dict_cn.keys())]
    if start_ind == 0:
        mask_arr = [([0]+first_dict_cn[samp][chrom][start_ind:end_ind-1]) for samp in list(first_dict_cn.keys())]
    else:
        mask_arr_shifted = [first_dict_cn[samp][chrom][start_ind-1:end_ind-1] for samp in list(first_dict_cn.keys())]
    bp_density = np.mean([len(np.where((mask_arr[i]-mask_arr_shifted[i]) != 0)[0]) for i in np.arange(len(mask_arr))])
    return bp_density

    mask_arr = [np.array(first_dict_cn[samp][chrom][start_ind:end_ind]) for samp in list(first_dict_cn.keys())]
    mask_arr_shifted = [np.append(np.array([0]),np.array(first_dict_cn[samp][chrom][start_ind:end_ind-1])) if start_ind == 0 else np.array(first_dict_cn[samp][chrom][start_ind-1:end_ind-1]) for samp in list(first_dict_cn.keys())]
    bp_density = np.mean([len(np.where((mask_arr[i]-mask_arr_shifted[i]) != 0)[0]) for i in np.arange(len(mask_arr))])
    return bp_density
    """

def calc_pbs(start_ind,end_ind,chrom,first_region,second_region,third_region,first_dict_cn, second_dict_cn, third_dict_cn,num_first,num_second,num_third):
    first_cns = [first_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(first_dict_cn.keys())]
    first_cns = np.mean(np.array(first_cns),axis=1)
    first_var = np.var(first_cns)

    second_cns = [second_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(second_dict_cn.keys())]
    second_cns = np.mean(np.array(second_cns),axis=1)
    second_var = np.var(second_cns)

    third_cns = [third_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(third_dict_cn.keys())]
    third_cns = np.mean(np.array(third_cns),axis=1)
    third_var = np.var(third_cns)

    first_second_cns = np.concatenate((first_cns,second_cns))
    first_second_var = np.var(first_second_cns)
    first_third_cns = np.concatenate((first_cns,third_cns))
    first_third_var = np.var(first_third_cns)
    second_third_cns = np.concatenate((second_cns,third_cns))
    second_third_var = np.var(second_third_cns)

    return pbs_formula(first_var, num_first, second_var, num_second, third_var, num_third, first_second_var, first_third_var, second_third_var)

def calc_dcgh(start_ind,end_ind,chrom,first_region,dcgh_ind,zero_replacement,first_dict_cn):
    cns = [first_dict_cn[samp][chrom][start_ind:end_ind] for samp in list(first_dict_cn.keys()) if samp != dcgh_ind]
    cns = np.mean(np.array(cns), axis=1)
    dcgh_ind_cns = first_dict_cn[dcgh_ind][chrom][start_ind:end_ind]
    dcgh_ind_cns = np.mean(dcgh_ind_cns)
    if dcgh_ind_cns == 0:
        dcgh_ind_cns = zero_replacement
    cns[np.where(cns == 0)[0]] = zero_replacement
    return np.mean([(np.log2(ind) - np.log2(dcgh_ind_cns)) for ind in cns])



def vst_variance_bpdensity_pbs_dcgh(args):
    h5_mask_file = args.h5_mask_file
    json_file = args.json_file
    first_region = args.first_region
    second_region = args.second_region
    #ONLY THIRD REGION USED FOR PBS - OUTGROUP
    third_region = args.third_region
    h5_path = args.h5_path
    out_name = args.out_name
    chimp = args.chimp
    human = args.human
    hmm = args.hmm
    variance = args.variance
    #IF BP DENSITY: ONLY PERFORMS BP DENSITY ON FIRST_REGION
    bp_density = args.bp_density
    vst = args.vst
    pbs = args.pbs
    dcgh = args.dcgh
    dcgh_json = args.dcgh_json
    zero_replacement = float(args.zero_replacement)
    max_max = args.max_max
    max_window = int(args.max_window)
    chimp_chroms = args.chimp_chroms
    nobon = args.nobon
    bon = args.bon


    print('VST: ', vst)
    print("VARIANCE: ", variance)
    print("BP DENSITY: ", bp_density)
    print('dCGH: ', dcgh)
    print('PBS: ', pbs)
    print('MAX MAX: ', max_max)
    print('CHIMP: ', chimp)
    print('HUMAN: ', human)

    if vst:
        out_grp = 'vst'
    elif variance:
        out_grp = 'variance'
    elif bp_density:
        out_grp = 'bp_density'
    elif pbs:
        out_grp = 'pbs'
    elif dcgh:
        out_grp = 'dcgh'
    print(out_grp)

    if dcgh:
        with open(dcgh_json,'r') as dcgh_json_f:
            dcgh_json_file = json.load(dcgh_json_f)
        dcgh_ind = dcgh_json_file[first_region]
    """
    If cnhmm, then the end of the name should be '_cnhmm.h5'
    """

    if hmm:
        suffix_len = 9
    else:
        suffix_len = 3

    if bp_density:
        overall_grp = 'sv_masks'
        #cn_grp = 'mask'
        cn_grp = 'all_breakpoints'
        suffix_len = 17
    else:
        overall_grp = 'depth'
        if hmm:
            cn_grp = 'cn_hmm'
        else:
            cn_grp = 'LOWESS_smoo_med_cn'


    t = time.time()
    
    chimp_chrom_file = open(chimp_chroms,'r')
    chimp_chrom_list = [chimp_chrom_file.readline()[:-1] for i in np.arange(4346)]
    chimp_chrom_list.append('')
    print('CHIMP CHROM LIST: ', chimp_chrom_list)
    chimp_chrom_file.close()
    
    h5_mask = h5py.File(h5_mask_file,'r')

    with open(json_file,'r') as json_f:
        json_data = json.load(json_f)

    first_samp_list = []
    if first_region == '':
        for reg in list(json_data['REGIONS'].keys()):
            if (chimp and nobon and reg != 'PANISCUS') or (not nobon):
                for samp in json_data['REGIONS'][reg]:
                    first_samp_list.append(samp)
    else:
        first_samp_list = json_data['REGIONS'][first_region]
    print('LEN FIRST SAMP LIST: ', len(first_samp_list))

    second_samp_list = []
    if second_region == '':
        second_samp_list = []
    else:
        second_samp_list = json_data['REGIONS'][second_region]
    print('LEN SECOND SAMP LIST: ', len(second_samp_list))

    third_samp_list = []
    if third_region == '':
        third_samp_list = []
    else:
        third_samp_list = json_data['REGIONS'][third_region]
    print('LEN THIRD SAMP LIST: ', len(third_samp_list))

    first_dict_cn = {}
    second_dict_cn = {}
    third_dict_cn = {}
    first = True
    h5s = os.listdir(h5_path)
    t = time.time()
    for s in h5s:
        if (s[-3:] == '.h5' and s[:-suffix_len] in first_samp_list):
            if first:
                example = s[:-suffix_len]
                not first
            print(s[:-3])
            samp_path = os.path.join(h5_path,s)
            samp_h5 = h5py.File(samp_path,'r')
            first_dict_cn[s[:-suffix_len]] = {}
            for chrom in list(samp_h5[overall_grp].keys()):
                first_dict_cn[s[:-suffix_len]][chrom] = samp_h5[overall_grp][chrom][cn_grp][...]
            samp_h5.close()
        elif (s[-3:] == '.h5' and s[:-suffix_len] in second_samp_list):
            print(s[:-3])
            samp_path = os.path.join(h5_path,s)
            samp_h5 = h5py.File(samp_path,'r')
            second_dict_cn[s[:-suffix_len]] = {}
            for chrom in list(samp_h5[overall_grp].keys()):
                second_dict_cn[s[:-suffix_len]][chrom] = samp_h5[overall_grp][chrom][cn_grp][...]
            samp_h5.close()
        elif (s[-3:] == '.h5' and s[:-suffix_len] in third_samp_list):
            print(s[:-3])
            samp_path = os.path.join(h5_path,s)
            samp_h5 = h5py.File(samp_path,'r')
            third_dict_cn[s[:-suffix_len]] = {}
            for chrom in list(samp_h5[overall_grp].keys()):
                third_dict_cn[s[:-suffix_len]][chrom] = samp_h5[overall_grp][chrom][cn_grp][...]
            samp_h5.close()
    print('TIME TO LOAD INDS: ', (time.time()-t)/60.0)

    h5_mask_dict = {}
    for chrom in h5_mask.keys():
        h5_mask_dict[chrom] = {}
        h5_mask_dict[chrom]['mask'] = h5_mask[chrom]['mask'][...]
        h5_mask_dict[chrom]['human_start'] = h5_mask[chrom]['human_start'][...]
        h5_mask_dict[chrom]['human_end'] = h5_mask[chrom]['human_end'][...]
        h5_mask_dict[chrom]['chimp_start'] = h5_mask[chrom]['chimp_start'][...]
        h5_mask_dict[chrom]['chimp_end'] = h5_mask[chrom]['chimp_end'][...]
        h5_mask_dict[chrom]['chimp_chrom'] = h5_mask[chrom]['chimp_chrom'][...]
        
    num_first = len(first_samp_list)
    num_second = len(second_samp_list)
    num_third = len(third_samp_list)

    vst_dict = {}

    overall_t = time.time()

    for human_chrom in list(h5_mask.keys()):
        t = time.time()
        vst_dict[human_chrom] = np.ones(len(h5_mask_dict[human_chrom]['human_start']))*-100
        matching_inds = np.where(h5_mask_dict[human_chrom]['mask'] == 1)[0]
        print('CHROM: %s LEN CHROM: %s NUM MATCHING INDS: %s LEN VST DICT: %s' % (human_chrom, len(h5_mask_dict[human_chrom]['human_start']),len(matching_inds),len(vst_dict[human_chrom])))
        chimp_chroms = [chimp_chrom_list[c] for c in h5_mask_dict[human_chrom]['chimp_chrom']]
        print("LENGTH CHIMP CHROMS: ", len(chimp_chroms))
        
        for ind in matching_inds:
            human_start = h5_mask_dict[human_chrom]['human_start'][ind]
            human_end = h5_mask_dict[human_chrom]['human_end'][ind]
            chimp_start = h5_mask_dict[human_chrom]['chimp_start'][ind]
            chimp_end = h5_mask_dict[human_chrom]['chimp_end'][ind]
            chimp_chrom = chimp_chroms[ind]

            print('HUMAN CHROM: %s CHIMP CHROM: %s IND: %s ALL MATCHING INDS: %s' % (human_chrom,chimp_chrom,ind,len(matching_inds)))
            
            human_start_ind = int(np.floor(human_start/1000))
            human_end_ind = int(np.floor(human_end/1000))
            chimp_start_ind = int(np.floor(chimp_start/1000))
            chimp_end_ind = int(np.floor(chimp_end/1000)) + 1

            if human:
                start_ind = human_start_ind
                end_ind = human_end_ind
                chrom = human_chrom
            elif chimp:
                start_ind = chimp_start_ind
                end_ind = chimp_end_ind
                chrom = chimp_chrom

            if max_max:
                overall_list = []
                human_start_inds = np.arange(human_start_ind,human_end_ind,max_window)
                human_end_inds = human_start_inds + max_window
                if human:
                    start_chunks = human_start_inds
                    end_chunks = human_end_inds
                elif chimp:
                    if len(human_start_inds) == 1:
                        start_chunks = start_ind
                        end_chunks = end_ind
                    else:
                        start_chunks = np.around(np.linspace(chimp_start_ind,chimp_end_ind,len(human_start_inds)))[:-1]
                        end_chunks = np.around(np.linspace(chimp_start_ind,chimp_end_ind,len(human_start_inds)))[1:]
                print('HUMAN START IND: %s HUMAN END IND: %s START CHUNKS: %s END CHUNKS: %s' % (human_start_ind,human_end_ind,start_chunks,end_chunks))
                for i in np.arange(len(start_chunks)):
                    if vst:
                        overall_list.append(calc_vst(int(start_chunks[i]),int(end_chunks[i]),chrom,first_region,second_region,first_dict_cn,second_dict_cn,num_first,num_second))
                    elif variance:
                        overall_list.append(calc_variance(int(start_chunks[i]),int(end_chunks[i]),chrom,first_region,second_region,first_dict_cn,second_dict_cn))
                    elif bp_density:
                        overall_list.append(calc_bp_density(int(start_chunks[i]),int(end_chunks[i]),chrom,first_region,second_region,first_dict_cn,second_dict_cn,example))
                    elif pbs:
                        overall_list.append(calc_pbs(int(start_chunks[i]),int(end_chunks[i]),chrom,first_region,second_region,third_region,first_dict_cn,second_dict_cn,third_dict_cn,num_first,num_second,num_third))
                    elif dcgh:
                        overall_list.append(calc_dcgh(int(start_chunks[i]),int(end_chunks[i]),chrom,first_region,dcgh_ind,zero_replacement,first_dict_cn))
                vst_dict[human_chrom][ind] = max(overall_list)
            else:
                if vst:
                    vst_dict[human_chrom][ind] = calc_vst(start_ind,end_ind,chrom,first_region,second_region,first_dict_cn,second_dict_cn,num_first,num_second)
                elif variance:
                    vst_dict[human_chrom][ind] = calc_variance(start_ind,end_ind,chrom,first_region,second_region,first_dict_cn,second_dict_cn)
                elif bp_density:
                    vst_dict[human_chrom][ind] = calc_bp_density(start_ind,end_ind,chrom,first_region,second_region,first_dict_cn,second_dict_cn,example)
                elif pbs:
                    print('PBS!')
                    vst_dict[human_chrom][ind] = calc_pbs(start_ind,end_ind,chrom,first_region,second_region,third_region,first_dict_cn,second_dict_cn,third_dict_cn,num_first,num_second,num_third)
                elif dcgh:
                    vst_dict[human_chrom][ind] = calc_dcgh(start_ind,end_ind,chrom,first_region,dcgh_ind,zero_replacement,first_dict_cn)
            print('OUT VAL: ', vst_dict[human_chrom][ind])
        print('TIME FOR CHROM: ', (time.time()-t)/60.0)

    print('TIME FOR ALL CHROMS: ', (time.time()-overall_t)/60.0)

    out_file = h5py.File(out_name,'w')


    for human_chrom in list(h5_mask.keys()):
        chrom_grp = out_file.create_group(human_chrom)
        chrom_grp.create_dataset('human_start',data=h5_mask_dict[human_chrom]['human_start'])
        chrom_grp.create_dataset('human_end',data=h5_mask_dict[human_chrom]['human_end'])
        chrom_grp.create_dataset('chimp_start',data=h5_mask_dict[human_chrom]['chimp_start'])
        chrom_grp.create_dataset('chimp_end',data=h5_mask_dict[human_chrom]['chimp_end'])
        chrom_grp.create_dataset('chimp_chrom',data=h5_mask_dict[human_chrom]['chimp_chrom'])
        chrom_grp.create_dataset(out_grp,data=vst_dict[human_chrom])

    out_file.close()
    h5_mask.close()
    json_f.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--h5_mask_file',required=True,help='h5 mask file with alignment info')
    parser.add_argument('--chimp_chroms',required=True,help='Chimp chrom txt to decode .h5 file')
    parser.add_argument('--json_file','-j',required=True,help='Json file with regions')
    parser.add_argument('--first_region','-r1',required=False,default='',help='First region of interest')
    parser.add_argument('--second_region','-r2',required=False,default='',help='Second region of interest')
    parser.add_argument('--third_region','-r3',required=False,default='',help='Third region of interest - only used for PBS')
    parser.add_argument('--h5_path','-h5',required=True,help='Path to h5')
    parser.add_argument('--out_name','-o',help='Name of out file')
    parser.add_argument('--chimp',action='store_true')
    parser.add_argument('--human',action='store_true')
    parser.add_argument('--hmm',action='store_true')
    parser.add_argument('--vst',action='store_true')
    parser.add_argument('--variance',action='store_true')
    parser.add_argument('--bp_density',action='store_true')
    parser.add_argument('--pbs',action='store_true')
    parser.add_argument('--dcgh',action='store_true')
    parser.add_argument('--nobon',action='store_true')
    parser.add_argument('--bon',action='store_true')
    parser.add_argument('--dcgh_json',required=False,default='',help='Example individual for dCGH')
    parser.add_argument('--zero_replacement',required=False,default=0.1,help='What to replace zeros with in dcgh (no log of 0)')
    parser.add_argument('--max_max',action='store_true')
    parser.add_argument('--max_window',required=False,default=0,help='Only used for max max - the smaller window within the larger window to perform the action, in kb (ex. 10 for 10kb within 100kb)')
    o = parser.parse_args()
    vst_variance_bpdensity_pbs_dcgh(o)

