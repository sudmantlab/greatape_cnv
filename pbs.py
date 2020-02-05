import h5py
import sys
import time
import numpy as np
import argparse
import os

#PBS: calculates the PBS values per sliding window in the three input branches
#Inputs: path to folder with the .h5 CN files of subspecies samples
#Outputs: .h5 files with the variance of individual subspecies as well as joint variances;
#.h5 file with /chrom/{Start, End, 3 branches}



def pbs(args):
    
    first_name = args.first_name
    first_path = args.first_path
    second_name = args.second_name
    second_path = args.second_path
    third_name = args.third_name
    third_path = args.third_path
    output_name = args.output_name

    print('First path: ', first_path)
    print('Second path: ', second_path)
    print('Third path: ', third_path)

    first_second_name = first_name + '_' + second_name
    first_third_name = first_name + '_' + third_name
    second_third_name = second_name + '_' + third_name

    t = time.time()

    pop_variance(first_name, first_path)
    pop_variance(second_name, second_path)
    pop_variance(third_name, third_path)


    first_var = h5py.File(first_name + '_variance.h5', 'r')
    second_var = h5py.File(second_name + '_variance.h5', 'r')
    third_var = h5py.File(third_name + '_variance.h5', 'r')
    print('Done with single variances! This took: ', time.time()-t)

    num_first = first_var['Num_Samples']['Num'][()]
    num_second = second_var['Num_Samples']['Num'][()]
    num_third = third_var['Num_Samples']['Num'][()]

    print('First num: ', num_first)
    print('Second num: ', num_second)
    print('Third num: ', num_third)

    t = time.time()


    pop_variance(first_second_name, first_path + second_path)
    pop_variance(first_third_name, first_path + third_path)
    pop_variance(second_third_name, second_path + third_path)
    
    first_second_var = h5py.File(first_second_name + '_variance.h5', 'r')
    first_third_var = h5py.File(first_third_name + '_variance.h5', 'r')
    second_third_var = h5py.File(second_third_name + '_variance.h5', 'r')
    print('Done with double variances! This took: ', time.time()-t)

    out_file = h5py.File(output_name, 'w')

    #first_out = h5py.File(first_name + '_pbs.h5', 'w')
    #second_out = h5py.File(second_name + '_pbs.h5', 'w')
    #third_out = h5py.File(third_name + '_pbs.h5', 'w')

    t_overall = time.time()

    for chrom in list(first_var['Variances'].keys()):
        
        t = time.time()

        var1 = np.array(first_var['Variances'][chrom]['Variance'][...])
        var2 = np.array(second_var['Variances'][chrom]['Variance'][...])
        var3 = np.array(third_var['Variances'][chrom]['Variance'][...])
        var12 = np.array(first_second_var['Variances'][chrom]['Variance'][...])
        var13 = np.array(first_third_var['Variances'][chrom]['Variance'][...])
        var23 = np.array(second_third_var['Variances'][chrom]['Variance'][...])

        grp = out_file.create_group(chrom)

        if (chrom == 'chr1'):
            print('Loading variance arrays for chr1: ', time.time()-t)
            print(len(var1))
            print(len(var2))
            print(len(var3))
            print(len(var12))
            print(len(var13))
            print(len(var23))

        #first_grp = first_out.create_group(chrom)
        #second_grp = second_out.create_group(chrom)
        #third_grp = third_out.create_group(chrom)

        #first_pbs = [vst(var1[i], num_first, var

        t = time.time()
    
        vst12 = vst(var1, num_first, var2, num_second, var12, num_first+num_second)
        vst13 = vst(var1, num_first, var3, num_third, var13, num_first+num_third)
        vst23 = vst(var2, num_second, var3, num_third, var23, num_second+num_third)

        pbs1 = 0.5*(vst12 + vst13 - vst23)
        pbs2 = 0.5*(vst12 + vst23 - vst13)
        pbs3 = 0.5*(vst13 + vst23 - vst12)

        if (chrom == 'chr1'):
            print('Time to calculate pbs_vals and slice: ', time.time()-t)
            print('Length of PBS1: ', len(pbs1))
            print('First few values of PBS1: ', pbs1[0:100])
            print('Length of PBS2: ', len(pbs2))
            print('Length of PBS3: ', len(pbs3))

        start = grp.create_dataset('Start', dtype='uint64', data=first_var['Variances'][chrom]['Start'][...])
        end = grp.create_dataset('End', dtype='uint64', data=first_var['Variances'][chrom]['End'][...])
        first_pbs = grp.create_dataset(first_name, data=pbs1)
        second_pbs = grp.create_dataset(second_name, data=pbs2)
        third_pbs = grp.create_dataset(third_name, data=pbs3)


    print("Time to calculate pbs for all chroms: ", time.time()-t_overall)

    out_file.close()


def pbs_calc(v1, n1, v2, n2, v3, n3, v12, v13, v23):
    
    vst12 = vst(v1, n1, v2, n2, v12, n1+n2)
    vst13 = vst(v1, n1, v3, n3, v13, n1+n3)
    vst23 = vst(v2, n2, v3, n3, v23, n2+n3)

    pbs1 = 0.5*(vst12 + vst13 - vst23)
    pbs2 = 0.5*(vst12 + vst23 - vst13)
    pbs3 = 0.5*(vst13 + vst23 - vst12)

    return np.array([pbs1, pbs2, pbs3])


def num_samples(*args):
    cwd = os.getcwd()
    samp_list = []

    for m in args:
        dirs = os.listdir(m)
        for d in dirs:
            samp = os.path.join(m, d)
            if (samp[-3:] == '.h5'):
                samp_list.append(d[-3:])
    return len(samp_list)                


def vst(v1, n1, v2, n2, v_tot, n_tot):
    zeroes = np.where(v_tot == 0)[0]
    if (len(zeroes) > 0):
        print('THERE ARE ZEROES IN V_TOT!')
        print('Zero indices: ', zeroes)
        v_tot[zeroes] = 0.0001
        #return 1
    result = 1 - (((v1*n1) + (v2*n2))/(v_tot*n_tot))
    result[zeroes] = 0
    return result



def pop_variance(output_name, path_list):
    out_file = h5py.File(output_name + '_variance.h5', "w")
    num_group = out_file.create_group('Num_Samples')
    var_group = out_file.create_group('Variances')


    cwd = os.getcwd()
    samp_list = []

    #print("input mappings_path", mappings_path)
    for m in path_list:
        print("First path: ", m)
        dirs = os.listdir(m)
        for d in dirs:
            #print("SAMPLE NAME: ", d)
            samp = os.path.join(m, d)
            if (samp[-3:] == ".h5"):
                samp_list.append(d[:-3])
                print("Name used: ", samp_list[-1])
                samp_list[-1] = h5py.File(samp, "r")

    print(samp_list)
    print("LENGTH OF SAMP_LIST: ", len(samp_list))

    num_group.create_dataset('Num', data=len(samp_list))

    chroms = list(samp_list[0]['depth'].keys())
    print("NUMBER OF CHROMOSOMES: ", len(chroms))
    

    for chr in chroms:
     #   print(chr)
        #concatenates arrays alpng the long axis
        #ex: [[1,2,3],[4,5,6]] -> [[1,4],[2,5],[3,6]]
        #putting all the cnv values from each position into one array in order to calculate the variance
        vals = [samp['depth'][chr]['copy_number'][...] for samp in samp_list]
        variances = np.var(vals, axis=0)
        """
    #tot = np.c_[[samp['depth'][chr]['copy_number'][...] for samp in samp_list]]
    #if (chr == 'chr1'):
    print("IN CHR 1")
    print('Length of copy number arrays: ', len(vals[0]))
    print('Length of copy number arrays: ', len(vals[8]))
    print('beginning of array with all hte cn: ', len([v[0:30] for v in vals]))
    print('the vals array: ', vals)
    print('Length of variances: ', len(variances))
    print('Wrong way for variances: ', len(np.var(vals, axis=1)))
      
        `
        #CHANGED: variance from 0 to 1 - calculating the variance per index
        variances = np.var(tot, axis=1)
        print("First row: ", tot[0:,0:10])
        print("First variances: ", variances[0:10])
        """    
        grp = var_group.create_group(chr)

        start = samp_list[0]['depth'][chr]['start'][...]
        end = samp_list[0]['depth'][chr]['end'][...]

        if (len(start) != len(variances)):
            print("SOMETHING IS WRONG")

        start_dset = grp.create_dataset('Start', dtype='uint64', data=start)
        end_dset = grp.create_dataset('End', dtype='uint64', data=end)
        var_dset = grp.create_dataset('Variance', data=variances)
    

    out_file.close()

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--first_name", required=True, help="Group 1 name e.g. verus_ellioti")
    parser.add_argument("--first_path", required=True, nargs='+', help="Path to group 1: ../verus ../ellioti")
    parser.add_argument("--second_name", required=True, help="Group 2 name e.g. schweinfurthii_troglodytes")
    parser.add_argument("--second_path", required=True, nargs='+', help="Path to group 2: ../schweinfurthii ../troglodytes")
    parser.add_argument("--third_name", required=True, help="Group 3 name e.g. bonobos")
    parser.add_argument("--third_path", required=True, nargs='+', help="Path to group 3: ../bonobos")
    parser.add_argument("--output_name", required=True, help="chimps_bonobos_pbs.hdf5")
    o = parser.parse_args()

    pbs(o)
    


