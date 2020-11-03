import h5py
import numpy as np
import os
import time
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.lines as lines
import json
import argparse
import pandas as pd

### CLUSTERING FUNCTIONS
class Tree:
    def __init__(self, id, indices=None, ro=False, left=False, right=False):
        self.id = id
        self.indices = indices
        self.ro = ro
        self.left  = left
        self.right = right

def make_tree_list(svs, ids): ### intitializes each sv call as a Tree instance
    tree_list = []
    for i in np.arange(len(svs)):
        tree_list.append(Tree(id=ids[i], indices=svs[i]))
    return tree_list

def overlap(indices1, indices2): ### distance function
    RO_per_pair = []
    pair_count = len(indices1) * len(indices2)
    for i in indices1:
        for j in indices2:
            S1, E1 = i[0], i[1]
            S2, E2 = j[0], j[1]
            L1, L2 = E1 - S1, E2 - S2
            O = max(0, min(E1, E2) - max(S1, S2))
            FO1, FO2 = O/L1, O/L2
            RO = np.mean([FO1, FO2])
            RO_per_pair.append(RO)
    similarity = sum(RO_per_pair)/pair_count
    return similarity

def cluster(tree_list, similarity_matrix=False):
    if len(tree_list) == 1: ### End Case
        return tree_list[0]
    else:
        ### On first iteration, complete similarity matrix is constructed
        if not similarity_matrix:
            similarity_matrix = []
            for i in tree_list:
                similarity_array = []
                for j in tree_list:
                    similarity = overlap(i.indices,j.indices)
                    similarity_array.append(similarity)
                similarity_matrix.append(similarity_array)
        ### Iterate through matrix to find pair with maximum reciprocal overlap
        max_index = [1,0]
        max_ro = 0
        length = np.arange(len(similarity_matrix))
        for i in length:
            for j in length:
                ro = similarity_matrix[i][j]
                if j != i and ro > max_ro:
                    max_ro = ro
                    max_index = [i,j]
        ### Identify corresponding nodes from tree_list and combine into new parent tree
        imax, jmax = max_index[0], max_index[1]
        combined_indices = tree_list[imax].indices + tree_list[jmax].indices
        combined_id = tree_list[imax].id + tree_list[jmax].id
        new_tree = Tree(id=combined_id, indices=combined_indices, ro=max_ro,
                        left=tree_list[imax], right=tree_list[jmax])
        ### Remove columns corresponding to clustered pair from tree_list and similarity_matrix
        ### Updating the similarity matrix instead of generating a new one lets this run MUCH faster
        print(imax, jmax, len(tree_list))
        similarity_matrix.pop(max(imax,jmax))
        similarity_matrix.pop(min(imax,jmax))
        for array in similarity_matrix:
            array.pop(max(imax,jmax))
            array.pop(min(imax,jmax))
        tree_list.pop(max(imax,jmax))
        tree_list.pop(min(imax,jmax))
        ### Add new parent tree to tree_list and similarity_matrix.
        tree_list.append(new_tree)
        similarity_matrix.append([])
        for i in np.arange(len(tree_list)):
            new_similarity = overlap(tree_list[i].indices, new_tree.indices)
            similarity_matrix[-1].append(new_similarity)
            similarity_matrix[i].append(new_similarity)
        similarity_matrix[-1].pop(-1)
        ### Recursively call cluster on updated tree_list and similarity_matrix
        cluster(tree_list, similarity_matrix)



### TREE TO JSON FUNCTIONS
# Recursively adds indices to list once a cluster is identified
def extract_cluster_indices(tree, indices):
    if len(tree.indices) == 1:
        indices.append(tree.indices[0])
    else:
        extract_cluster_indices(tree.left, indices)
        extract_cluster_indices(tree.right, indices)

# Identifies cluster in tree based on RO or terminal status and calls extract_cluster_indices()
def make_indices_list(tree, cutoff, indices_list):
    forms_cluster = (tree.ro > cutoff or tree.left==False)
    if forms_cluster:
        cluster_indices = []
        extract_cluster_indices(tree, cluster_indices)
        indices_list.append(cluster_indices)
    else:
        make_indices_list(tree.left, cutoff, indices_list)
        make_indices_list(tree.right, cutoff, indices_list)

# Builds the index_dict, which contains the SV name as the kay and average indices for each SV cluster as the value
# Also builds the helper_dict, which maps each individual SV indices back to its assigned SV
def make_index_and_helper_dict(indices_list, index_dict, helper_dict, svtype):
    avgs_list = []
    for cluster in indices_list:
        start_indices = [i[0] for i in cluster]
        end_indices = [i[1] for i in cluster]
        start_avg = np.mean(start_indices)
        end_avg = np.mean(end_indices)
        avgs_list.append([start_avg, end_avg])
    zipped_lists = zip(avgs_list, indices_list)
    sorted_pairs = sorted(zipped_lists, key=lambda x: x[0])
    tuples = zip(*sorted_pairs)
    avgs_list, indices_list = [ list(tuple) for tuple in  tuples]
    for i in np.arange(len(avgs_list)):
        key = svtype + '_' + str(i + 1)
        value = avgs_list[i]
        index_dict[key] = value
        for j in indices_list[i]:
            key2 = tuple(j)
            value2 = key
            helper_dict[key2] = value2

# Recursivley assigns SVs to individuals once a cluster is identified
def extract_samples(tree, sample_dict, helper_dict):
    if len(tree.id) == 1:
        key = tree.id[0]
        value = helper_dict[tuple(tree.indices[0])]
        if key in sample_dict.keys():
            sample_dict[key].append(value)
        else:
            sample_dict[key] = [value]
    else:
        extract_samples(tree.left, sample_dict, helper_dict)
        extract_samples(tree.right, sample_dict, helper_dict)

# Builds out sample dict by identifying cluster in tree based on RO or terminal status, then calls extract_samples()
def make_sample_dict(tree, cutoff, sample_dict, helper_dict):
    forms_cluster = (tree.ro > cutoff or tree.left==False)
    if forms_cluster:
        extract_samples(tree, sample_dict, helper_dict)
    else:
        make_sample_dict(tree.left, cutoff, sample_dict, helper_dict)
        make_sample_dict(tree.right, cutoff, sample_dict, helper_dict)



### TABLE FUNCTIONS
def fill_sample_df(tree, cutoff, sample_df, helper_dict, new_col1, new_col2):
    forms_cluster = (tree.ro > cutoff or tree.left==False)
    if forms_cluster:
        sample_df_helper(tree, sample_df, helper_dict, new_col1, new_col2)

    else:
        fill_sample_df(tree.left, cutoff, sample_df, helper_dict, new_col1, new_col2)
        fill_sample_df(tree.right, cutoff, sample_df, helper_dict, new_col1, new_col2)

def sample_df_helper(tree, sample_df, helper_dict, new_col1, new_col2):
    if len(tree.id) == 1:
        index = tree.indices[0]
        sv = helper_dict[tuple(index)]
        if not sample_df.loc[tree.id[0], new_col1]:
            sample_df.loc[tree.id[0], new_col1] = str(index)
            sample_df.loc[tree.id[0], new_col2] = sv
        else:
            sample_df.loc[tree.id[0], new_col1] += ', ' + str(index)
            sample_df.loc[tree.id[0], new_col2] += ', ' + sv
    else:
        sample_df_helper(tree.left, sample_df, helper_dict, new_col1, new_col2)
        sample_df_helper(tree.right, sample_df, helper_dict, new_col1, new_col2)



### TREE PLOTTING FUNCTIONS
color_list_raw = [[255,255,153],[255,128,0],[165,0,0],[255,0,0],[255,0,222],[198,255,219],[171,0,198],[242,148,255],
                [0,56,210],[129,239,255],[28,168,158],[104,229,200],[0,255,60],[0,177,27],[255,178,102],[236,207,250],
                [255,0,137],[192,192,0],[255,139,158],[205,75,0],[165,32,130],[53,114,0],[216,216,216],[159,159,159],
                [108,108,108],[250,216,78],[255,215,198],[153,0,76],[255,102,178],[102,0,102],
                [7.700000000000000000e+01,1.780000000000000000e+02,2.460000000000000000e+02]]

color_list = []
for color in color_list_raw:
    rgb = [number/255.0 for number in color]
    color_list.append(rgb)


def color_cycle(color_list):
    next = color_list.pop(0)
    color_list.append(next)
    return next

def plot_tree(tree, ax, window=[0,1], prior=False, cutoff=False, color='black'):
    if cutoff and prior:
        if tree.ro > cutoff and prior[2] < cutoff:
            color=color_cycle(color_list)
        elif tree.left == False and prior[2] < cutoff:
            color=color_cycle(color_list)
    if tree.left == False:
        y = np.mean(window)
        x=1.1
        ax.add_patch(matplotlib.patches.Circle((x, y), 0.005, facecolor=color))
        label = tree.id[0] + ';    ' + str(tree.indices[0])
        ax.text(x+0.01, y, label, fontsize=8)
    else:
        above = len(tree.left.id)
        below = len(tree.right.id)
        window_size = window[1]-window[0]
        y = (below/(above+below)) * window_size + window[0]
        x = tree.ro
        #ax.add_patch(matplotlib.patches.Circle((x, y), 0.005, facecolor=color))
        plot_tree(tree.left, ax, [y, window[1]], [x,y,tree.ro], cutoff, color)
        plot_tree(tree.right, ax, [window[0], y], [x,y,tree.ro], cutoff, color)
    if prior:
        horizontal = lines.Line2D([x,prior[0]], [y, y], color=color)
        vertical = lines.Line2D([prior[0], prior[0]], [y, prior[1]], color=color)
        ax.add_line(horizontal)
        ax.add_line(vertical)


### REGION INDICES (update when you want to use script on new regions)
region_dict = {'PSG':{'CHROM': 'chr19','START':42676412,'END':43311634},
               'GYP':{'CHROM': 'chr4','START':143832638,'END':144180633},
               'GSTM':{'CHROM': 'chr1','START':109672875,'END':109711913}}


### ALL THE FUNCTIONS ABOVE ARE CALLED BY hierarchical_cluster()
def hierarchical_cluster(args):
    #human, chimp, or ancient
    species = args.species
    #gstm, pbs, gyp, psg, q21, covid
    locus = args.locus
    #hmm, hmm_8states, pure_cn, sv10kb, sv25kb, sv50kb
    min_sv_size = args.min_sv_size
    #cluster cutoff
    cluster_cutoff = args.cluster_cutoff
    if cluster_cutoff == "all":
        cluster_cutoff = [0.75, 0.9, 0.95]
    else:
        cluster_cutoff = [float(cluster_cutoff)]
    #png_path
    png_path = args.png_path
    #csv_path
    csv_path = args.csv_path
    #json_path
    json_path = args.json_path
    #sv_type - either del or dup
    sv_type = args.sv_type
    #desired outputs
    outputs = args.outputs

    chrom = region_dict[locus]['CHROM']
    start = int(np.floor(region_dict[locus]['START']/1000.0))
    end = int(np.floor(region_dict[locus]['END']/1000.0))

    ### FETCHING DATA
    sv_list = []
    samp_list = []
    all_samples = []
    print('Fetching data. This may take a few minutes.')
    t = time.time()
    samps_with_svs = 0
    total_svs = 0
    h5_path = '/global/scratch2/almahalgren/' + species + '/' + min_sv_size + '_large_SVs/'
    h5s = os.listdir(h5_path)

    for h5 in h5s:
        all_samples.append(h5)
        h5_file_path = os.path.join(h5_path,h5)
        h5_file = h5py.File(h5_file_path,'r')
        samp_svs = np.where((h5_file['sv_masks'][chrom][sv_type + '_pairs'][...].T[0] >= start)*(h5_file['sv_masks'][chrom][sv_type + '_pairs'][...].T[1] <= end))[0]
        if len(samp_svs) > 0:
            samps_with_svs += 1
            total_svs += len(samp_svs)
            for i in np.array(h5_file['sv_masks'][chrom][sv_type + '_pairs'])[samp_svs]:
                sv_list.append(i)
                samp_list.append([h5])

        h5_file.close()

    for i in np.arange(len(sv_list)):
        sv_list[i] = [sv_list[i].tolist()]

    samp_svs = np.concatenate(sv_list)
    print('TIME TO LOAD SV PAIRS: %s MINUTES' % (np.around((time.time()-t)/60.0,decimals=2)))
    print('# SAMPLES WITH SVS: %s TOTAL SVS: %s' % (samps_with_svs,total_svs))

    ### PERFORMING HIERARCHICAL CLUSTERING
    print('Performing hierarchical clustering')
    tree_list = make_tree_list(sv_list, samp_list)
    t = time.time()
    print('Node_A', 'Node_B', 'Nodes_Remaining')
    cluster(tree_list)
    print('TIME TO RUN CLUSTERING: %s SECONDS' % (np.around((time.time()-t),decimals=2)))
    completed_tree = tree_list[0]

    ### GENERATING TREE FIGURES
    if 'p' in outputs:
        print('Generating Tree Figures')
        for cutoff in cluster_cutoff:
            fig_height = len(completed_tree.id)/4
            fig = plt.figure(figsize=(12,fig_height))
            ax = fig.add_subplot(111)
            plt.xlim([-.1,1.1])
            plt.xticks(np.arange(0, 1.1, 0.1))
            plt.xlabel('Reciprocal Overlap')
            plt.yticks([])
            ax.axvline(x=cutoff, color='black', ls="--")
            plot_tree(completed_tree, ax, cutoff=cutoff)
            fig_name = png_path + locus + '_' + sv_type + '_' + 'TreePlot' + '_' + 'cutoff=' + str(cutoff) + '.png'
            plt.savefig(fig_name, bbox_inches='tight')
        print('Tree Figures Complete')

    ### CREATING JSON FILES
    if 'j' in outputs:
        print('Creating JSON files')
        for cutoff in cluster_cutoff:
            indices_list = []
            helper_dict = {}
            index_dict = {}
            sample_dict = {}
            make_indices_list(completed_tree, cutoff, indices_list)
            make_index_and_helper_dict(indices_list, index_dict, helper_dict, sv_type)
            make_sample_dict(completed_tree, cutoff, sample_dict, helper_dict)
            output_json = {'INDEX_DICT':index_dict, 'SAMPLE_DICT':sample_dict}
            json_name = json_path + str.lower(locus) + '_' + species + '_' + sv_type + '_' + str(cutoff) + '_' + 'cluster.json'
            FOUT = open(json_name,'w')
            FOUT.write(json.dumps(json_name, indent=4, separators=(",", ": ")))
            FOUT.close()
        print('JSON files complete')

    ### CREATING A DATAFRAME OF SVs PER INDIVIDUAL PER CUTOFF & SAVING ANS A CSV
    if 'c' in outputs:
        print('Creating CSV file')
        sample_df = pd.DataFrame({'sample':all_samples})
        sample_df.set_index('sample', inplace=True)
        for cutoff in cluster_cutoff:
            indices_list = []
            index_dict = {}
            helper_dict = {}
            col1_name = locus + '_' + sv_type + '_' + 'indices' + '_' + 'cutoff='
            col2_name = locus + '_' + sv_type + '_' + 'SVs' + '_' + 'cutoff='
            new_col1  = col1_name + str(cutoff)
            new_col2 = col2_name + str(cutoff)
            sample_df[new_col1] = False
            sample_df[new_col2] = False
            make_indices_list(completed_tree, cutoff, indices_list)
            make_index_and_helper_dict(indices_list, index_dict, helper_dict, sv_type)
            fill_sample_df(completed_tree, cutoff, sample_df, helper_dict, new_col1, new_col2)
        csv_name = csv_path + locus + '_' + species + '_' + sv_type + str(cutoff) + 'cluster.csv'
        sample_df.to_csv(csv_name)
        print('CSV file complete')

    print('All complete')

    ### ArgumentParser

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--species','-s',required=True,help='human, chimp, or ancient')
    parser.add_argument('--locus','-l',required=True,help='psg, gyp, gstm, covid, none, or q21')
    parser.add_argument('--min_sv_size','-m',required=False,default='10kb',help='10kb, 25kb, or 50kb')
    parser.add_argument('--cluster_cutoff','-cc',required=False,default='all',help='Cutoff fraction on clustering tree. Fraction between 1 and 0, or "all" for [0.75, 0.9, 0.95]')
    parser.add_argument('--sv_type',required=True,help='either del or dup')
    parser.add_argument('--png_path','-p',required=False,default='/global/scratch/benjaminhyams/plotting_figures/cluster_tree_pngs/',help='path for Tree Plot')
    parser.add_argument('--json_path','-j',required=False,default='/global/scratch/benjaminhyams/plotting_figures/cluster_jsons/',help='path for json')
    parser.add_argument('--csv_path','-c',required=False,default='/global/scratch/benjaminhyams/plotting_figures/cluster_csvs/',help='path for csv file')
    parser.add_argument('--outputs','-o',required=False,default='pjc',help='any combo of p, j and c')
    o = parser.parse_args()
    hierarchical_cluster(o)
