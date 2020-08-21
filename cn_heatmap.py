import numpy as np
#from hmmlearn import hmm
import matplotlib
import matplotlib.pyplot as plt
import h5py
import os
import time
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#import pandas
from scipy.spatial.distance import pdist, cdist,squareform
from matplotlib.ticker import MaxNLocator
import json
import argparse

"""
cn_heatmap returns a png figure with a heatmap of the population in the desired region.

inputs:
    start: start coordinate in chrom, ex. 1400000
    end: end coordinate in chrom, ex. 1800000
    chrom: chromosome, ex. 'chr1'
    title: Name of locus, ex. PSG or GYP
    png_path: path to output folder for png file
    humans, ancients, chimps: True/false for each
    hmm: true/false

if gene_json_file:
    config json file with the following info:
    start = 40687873 - 15000
    end = 41186467 + 15000
    chrom = 'chr19'
    pos = str(chrom) + ':' + str(start) + '-' + str(end)
    suptitle = 'Chimp PSG region: (' + pos + ')'
    gap = ''
    gene_dict = {'PSG1':'chr19:40864573-40877939','PSG2': 'chr19:40961042-40977614','PSG3': 'chr19:40702873-40719468','PSG4': 'chr19:41094922-41108008','PSG5':'chr19:41072271-41088758','PSG6': 'chr19:41160100-41171467','PSG7': 'chr19:40756981-41107872','PSG8': 'chr19:40731646-40744651','PSG9': 'chr19:41155284-41171467','PSG11': 'chr19:40904234-40920781'}
    top_tier = ['PSG3','PSG6','PSG11','PSG5']
    bottom_tier = ['PSG8','PSG1','PSG7','PSG2','PSG4','PSG9']
    gene_order = ['PSG1','PSG2','PSG3','PSG4','PSG5','PSG6','PSG7','PSG8','PSG9','PSG11']
    num_genes = len(gene_order)
    geneblocks_colors = ['xkcd:wine']

outputs:
    .png file

"""

def cn_heatmap(args):
    gene_json_file = args.gene_json_file
    title = str(args.title)
    png_path = str(args.png_path)
    humans = args.humans
    ancients = args.ancients
    chimps = args.chimps
    hmm = args.hmm

    window = 1000
    rgb_factor = 1/255.0

    if (humans):
        h5_path = '/global/scratch2/almahalgren/humans/repmasked_smoothed_outputs/'
        reg_json_file = '/global/scratch2/almahalgren/humans/configs/human_heatmap_config.json'
        trace_per_ind = 0.1
        fig_len = 15
        png_name  = 'humans_' + str(title)
    elif (ancients):
        h5_path = '/global/scratch2/almahalgren/ancient_dna/hg38_data/repmasked_smoothed_outputs/'
        reg_json_file = '/global/scratch2/almahalgren/ancient_dna/hg38_data/configs/ancient_heatmap_config.json'
        trace_per_ind = 0.25
        fig_len = 12
        png_name  = 'ancients_' + str(title)
    elif (chimps):
        h5_path = '/global/scratch2/almahalgren/chimps/repmasked_smoothed_outputs/'
        reg_json_file = '/global/scratch2/almahalgren/chimps/configs/chimp_heatmap_config.json'
        trace_per_ind = 0.25
        fig_len = 8
        png_name  = 'chimps_' + str(title)

    if (hmm):
        cn_color_file = '/global/scratch2/almahalgren/chimps/cnv/hmm_cn_colors.txt'
    else:
        cn_color_file = '/global/scratch2/almahalgren/chimps/cnv/pure_cn_colors.txt'

    #COLORS-----------------------------------

    #LOADING COLORS FROM FILE:
    color_array = np.loadtxt(cn_color_file)
    #Each color needs to be a fraction out of 255.0
    color_array = [rgb_factor*np.array(col) for col in color_array]
    colormaps = [ListedColormap(color_array)]
    n = len(colormaps)
    max_colors = len(color_array)
    print('MAX COLORS: ', max_colors)

    
    with open(reg_json_file) as reg_json_f:
        reg_json = json.load(reg_json_f)
    
    reg_list = reg_json['REGIONS'].keys()
    nice_titles = reg_json['TITLE_ABBREV']
    color_dict = reg_json['COLOR_DICT']
    
    if gene_json_file:
        with open(gene_json_file) as gene_json_f:
            gene_json = json.load(gene_json_f)
        start = gene_json['START'] - 15000
        end = gene_json['END'] + 15000
        chrom = gene_json['CHROM']
        gene_dict = gene_json['GENE_DICT']
        top_tier = gene_json['TOP_TIER']
        bottom_tier = gene_json['BOTTOM_TIER']
        gene_order = gene_json['GENE_ORDER']
        num_genes = len(gene_order)
        geneblocks_colors = gene_json['GENE_COLORS']
        title = gene_json['TITLE']
    else:
        start = int(args.start) - 15000
        end = int(args.end) + 15000
        chrom = str(args.chrom)
        title = str(args.title)

    pos = str(chrom) + ':' + str(start) + '-' + str(end)
    suptitle = title + ': (' + pos + ')'

    start_bp = int(np.around(start/window)*window)
    end_bp = int(np.around(end/window)*window)

    start_cn = int(start_bp/window)
    end_cn = int(end_bp/window)

    print(start_cn)
    print(end_cn)
    
    cn_dict = {}
    
    h5s = os.listdir(h5_path)
    for i in h5s:
        samp = os.path.join(h5_path,i)
        if (hmm):
            truth_statement = (i[-6:] == 'hmm.h5')
        else:
            truth_statement = (i[-3:] == '.h5' and i[-6:] != 'hmm.h5')
        if (truth_statement):
            h5_file = h5py.File(samp,'r')
            if (hmm):
                cn_dict[i[:-9]] = {}
                lowess = h5_file[chrom]['cn_hmm'][start_cn:end_cn]
                lowess[np.where(lowess >= max_colors)[0]] = max_colors
                cn_dict[i[:-9]]['LOWESS'] = lowess
            else:
                if ('LOWESS_smoo_med_cn' not in h5_file['depth'][chrom].keys()):
                    print('SAMPLE WITHOUT LOWESS_smoo_med_cn: ', i)
                else:
                    cn_dict[i[:-3]] = {}
                    lowess = np.around(h5_file['depth'][chrom]['LOWESS_smoo_med_cn'][start_cn:end_cn]).astype(int)
                    lowess[np.where(lowess >= max_colors)[0]] = max_colors
                    cn_dict[i[:-3]]['LOWESS'] = lowess
            #lowcut = np.around(h5_file['depth'][chrom]['LOWESS_cutoff_smoo_med_cn'][start_cn:end_cn]).astype(int)
            #lowcut[np.where(lowcut >= max_colors)[0]] = max_colors
            #cn_dict[i[:-3]]['LOWCUT'] = lowcut
            h5_file.close()
        
    print('LENGTH OF CN DICT: ', len(cn_dict))
    print('LENGTH OF CN KEYS: ', len(cn_dict.keys()))
          
    reg_dict = {}
    heights = []
    for reg in reg_list:
        reg_dict[reg] = {}
        reg_dict[reg]['CNs'] = []
        heights.append(len(reg_json[reg])*trace_per_ind)
        reg_dict[reg]['Height'] = len(reg_json[reg])*trace_per_ind
 
        for ind in reg_json[reg]:
            if (ind in cn_dict.keys()):
                reg_dict[reg]['CNs'].append(list(cn_dict[ind]['LOWESS'][...]))
    first_heights = max(heights)
    if (gene_file):
        heights.insert(0,first_heights*0.5)
    """
    if (pbs):
        heights.append(first_heights*1.2)
        heights.append(first_heights*1.2)
        heights.append(first_heights*1.2)
        heights.append(first_heights*1.2)
    #heights.append(first_heights*1.2)
    """
    print('HEIGHTS: ', heights)
    
    fig = plt.figure(figsize=(12,fig_len))
    widths = [14,0.5]
    spec = fig.add_gridspec(ncols=2, nrows=len(heights), width_ratios=widths,height_ratios=heights)
    gene_row = 0
    if (gene_json_file):
        ax0 = fig.add_subplot(spec[0,0])
        color_ind = 0
        gene_height = 200
        gene_slot = 400
        bottom_cushion = 200
        total_y = 1600
        width_bp = end-start
        total_x = 400
        for i,gene in enumerate(gene_order):
            if (gene != ''):
            #if len(gene_dict[gene] > 1):
                if (True):
                    start_x_bp = int(gene_dict[gene].split(':')[1].split('-')[0])
                    end_x_bp = int(gene_dict[gene].split(':')[1].split('-')[1])
                    start_x_rect = ((start_x_bp - start)/width_bp)*total_x
                    end_x_rect = ((end_x_bp - start)/width_bp)*total_x
                    if (gene in bottom_tier):
                        start_y_rect = 200
                    elif (gene in top_tier):
                        start_y_rect = 800
                    rect_width = end_x_rect - start_x_rect
                    ax0.add_patch(matplotlib.patches.Rectangle((start_x_rect,start_y_rect), rect_width, gene_height, color=geneblocks_colors[color_ind]))
                    if (humans and psg):
                        ax0.text(start_x_rect + (rect_width/2) - 8, start_y_rect + gene_height + 18, gene, fontsize = 12)
                    else:
                        ax0.text(start_x_rect + (rect_width/2) - 15, start_y_rect + gene_height + 18, gene, fontsize = 12)
                color_ind = 0
        #plt.tick_params(axis='y',which='both',bottom=False,top=False,labelbottom=False)
        ax0.set_xlim([0, total_x])
        ax0.set_ylim([0, total_y])
        ax0.set_xticks([])
        ax0.set_yticks([])
        ax0.spines['bottom'].set_color('white')
        ax0.spines['top'].set_color('white') 
        ax0.spines['right'].set_color('white')
        ax0.spines['left'].set_color('white')
        
        print('SPEC[0,1]:',spec[0,1])
        ax1 = fig.add_subplot(spec[0,1])
        ax1.set_xlim([0,10])
        ax1.set_ylim([0,10])
        col = 'white'
        ax1.add_patch(matplotlib.patches.Rectangle((0,0),10,10,color=col))
        ax1.spines['bottom'].set_color(col)
        ax1.spines['top'].set_color(col) 
        ax1.spines['right'].set_color(col)
        ax1.spines['left'].set_color(col)
        ax1.set_xticks([])
        ax1.set_yticks([])
        gene_row = 1
    
    
    for row in np.arange(num_regions):   
        ax0 = fig.add_subplot(spec[row+gene_row,0])
        reg = reg_list[row]
        col = color_dict[reg]
        psm = ax0.pcolormesh(reg_dict[reg]['CNs'], cmap=colormaps[0], rasterized=True, vmin=0, vmax=max_colors)
        if (row == num_regions-1):
            num_ticks = 5    
            tick_spread = np.linspace(0,len(reg_dict[reg]['CNs'][0]),num_ticks)
            ax0.set_xticks(tick_spread)
            ax0.set_xticklabels([str(str(np.around(i,decimals=3)) + ' Mb') for i in np.linspace(start/10**6,end/10**6,num_ticks + 1)],fontsize=12)
        else:
            ax0.set_xticks([])
        if (row == 2):
            ax0.set_ylabel('Copy Number',fontsize=14)
        ax0.set_yticks([])
        #plt.title(nice_titles[reg],fontsize=15)

        ax1 = fig.add_subplot(spec[row+1,1])
        ax1.set_xlim([0,10])
        ax1.set_ylim([0,10])
        ax1.add_patch(matplotlib.patches.Rectangle((0,0),10,10,color=col))
        if (chimps):
            ax1.set_ylabel(nice_titles[reg],fontsize=14,rotation=90)
            #ax1.set_ylabel(nice_titles[reg],fontsize=15,style='italic',rotation='horizontal')
        else:
            ax1.set_ylabel(nice_titles[reg],fontsize=14,rotation=90)
            #ax1.set_ylabel(nice_titles[reg],fontsize=15,rotation='horizontal')
        ax1.yaxis.set_label_position('right')
        ax1.spines['bottom'].set_color(col)
        ax1.spines['top'].set_color(col) 
        ax1.spines['right'].set_color(col)
        ax1.spines['left'].set_color('black')
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        
    plt.suptitle(suptitle,fontsize=15,y=1)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig(png_path + png_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", '-s', required=False,help="Start index")
    parser.add_argument("--end", '-e', required=False,help="End index")
    parser.add_argument("--chrom", '-c', required=False,help="Chromosome")
    parser.add_argument("--gene_json_file", '-g', required=False,default='',help="json file with info about genes in this region ex. human_regions.json")
    parser.add_argument('--title','-t',required=True,help='Title of the region ex. PSG, GYP')
    parser.add_argument('png_path','-png',required=True,help='png folder to write heatmap: /global/scratch2/almahalgren/humans/figures/')
    parser.add_argument('--humans',action='store_true')
    parser.add_argument('--ancients',action='store_true')
    parser.add_argument('--chimps',action='store_true')
    parser.add_argument('--hmm',action='store_true')
    o = parser.parse_args()
    cn_heatmap(o)
    
