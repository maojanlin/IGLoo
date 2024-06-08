import pandas as pd
import matplotlib.pyplot as plt
import argparse

import seaborn as sns
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap


def plot_asm(gene_report, sample_id, output_name):
    my_colors = ['#000000', '#751C1C', '#E93838', '#F49A90', '#FACBBC', '#FFECDA', '#FFFCE8']
    my_cmap = ListedColormap(my_colors)
    bounds = [0, 0.5, 1, 1.5, 2, 3, 10, 100]
    my_norm = BoundaryNorm(bounds, ncolors=len(my_colors))
    
    ###################################### HAP1 ######################################
    df_sample = pd.read_csv(gene_report, sep=',', header=0, index_col=0)
    
    grid_kws = {"width_ratios": (.5, .01), "hspace": .1}
    fig, (ax, cbar_ax) = plt.subplots(ncols=2, figsize=(25,df_sample.shape[0]/1.5), gridspec_kw=grid_kws)
    plt.subplots_adjust(wspace=0.05)
    sns.heatmap(df_sample,
            ax=ax,
            cmap=my_cmap,
            norm=my_norm,
            cbar_ax=cbar_ax,
            cbar_kws={"orientation": "vertical"})
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    ax.invert_xaxis()
    title = ax.set_title(sample_id, y=1.08)
    title.set_position([0.5, 1.05])
    
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([(b0+b1)/2 for b0, b1 in zip(bounds[:-1], bounds[1:])])
    colorbar.set_ticklabels(['0', '<=0.5', '<=-1', '<=1.5', '<=2', '<=3', '>3'])
    
    plt.savefig(output_name, bbox_inches='tight')



def main(arguments=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-csv', '--contig_genes_report', help='the report listing number of called IG genes.', required=True)
    parser.add_argument('-id',  '--sample_ID', help='sample ID', default="sample", required=False)
    parser.add_argument('-out', '--output_name', help='output file name')
    args = parser.parse_args(arguments)
    
    gene_report = args.contig_genes_report
    sample_id = args.sample_ID
    output_name = args.output_name
    if output_name == None:
        output_name = gene_report + '.pdf'

    plot_asm(gene_report, sample_id, output_name)



if __name__ == "__main__":
    main()
