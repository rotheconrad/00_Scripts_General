 #!/usr/bin/env python

'''Create stacked bar plot from data table is tsv format.

This script reads in a datatable in tab separated format and builds a
stacked bar plot output as a publication ready PDF.

Each column will be a separate bar on the x-axis and each
row value will be a separate stack in the bar along the y-axis.

Requires two column tab separated color table:

    row1 label  #FFFFFF
    row2 label  #FFF45E
    row3 label  #5EFFF3
    row4 label  #FF5EC5
    row5 label  #6A36A5

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Jan 2024
License :: GNU GPLv3
Copyright 2024 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd; import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import seaborn as sns
from statsmodels.stats.multitest import multipletests


def parse_colors(incolors):

    colors = {}

    with open(incolors, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            label = X[0]
            value = X[1]
            colors[label] = value

    return colors


def post_hoc_test(adf):
    ''' loops through individual rows and performs chi2 post hoc '''
    pvals = []

    bdf = adf.T
    for name in bdf.columns:
        xdf = bdf.drop(name, axis=1)
        xdf['OTHERs'] = xdf.sum(axis=1)
        xdf[name] = bdf[name]
        ydf = xdf[['OTHERs', name]].T
        c, p, d, x = stats.chi2_contingency(ydf, correction=True)
        # create expected frequency table
        extab = pd.DataFrame(x, index=ydf.index, columns=ydf.columns)
        extab['Total'] = extab.sum(axis=1)
        extab.loc['Total'] = extab.sum(axis=0)
        # create post hoc test contingency table
        ydf['Total'] = ydf.sum(axis=1)
        ydf.loc['Total'] = ydf.sum(axis=0)
        # print post hoc info
        print(f'\nPost hoc Chi2 test contingency table for {name}:\n')
        print(ydf)
        print(f'\nChi2 expected frequency table:\n\n{extab}')
        print(f'\nchi2 statistic: {c:.4f}, dof: {d}, chi2 pvalue: {p:.6f}')

        pvals.append(p)

        reject_list, pvals_corrected = multipletests(pvals, method='fdr_bh')[:2]

    return pvals, pvals_corrected


def chi2_hypothesis_test(adf):
    ''' performs chi square and post hoc tests between annotation
    categorgies for recombining vs non-recombining genes.
    '''
    
    # create and print contingency table
    ctab = adf.copy(deep=True)
    ctab['Total'] = ctab.sum(axis=1)
    ctab.loc['Total'] = ctab.sum(axis=0)
    #print(f'\nInitial Chi2 test contingency table:\n\n{ctab}')

    # chi2 test on the full data
    chi2, chip, dof, ex = stats.chi2_contingency(adf, correction=True)
    # create expected frequency table
    efreq = pd.DataFrame(ex, index=adf.index, columns=adf.columns)
    efreq['Total'] = efreq.sum(axis=1)
    efreq.loc['Total'] = efreq.sum(axis=0)
    print(f'\nChi2 expected frequency table:\n\n{efreq}')
    # print test statitic and p value
    print(f'\nchi2 statistic: {chi2:.4f}, dof: {dof}, chi2 pvalue: {chip:.6f}')

    # perform post hoc test on combinations if significant (< 0.05)
    if chip < 0.05:
        pvals, pvals_corrected = post_hoc_test(adf)
        print('\nPost hoc p values:\n', pvals)
        print('\nBenjamini/Hochberg corrected p values:\n', pvals_corrected)

    else:
        pvals_corrected = [1] * len(adf)

    return chip, pvals_corrected


def plot_stacked_barplot(df, title, colors, outfile, W, H):

    '''Plots annotation by group name'''

    print('\nBuilding annotation plot and tests ...')

    # categorical hypothesis testing with raw counts
    chip, pvals_corrected = chi2_hypothesis_test(df)

    # calculate the percents of total per category
    total = df.sum().to_list()
    ptots = [f'({round(i/sum(total) * 100, 2)}%)' for i in total]
    adf = df.div(df.sum(axis=0), axis=1)

    # initiate plot
    fig, ax = plt.subplots(figsize=(W,H))
    # plot data
    ax = adf.T.plot.bar(stacked=True, ax=ax, color=colors, width=.7)
    # set plot title
    ax.set_title(title)
    # change axis labels
    ax.set_xlabel('')
    ax.set_ylabel("Gene fraction", fontsize=12)
    # set the axis parameters / style
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=12)
    ax.tick_params(axis='x', labelrotation=45)
    ax.tick_params(axis='x', which='minor', bottom=False)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    ax.set_axisbelow(True)

    # annotate percent of total gene difference
    for i in range(len(ptots)):
        ax.annotate(ptots[i], (i, 1.02), transform=ax.transAxes, ha='center')
    

    # annotate individual percents
    for i, p in enumerate(ax.patches):
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        # don't print categories with 0 genes on the bar plot
        if round(height, 2) == 0:
            continue
        line = f'{height:.2f}'
        ax.text(x+width/2, 
                y+height/2, 
                line, 
                horizontalalignment='center', 
                verticalalignment='center')



    # add asterisk if significant post hoc
    # double the corrected p value array since our plot has two columns

    # create legend to separate file
    handles, labels = ax.get_legend_handles_labels()

    new_labels = []

    for i, l in enumerate(labels):
        nl = f'* {l} *' if pvals_corrected[i] <= 0.05 else l
        new_labels.append(nl)

    ax.legend(
                reversed(handles),
                reversed(new_labels),
                ncol=1,
                loc='center left',
                frameon=False,
                fontsize=12,
                bbox_to_anchor=(1, 0.5)
                )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    fig.savefig(outfile)
    plt.close() 
    
    return True


###############################################################################
##### MAIN MAIN MAIN MAIN MAIN MAIN MAIN ##############################
###############################################################################

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_datatable_tsv',
        help='Please specify the tsv formated data table.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--row_value_colors',
        help='Please specify the tsv formated color table.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name (Use .pdf).',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-t', '--plot_title',
        help='(OPTIONAL) Please specify the plot title (default: My Stacked Bar Plot).',
        metavar='',
        type=str,
        nargs='+',
        required=False,
        default=["My", "Stacked", "Bar", "Plot"]
        )
    parser.add_argument(
        '-x', '--figure_width',
        help='(OPTIONAL) Specify the figure width in inches.',
        metavar='',
        type=float,
        required=False,
        default=None
        )
    parser.add_argument(
        '-y', '--figure_height',
        help='(OPTIONAL) Specify the figure height in inches.',
        metavar='',
        type=float,
        required=False,
        default=None
        )
    args=vars(parser.parse_args())

    # define input params
    infile = args['input_datatable_tsv']
    outfile = args['output_file_name']
    title = ' '.join(args['plot_title'])
    incolors = args['row_value_colors']
    W = args['figure_width']
    H = args['figure_height']

    # Do what you came here to do:
    print('\n\nRunning Script ...')

    # read in the data table
    df = pd.read_csv(infile, sep='\t', index_col=0)

    colors = parse_colors(incolors)

    # define width and height based on number of columns and rows
    if not W:
        W = 2 * len(df.columns)

    if not H:
        H = 1 * len(df)

    # build the plotty plot
    _ = plot_stacked_barplot(df, title, colors, outfile, W, H)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
    