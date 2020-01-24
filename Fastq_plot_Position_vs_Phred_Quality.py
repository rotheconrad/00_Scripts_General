#!/usr/bin/env python

'''Evaluate Position(bp) vs Phred Score

Calculates the number of basepairs >= and < the user defined quality
threshold across the dataset for each bp position of the all reads.

Write 2 output files for plots as *.png:
1) * Heatmap of sequence position (x), quality (y), and count (heat).
   * Histogram of quality scores.
2) * Line plot of sequence position (x), and count (y) above Q score
   * Line plot of sequence position (x), and count (y) below Q score

!!! Requires Phred.tsv file - located in github repo !!!

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: November 22nd, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, os
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def make_phred_table(phred):
    d = {}
    with open(phred, 'r') as p:
        for l in p:
            X = l.rstrip().split('\t')
            Q = int(X[0])
            C = X[1]
            d[C] = Q

    return d

def evaluate_data(fdir, phred, Q): 

    daboveQ = defaultdict(int)
    dbelowQ = defaultdict(int)
    dAll = defaultdict(lambda: defaultdict(int))
    Qcounts = defaultdict(int)
    total = 0
    totalbelowQ = 0
    totalaboveQ = 0

    f_list = [f for f in os.listdir(fdir) if os.path.isfile(f'{fdir}/{f}')]

    for fastq in f_list:
        line_count = 0

        with open(f'{fdir}/{fastq}', 'r') as f:
            for l in f:
                line_count += 1
                if line_count%4 == 0:
                    seq = l.rstrip()
                    for i, b in enumerate(seq):
                        i = i+1
                        score = phred[b]
                        total += 1
                        dAll[i][score] += 1
                        Qcounts[score] += 1

                        if score < Q:
                            dbelowQ[i] += 1
                            totalbelowQ += 1

                        elif score >= Q:
                            daboveQ[i] += 1
                            totalaboveQ += 1

                        else:
                            print(b)

    percentQ = totalbelowQ / total * 100

    totals = [
                f'Total Base Pairs: {total}',
                f'Total BPs >= Q{Q}: {totalaboveQ}',
                f'Total BPs < Q{Q}: {totalbelowQ}',
                f'Percent BPs < Q{Q}: {percentQ:.2f}%'
                ]

    return daboveQ, dbelowQ, dAll, Qcounts, totals

def plot_data(daboveQ, dbelowQ, totals, Q, out):

    # Set Colors
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    
    # prepare to plot data
    xBQ, xAQ = [], []
    yBQ, yAQ = [], []

    for i in range(1,152):
        xBQ.append(i)
        xAQ.append(i)
        yBQ.append(dbelowQ.get(i, 0))
        yAQ.append(daboveQ.get(i,0))

    # Build the Plot
    fig, (ax1, ax2) = plt.subplots(
                        2, 1, figsize=(20,10), sharex=True, sharey=False
                        )

    # Plot titles and labels
    ax1.set_title(f'Sequence Position vs bps above Q{Q}')
    ax1.set_xlabel('Sequence Position')
    ax1.set_ylabel(f'Number of bps above Q{Q}')

    ax2.set_title(f'Sequence Position vs bps equal to Q{Q}')
    ax2.set_xlabel('Sequence Position')
    ax2.set_ylabel(f'Number of bps equal to Q{Q}')
    ax2.text(
        0.50, 0.98, '\n'.join(totals), 
        fontsize=18, color='#a50f15',
        verticalalignment='top', horizontalalignment='right',
        transform=ax2.transAxes
        )

    # Setup up plot grids, spines and ticks
    for ax in fig.axes:
        ax.minorticks_on()
        ax.tick_params(
            which='minor', axis='both', left=False, bottom=False
            )
        ax.tick_params(
                    which='major', axis='both',
                    left=False, bottom=True,
                    size=8, width=5, tickdir='in',
                    labelsize=16, zorder=10
                    )
        ax.yaxis.grid(
            which="minor", color=gridm, linestyle='--',
            linewidth=1, alpha=0.6, zorder=1
            )
        ax.yaxis.grid(
            which="major", color=gridM, linestyle='--',
            linewidth=1.5, alpha=0.4, zorder=1
            )
        for spine in ax.spines.values(): spine.set_linewidth(2)
        ax.set_axisbelow(True)

    # Plots
    ax1.plot(xAQ, yAQ)
    ax2.plot(xBQ, yBQ)

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{out}_lines.png')
    plt.close()


def plot_all(dAll, Qcounts, out):

    # Set Colors
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    mc = '#a50f15'
    
    # prepare to plot data
    All = {}
    Qx = []
    Qy = []

    for col in range(1,152):
        vals = []
        for v in range(0,43):
            vals.append(dAll[col].get(v, 0))
        All[col] = vals

    df = pd.DataFrame.from_dict(All)
    #df = df.div(df.to_numpy().sum())

    for i in range(42, -1, -1):
        Qx.append(i)
        Qy.append(Qcounts.get(i, 0))


    # Build the Plot
    fig, (ax1, ax2) = plt.subplots(
                        2, 1, figsize=(20,10), sharex=False, sharey=False
                        )

    # Setup up plot grids, spines and ticks
    ax1.minorticks_on()
    ax1.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax1.tick_params(
                which='major', axis='both',
                left=False, bottom=True,
                size=8, width=5, tickdir='in',
                labelsize=16, zorder=10
                )
    ax1.yaxis.grid(
        which="minor", color=gridm, linestyle='--',
        linewidth=1, alpha=0.6, zorder=1
        )
    ax1.yaxis.grid(
        which="major", color=gridM, linestyle='--',
        linewidth=1.5, alpha=0.4, zorder=1
        )
    for spine in ax1.spines.values(): spine.set_linewidth(2)
    ax1.set_axisbelow(True)

    # Plots
    ax1 = sns.heatmap(df, ax=ax1)
    ax1.invert_yaxis()
    ax2.bar(Qx, Qy, width=0.8, align='center')

    # Plot titles and labels
    ax1.set_title('Sequence Position vs Quality Score')
    ax1.set_xlabel('Sequence Position')
    ax1.set_ylabel('Quality Score')

    ax2.set_title(f'Histogram of Quality Scores')
    ax2.set_xlabel('Quality Score')
    ax2.set_ylabel('Counts')

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{out}_HeatmapHist.png')
    plt.close()


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-d', '--file_directory',
        help='Please specify the directory name of fastq files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='Please specify the output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-q', '--Phred_Quality_Threshold',
        help='Please specify Phred Quality (ie: 20)!',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-p', '--Phred_tsv_file',
        help='Please specify the Phred.tsv file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Retrieve Phred Score ASCII characters
    phred = make_phred_table(args['Phred_tsv_file'])

    # check fdir path does not end with '/'
    fdir = args['file_directory']
    if fdir[-1] == '/': fdir = fdir[:-1]

    daboveQ, dbelowQ, dAll, Qcounts, totals, = evaluate_data(
                                            fdir,
                                            phred,
                                            args['Phred_Quality_Threshold'],
                                            )

    _ = plot_data(
                daboveQ,
                dbelowQ,
                totals,
                args['Phred_Quality_Threshold'],
                args['out_file']
                )

    _ = plot_all(dAll, Qcounts, args['out_file'])

if __name__ == "__main__":
    main()