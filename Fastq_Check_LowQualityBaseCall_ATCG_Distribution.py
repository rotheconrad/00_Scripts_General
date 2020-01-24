#!/usr/bin/env python

'''Low Quality Base Call ATCGN Distribution

Reads through directory of fastq files and counts the total number and
the number <= the user defined quality score of ATCGN in the dataset.

Writes bar plots to outfile.

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
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def make_phred_table(phred):
    d = {}
    with open(phred, 'r') as p:
        for l in p:
            X = l.rstrip().split('\t')
            Q = int(X[0])
            C = X[1]
            d[C] = Q

    return d

def count_ATCGN(fdir, phred, Q):
    
    total = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
    low = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}

    f_list = [f for f in os.listdir(fdir) if os.path.isfile(f'{fdir}/{f}')]

    for fastq in f_list:

        line_count = 0
        seq, qal = None, None

        with open(f'{fdir}/{fastq}', 'r') as f:
            for l in f:
                line_count += 1

                if line_count%4 == 0:
                    qal = l.rstrip()

                elif line_count%2 == 0:
                    seq = l.rstrip().upper()
                    for b in seq: total[b] += 1

                if qal:
                    for i, b in enumerate(qal):
                        score = phred[b]
                        if score <= Q:
                            low[seq[i]] += 1

                    seq, qal = None, None

    return low, total


def plot_data(low, total, Q, out):

    # Set Colors
    gridM = '#bdbdbd'
    gridm = '#d9d9d9'
    
    # prepare to plot data
    x = []
    y1, y2 = [], []

    for i in total.keys():
        x.append(i)
        y1.append(low[i])
        y2.append(total[i])

    t1, t2 = sum(y1), sum(y2)

    # Build the Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15,10), sharex=True)

    # Plot titles and labels
    ax1.set_title(f'Nucleotides with Q <= {Q}', fontsize=28)
    ax1.set_ylabel(f'Base Pair Count', fontsize=24)
    ax2.set_title(f'Total Nucleotide Count', fontsize=28)
    ax2.set_xlabel('Nucleotide', fontsize=24)
    ax2.set_ylabel(f'Base Pair Count', fontsize=24)

    # Plots
    ax1.bar(x, y1, width=0.8, align='center')
    ax2.bar(x, y2, width=0.8, align='center')

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

    # percent of total to top of bar
    for p in ax1.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy() 
        ax1.annotate(
                    f'{height/t1*100:.2f}%',
                    (x+width/2, y + height + 0.05),
                    fontsize=12,
                    horizontalalignment='center'
                    )
    for p in ax2.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy() 
        ax2.annotate(
                    f'{height/t2*100:.2f}%',
                    (x+width/2, y + height + 0.05),
                    fontsize=12,
                    horizontalalignment='center'
                    )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    plt.savefig(f'{out}')
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

    # read all files in directory and count ATCGN
    low, total = count_ATCGN(
                            fdir,
                            phred,
                            args['Phred_Quality_Threshold']
                            )

    # Plot the count data
    _ = plot_data(
                low,
                total,
                args['Phred_Quality_Threshold'],
                args['out_file']
                )

if __name__ == "__main__":
    main()