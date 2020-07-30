#!/usr/bin/env python

''' Collect MAGs in fasta format from GraphBin output.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 23rd, 2020d
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def get_graphbin_fasta(Assembly, Graphbin, Prefix):
    ''' Collects contig and bin names from graphbin output file and
        writes out fasta files for each MAG. Collects fastas sequence
        from the assembly fasta file'''

    MAGs = {}

    with open(Graphbin, 'r') as bins:
        for b in bins:
            X = b.rstrip().split(',')
            contig = X[0]
            mag = X[1]
            MAGs[contig] = mag

    data = defaultdict(list)

    with open(Assembly, 'r') as fasta:
        for name, seq in read_fasta(fasta):
            node = '_'.join(name.split('_')[1:3])
            fa = f'{name}\n{seq}\n'
            
            if node in MAGs:
                bin_number = MAGs[node]
                data[bin_number].append(fa)

    for MAG, sequences in data.items():
        print(f'\nWriting MAG {MAG}...\n')
        outfile = f'{Prefix}_{MAG}.fasta'
        with open(outfile, 'w') as out:
            for sequence in sequences:
                out.write(sequence)


def main():

        # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-a', '--assembly_fasta_file',
        help='Please specify the assembly file in fasta format!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-g', '--graphbin_output_file',
        help='Please specify the graphbin output file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    get_graphbin_fasta(
                    args['assembly_fasta_file'],
                    args['graphbin_output_file'],
                    args['output_file_prefix']
                    )

if __name__ == "__main__":
    main()