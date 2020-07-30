#!/usr/bin/env python

'''Retrieves MiGA's ESS genes in nucleotide format.

You will need the following files from MiGA:

    - 07.annotation/01.function/01.essential/*.ess/proteins.aln
    - 06.cds/*.fna.gz # guznip this file

Give these two files to the script and it will write out the ESS genes
in as nucleotides in fasta format.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, July 15th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


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


def get_ess_nucs(aln, fna):

    ess = {}

    with open(aln, 'r') as file:
        for line in file:
            if line[0] == '#':
                X = line.split(' ')
                gene = X[3]
                prot = X[1]
                ess[gene] = f'>{gene}-{prot}\n'

    outfile = '.'.join(fna.split('.')[:-1]) + '.nucs.fa'

    with open(fna, 'r') as fa, open(outfile, 'w') as out:
        for name, seq in read_fasta(fa):
            if name[1:] in ess:
                out.write(f'{ess[name[1:]]}{seq}\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-aln', '--miga_proteins_aln',
        help='This is the proteins.aln file from MiGA!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-fna', '--miga_cds_fna',
        help='This is the CDS.fna file from MiGA!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')

    aln = args['miga_proteins_aln']
    fna = args['miga_cds_fna']

    _ = get_ess_nucs(aln, fna)


if __name__ == "__main__":
    main()
