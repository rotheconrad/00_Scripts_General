#!/usr/bin/env python

''' Gets Fastq Reads matching Filtered Blast Output.

For retrieving reads that map to a reference genome or MAG.

Takes query fasta file and tabular blast output and returns
the reads matching the blast file. The blast output should be
filtered prior to running this script for best hit and any
desired pIdent or match length cutoffs.

Should work with Blast+ or MagicBlast tabular outputs.

Returns reads in fasta format. Converts fastq to fasta.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: October 29th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def read_fastq(fq):
    linecount = 0
    while fq:
        name = fq.readline().rstrip().split(' ')[0][1:]
        if not name: break
        seq = fq.readline().rstrip()
        blank = fq.readline().rstrip()
        qal = fq.readline().rstrip()
        linecount += 4
        if linecount % 4 == 0: yield (name, seq, qal)

def get_matching_reads(blast, query, outfile):
    """Reads files and writes fasta output of matching reads"""

    blastmatch = {}
    outbase = blast.split('/')[-1]
    tag = outbase.split('.')[0]
    readcount = 0

    with open(blast, 'r') as b:
        for l in b:
            qry = l.split('\t')[0]
            blastmatch[qry] = ''

    with open(query, 'r') as q, open(outfile, 'w') as o:
        for name, seq, qal in read_fastq(q):
            if name in blastmatch:
                readcount += 1
                o.write(f'>{tag}_{readcount:07}\n{seq}\n')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--tabular_blast_file',
        help='Please specify the filtered tabular blast file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-q', '--query_fasta_file',
        help='Please specify the query fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--out_file',
        help='What do you want to name the output file?',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('Running Script...')
    get_matching_reads(
                    args['tabular_blast_file'],
                    args['query_fasta_file'],
                    args['out_file']
                    )


if __name__ == "__main__":
    main()