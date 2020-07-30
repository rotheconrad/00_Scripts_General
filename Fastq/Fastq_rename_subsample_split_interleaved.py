#!/usr/bin/env python

''' Rename reads, subsample, and split interleaved fastq file.

Blast and Magic Blast cut the read name at the first space character and
report the name before the space as the query ID.

However, Fastq files can be named as:
    @D00468:261:HYTMHBCX2:1:1101:9119:31637 1:N:0:CAGAGAGG+ACTGCATA
    @D00468:261:HYTMHBCX2:1:1101:9119:31637 2:N:0:CAGAGAGG+ACTGCATA
where the unique identifier can come after the space chacter.

Filtering for best hits and retrieving the fasta sequence of a magic
blast match becomes impossible once the unique identifier is lost.

This script renames the fastq sequences with a simple name:
@prefix_i

where prefix is supplied by the user and i = 1 through number of reads.

Sometimes fastq files come in interposed (interleaved) format and they
may be wanted as paired read files for trimming and assembly.

Sometimes reads sets may need to be subsampled for various reason.

Sometimes all three of these steps are need simultaneously.

This script does all three steps at once.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 6th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random
from collections import defaultdict


def read_fastq(fq):
    linecount = 0
    while fq:
        name = fq.readline().rstrip().split(' ')[0][1:]
        if not name: break
        seq = fq.readline().rstrip()
        blank = fq.readline().rstrip()
        qal = fq.readline().rstrip()
        linecount += 4
        if linecount % 4 == 0: yield (name, seq, blank, qal)


def interleaved_toDict(infile, prefix):
    '''
    Parses an interleaved (interposed) fastq file into a dictionary
    with readpair as the key and an f-string containing all 8 lines
    of both readsets
    '''

    data = defaultdict(list)
    rCount = 1 # Read Count
    pCount = 1 # Read Pair Count
    rPair = f'readPair_{pCount}' # name for read pair

    with open(infile, 'r') as f:

        for name, seq, blank, qal in read_fastq(f):

            if rCount % 2 == 0:

                entry = f'@{prefix}_{pCount:010}/2\n{seq}\n{blank}\n{qal}\n'
                data[rPair].append(entry)

                pCount += 1
                rPair = f'readPair_{pCount}'

                if pCount % 1000000 == 0:
                    print(f'\t... Read pairs processed: {pCount}')

            else:

                entry = f'@{prefix}_{pCount:010}/1\n{seq}\n{blank}\n{qal}\n'
                data[rPair].append(entry)

            rCount += 1

    totalReadPairs = len(data)

    print(f'\n\tNumber of read pairs processed: {totalReadPairs}\n\n')

    return data, totalReadPairs


def subsample_fastq(data, totalReadPairs, subsample, prefix):

    # for percent submitted, select random subsample and write to file
    for percent in subsample:
        print(f'\n\tSubsampling at {percent}% ...')
        subsample_percent = int(totalReadPairs * percent / 100)
        random.seed()
        subsampleDict = dict(random.sample(data.items(), subsample_percent))

        outFile1 = f'{prefix}_sbsmpl-{percent}_p1.fastq'
        outFile2 = f'{prefix}_sbsmpl-{percent}_p2.fastq'

        print(
            f'\tNumber of read pairs selected for {percent}% subsample: '
            f'{len(subsampleDict)}\n\tWritting reads to files: \n'
            f'\t\t{outFile1}\n\t\t{outFile2}\n'
            )

        with open(outFile1, 'w') as fout1, open(outFile2, 'w') as fout2:
            for readPair, entries in subsampleDict.items():
                fout1.write(entries[0])
                fout2.write(entries[1])

    return True


def main():
        # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-s', '--percents_to_subsample',
        help=
            'List of whole integer values as percents to subsample '
            '(ie: -s 20 40 60 80',
        metavar='',
        type=int,
        nargs='+',
        required=True
        )
    parser.add_argument(
        '-i', '--fastq_input_file',
        help='Please specify the fastq file to read!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-p', '--naming_prefix',
        help='Please specify the prefix to use for renaming sequences!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='Please specify the prefix to use for the output files!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    infile = args['fastq_input_file']
    prefix = args['naming_prefix']

    print('\n\nRunning Script ...')
    print('\n\nReading paired read interleaved fasta file ...\n')
    data, totalReadPairs = interleaved_toDict(infile, prefix)
    print('\n\nSubsampling data and writing output files ...')
    _ = subsample_fastq(
                    data,
                    totalReadPairs,
                    args['percents_to_subsample'],
                    args['output_prefix']
                    )


if __name__ == "__main__":
    main()
