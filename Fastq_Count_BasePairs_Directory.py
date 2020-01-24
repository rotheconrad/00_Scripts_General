#!/usr/bin/env python

'''Calculate base pairs for all fastq files in a directory

This is useful for counting the total number of base pairs sequenced
for a given directory (project). The directory given should contain only
fastq files to be counted.

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

def evaluate_data(fdir):

    f_list = [f for f in os.listdir(fdir) if os.path.isfile(f'{fdir}/{f}')]
    
    total = 0

    for i, file in enumerate(f_list):

        print(f'{file}: {i+1:03}')
        line_count = 0

        with open(f'{fdir}/{file}', 'r') as f:
            for l in f:
                line_count += 1
                if line_count%4 == 0:
                    total += len(l.rstrip())

    print(f'Total Bases\t{total}')

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-fdir', '--file_directory',
        help='Please specify the directory of fastq files!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    fdir = args['file_directory']

    if fdir[-1] == '/': fdir = fdir[:-1]

    evaluate_data(fdir)

if __name__ == "__main__":
    main()