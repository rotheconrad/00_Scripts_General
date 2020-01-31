#!/usr/bin/env python

''' Randomly Sample Genomic Fasta Files Above ANI Threshold.

This script takes a directory of genomic fasta files and selects a
seed genome at random. It then continues to randomly select genomes,
comparing ANI with the seed genome using the fastANI program. If the
ANI is within the user defined ANI range, the genome is kept. If the ANI
is outside of the ANI range, the next genome is selected at random until
the user defined number of genomes (n) is obtained or all of the genomes
in the provided directory have been sampled.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This script requires the fastANI program to be installed and  !!!!
!!!! availble through the system path. github.com/ParBLiSS/FastANI !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 20th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys, os, re, random, subprocess

def choice_continue():
    yes = {'yes','y', 'ye', 'yep', 'yeah', 'sure'}
    no = {'no','n', 'nope', 'not'}
    choice = input("Would you like to try again? [y/n]").lower().strip()
    if choice in yes:
        return True
    elif choice in no:
        return False
    else:
        sys.stdout.write("Please respond with (yes / y) or (no / n)\n\n")
        return choice_continue()

def if_gzip(new_dir, name):
    _ = subprocess.run(['gunzip', f'{new_dir}{name}'])
    name = '.'.join(name.split('.')[:-1])
    return name

def setup_the_bomb(t, n, gd, exp):
    ''' Initialize new folder and dictionary. Sets up the bomb '''

    # initialize a dictionary for {NewFileName: OldFileName}
    file_name_key = {} 
    # get the list of files from the genome directory
    if gd[-1] != '/': gdir = gdir + '/'
    g_list = [g for g in os.listdir(gd) if os.path.isfile(f'{gd}{g}')]
    # if directory exists, exit. else make a new directory for this trial.
    new_dir = f'{exp}_RndmGnms_{t:04}/'
    if os.path.isdir(new_dir):
        print(f'{new_dir} already exists!!')
        print('Please choose a directory name that does not exist\n\n')
        sys.exit()
    else:
        os.makedirs(new_dir)

    print('Selecting Random Seed Genome 1 to use as ANI reference.')
    print('#######################################################\n')
    # Select the first genome of this sample, copy, unzip, and rename it
    g_prime = random.choice(g_list) 
    _ = g_list.remove(g_prime)
    _ = subprocess.run(['cp', f'{gd}{g_prime}', new_dir])
    if g_prime.split('.')[-1] == 'gz': g_prime = if_gzip(new_dir, g_prime)
    new_name = f'{exp}_{t:04}_RndmGnm_0001.fna'
    _ = subprocess.run(['mv', f'{new_dir}{g_prime}', f'{new_dir}{new_name}'])
    # Update the file name dictionary
    file_name_key[new_name] = g_prime

    return file_name_key, g_list, new_dir, new_name

def select_genomes(t, n, gd, exp, atl, atu):
    '''
    Selects n genomes randomly from genome directory that are within the
    ANI threshold. Renames them iteratively and copies them to a new
    directory (n_organism_name_samples). 

    Also writes a key to match n sample to original genome file name.
    '''
    # Initialize random trial
    file_name_key, g_list, new_dir, g_prime = setup_the_bomb(t, n, gd, exp)
    genomes_avail = len(g_list)+1
    fastANIout = open(f'{new_dir}fastANIout.log', 'a')

    i, j = 1, 1
    while (j <= n):
        # if all genomes in directory have been tested against the seed genome
        # and the requested number of genome have not been found stop.
        # ask to try again or exit.
        if len(g_list) == 0:
            print('\n############# FAIL ####################################')
            print(
                f'Sampled all {genomes_avail} genomes in {gd} directory and '
                f'only found {j} of the {n} matches requested between {atl}% '
                f'and {atu}% ANI!\n\n'
                )
            # Check if the user wants to try again or quit
            if choice_continue():
                _ = subprocess.run(['rm', '-r', new_dir])
                select_genomes(t, n, gd, exp, atl, atu)
            else: sys.exit()

        # randomly select genome from the file list.
        i += 1
        print(f'\nTesting Random Genome {i} of {genomes_avail}')
        g_new = random.choice(g_list)
        # Remove genome from the file list before the next iteration.
        # ie random sample without replacement.
        print('Number of genomes remaining:', len(g_list))
        _ = g_list.remove(g_new)
        # copy random genome to new directory with new naming scheme.
        _ = subprocess.run(['cp', f'{gd}{g_new}', new_dir])
        # if the genome fasta is gzipped, unzip.
        if g_new.split('.')[-1] == 'gz': g_new = if_gzip(new_dir, g_new)
        # rename new random genome to fit scheme.
        new_name = f'{exp}_{t:04}_RndmGnm_{j+1:04d}.fna'
        _ = subprocess.run(['mv', f'{new_dir}{g_new}', f'{new_dir}{new_name}'])
        # run fastANI
        fastANI = (
            f"fastANI -r {new_dir}{g_prime} -q {new_dir}{new_name}"
            f" -o {new_dir}temp.ani"
            )
        _ = subprocess.run(
                            fastANI,
                            shell=True,
                            stderr=subprocess.STDOUT,
                            stdout=fastANIout
                            )
        # if fastANI successful, compare ANIs, print results, and continue.
        if os.stat(f'{new_dir}temp.ani').st_size > 0:
            with open(f'{new_dir}temp.ani', 'r') as f:
                ani = float(f.readline().split('\t')[2])

            if ani < atl:
                _ = subprocess.run(['rm', f'{new_dir}{new_name}'])
                print(f'FAIL: {ani}% ANI is below the {atl}% ANI threshold!')

            elif ani > atu:
                _ = subprocess.run(['rm', f'{new_dir}{new_name}'])
                print(f'FAIL: {ani}% ANI is above the {atu}% ANI threshold!')

            elif ani >= atl and ani <= atu:
                file_name_key[new_name] = g_new
                j += 1
                print(
                    f'PASS: {ani}% ANI is within the {atl}% to {atu}% threshold'
                    f' range!\nFound {j} genomes of the {n} requested.'
                    )
            # if ANI fails threshold tests
            else:
                print(ani, atl, atu)
                print('FREAK OUT!! Error during threshold test lines 129-151.')

        # once finished selecting genomes remove the temp.ani file
        _ = subprocess.run(['rm', f'{new_dir}temp.ani'])
    # close the fastANIout.log file
    fastANIout.close()
    # Write out genome name key
    with open(f'{new_dir}00_Genome_Name_Key.tsv', 'w') as o:
        for k,v in file_name_key.items():
            o.write(f'{k}\t{v}\n')
    # Print closing text
    print('\n\n############# SUCCESS ####################################')
    print(
        f'{j-1} of {n} Genomes between {atl}% and {atu}% ANI '
        f'sampled successfully!'
        )
    print(f'Genome samples copied, Genome name key written to {new_dir}')
    print(f'Congratulations and enjoy your research.')


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-t', '--trial_number',
        help='Trial number (to keep track of multiple random trials).',
        metavar=':',
        type=int,
        required=True
        )
    parser.add_argument(
        '-n', '--n_genomes',
        help='How many genomes would you like to select?',
        metavar=':',
        type=int,
        required=True
        )
    parser.add_argument(
        '-g', '--genomes_directory',
        help='Specify directory containing genomes. (fasta or gzipped fasta)',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-e', '--experiment_name',
        help='What is the name of the experiment?',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--ani_threshold_lower',
        help='Specify the Lower ANI Threshold <= Genome Selection. (eg 95)',
        metavar=':',
        type=float,
        required=True
        )
    parser.add_argument(
        '-u', '--ani_threshold_upper',
        help='Specify the Upper ANI Threshold >= Genome Selection. (eg 100)',
        metavar=':',
        type=float,
        required=True
        )
    args=vars(parser.parse_args())
    # Run this scripts main function
    print('\n\nRunning Script...\n\n')
    select_genomes(
        args['trial_number'],
        args['n_genomes'],
        args['genomes_directory'],
        args['experiment_name'],
        args['ani_threshold_lower'],
        args['ani_threshold_upper']
        )

if __name__ == "__main__":
    main()
