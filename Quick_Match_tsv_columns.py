# python

## for two tsv files
## Select lines from file2 that a match a column in file1

import sys

if len(sys.argv[1:]) != 5:
    print(
        '\n\nFor two *.tsv files with 1st line header:\n'
        'This script selects lines in file 2 that match a column in file 1\n'
        'and writes them to a new file.\n\n'
        'This script takes 5 input parameters in order:\n'
        '\tp1: name of file 1 to select column to find in file 2\n'
        '\tp2: column number to select unique id from (starting at 0)\n'
        '\tp3: name of file 2 to select lines matching column in file 1\n'
        '\tp4: ncolumn number to select unique id for file 2\n'
        '\tp5: name of the output file\n'
        '\n\nExample:\npython scriptname.py p1 p2 p3 p4 p5\n\n'
        )
    sys.exit('Please provide all 5 input parameters\n\n')

# first file name
file1 = sys.argv[1]
# column to select from first file (count columns starting at 0)
colnum1 = int(sys.argv[2])
# second file name
file2 = sys.argv[3]
# column to select from second file (count columns starting at 0)
colnum2 = int(sys.argv[4])
# Output file name
outfile = sys.argv[5]

# dict of unique column ids from file 1 to look for in file 2
match = {}

with open(file1, 'r') as f1:
    header = f1.readline()
    for line in f1:
        match_id = line.rstrip().split('\t')[colnum1]
        match[match_id] = ''

with open(file2, 'r') as f2, open(outfile, 'w') as o:
    header = f2.readline()
    o.write(header)
    for line in f2:
        X = line.rstrip().split('\t')[colnum2]
        if X in match:
            o.write(line)
