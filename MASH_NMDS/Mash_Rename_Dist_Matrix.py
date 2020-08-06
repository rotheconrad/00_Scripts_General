#!/usr/local/pacerepov1/python/2.7/bin/python

## USAGE :: python $scriptname $infile
## This script is adapted to rename duplicate row and column names from a data table.
## It reads through the rows and columns appending an "a", "b", "c" etc as needed.


import sys

def rename_table(infile):

	outfile = infile.split('.')[0] + '_rnmd' + '.' + '.'.join(infile.split('.')[1:])

	with open(infile, 'r') as f, open(outfile, 'w') as o:
		colNames = f.readline().rstrip().split('\t')[1:]
		newColNames = ['Samples',]

		for n in colNames:
			newName = n.split('.')[0] # sometimes this
#			newName = n.split('.')[0].split('_')[1] # othertimes this
			newColNames.append(newName)

		o.write('\t'.join(newColNames) + '\n')		

		for l in f:
			line = l.split('\t')
			newName = line[0].split('.')[0] # sometimes this
#			newName = line[0].split('.')[0].split('_')[1] # othertimes this
			line[0] = newName
			o.write('\t'.join(line))

	print("I think we're finished here.")

def main():
	Data_Table = sys.argv[1]
	rename_table(Data_Table)

if __name__ == "__main__":
	main()
