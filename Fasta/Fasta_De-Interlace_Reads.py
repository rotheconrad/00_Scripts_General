#!/usr/local/pacerepov1/python/2.7/bin/python

## USAGE :: python scriptname.py read1.fasta prefix
## De-Interlaces an interlaced fasta file separating the interlaced file into
## prefix_P1.fasta and prefix_P2.fasta

import sys

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

def Fasta_De_interlace_reads(infile, prefix):

	d = {}
	P1name = prefix + '_P1.fasta'
	P2name = prefix + '_P2.fasta'

	with open(infile, 'r') as f1:
		with open(P1name, 'w') as P1out:
			with open(P2name, 'w') as P2out:
				for name, seq in read_fasta(f1):
					if name in d:
						P2out.write(name + '\n' + seq + '\n')
					else:
						P1out.write(name + '\n' + seq + '\n')
						d[name] = ''

def main():
	infile = sys.argv[1]
	prefix = sys.argv[2]
	Fasta_De_interlace_reads(infile, prefix)
	
if __name__ == "__main__":
	main()

