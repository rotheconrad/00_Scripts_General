#!~/.conda/envs/rothon/bin/python

## USAGE :: python $scriptname $infile
## This script is used to reduce a distance matrix such as a mash distance matrix
## to n=2 and plot the matrix. The m and c arrays are for markers and colors.
## The current implementation is to assigne color and marker to each row of the matrix
## individually. This is simply done in sequential order.

# The stress value is calculated and normalized and printed on the plot.
# The stress value can be interpretted as follows:
# According to Kruskal (1964, p. 3): value 0 indicates "perfect" fit, 0.025 excellent,
# 0.05 good, 0.1 fair, and 0.2 poor.
# For more information cf. Kruskal (1964, p. 8-9) and Borg (2005, p. 41-43).

import sys, matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.manifold import MDS

def NMDS_Plot(df):

	outFile = sys.argv[1] + '.NMDS.pdf'

	m = ['.', 'v', '^', '<', '>',
		   '.', 'v', '^', '<', '>',
		   '.', 'v', '^', '<', '>',
		   '.', 'v', '<', '>']

	c = ['#e41a1c', '#e41a1c', '#e41a1c', '#e41a1c', '#e41a1c',
		  '#377eb8', '#377eb8', '#377eb8', '#377eb8', '#377eb8',
		  '#4daf4a', '#4daf4a', '#4daf4a', '#4daf4a', '#4daf4a',
		  '#984ea3', '#984ea3', '#984ea3', '#984ea3']

	seed = np.random.RandomState(seed=3)
	nmds = MDS(n_components=2, metric=False, max_iter=10000, eps=1e-9,
		dissimilarity="precomputed", random_state=seed, n_jobs=1, n_init=10000)

	npos = nmds.fit_transform(df)
	
	distances = []

	for i in range(len(npos)-1):
		j = i+1
		distances.append(np.sqrt( (npos[j,0]-npos[i,0])**2 + (npos[j,1]-npos[i,1])**2) )
	
	stress = np.sqrt(nmds.stress_ / sum(distances)**2)
	stext = 'Stress = %f' % (stress)

	for i in range(len(npos)):
		plt.scatter(npos[i, 0], npos[i, 1], c=c[i], marker=m[i], label=df.index[i])

	plt.subplots_adjust(right=0.7)
	ax = plt.gca()
	plt.axis('equal')
	plt.title('NMDS plot of Mash Distance')
	plt.legend(frameon=False, bbox_to_anchor=(1.04,1.05), loc="upper left")
	plt.text(1, .010, stext, fontsize=10, color='#737373', horizontalalignment='right', transform=ax.transAxes)
	plt.savefig(outFile)
	plt.close()

def main():
	df = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
	NMDS_Plot(df)

if __name__ == "__main__":
	main()
