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
	Treatment = []
	Samples = list(df.Samples)
	Treatment = [s[:2] for s in Samples]

#	for index, row in df.iterrows():
#		Treatment.append(row['Samples'][:2])

	df['Treatment'] = Treatment
	tcount = list(df.groupby('Treatment').count().Samples)
	t = list(df.Treatment.unique())
	tt = list(zip(t, tcount))

	print(t)
	print(tcount)
	
	del df['Samples']
	del df['Treatment']

	m = [ 'o', 'v', '^', 's', 'p',
	     'P', '*', 'X', 'd', 'x' ]

	c = [ '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
	      '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a' ]

	seed = np.random.RandomState(seed=3)
	nmds = MDS(n_components=2, metric=False, max_iter=10000, eps=1e-12,
		dissimilarity="precomputed", random_state=seed, n_jobs=1, n_init=10000)

	npos = nmds.fit_transform(df)
	
	distances = []

	for i in range(len(npos)-1):
		j = i+1
		distances.append(np.sqrt( (npos[j,0]-npos[i,0])**2 + (npos[j,1]-npos[i,1])**2) )
	
	stress = np.sqrt(nmds.stress_ / sum(distances)**2)
	stext = 'Stress = %f' % (stress)

	df2 = pd.DataFrame({'Sample': Samples, 'X': npos[:,0], 'Y': npos[:,1]})
	df_file = sys.argv[1] + '_NMDS_DataFrame.tsv'
	df2.to_csv(df_file, sep='\t')

	a = 0
	for i,j in enumerate(tt):
		b = a+j[1]
		plt.scatter(npos[a:b, 0], npos[a:b, 1], c=c[i], marker=m[i], label=j[0])
		a += j[1]

	plt.subplots_adjust(right=0.7)
	ax = plt.gca()
	plt.axis('equal')
	plt.title('NMDS plot of Mash Distance')
	plt.legend(frameon=False, bbox_to_anchor=(1.04,1.05), loc="upper left")
	plt.text(1, .010, stext, fontsize=10, color='#737373', horizontalalignment='right', transform=ax.transAxes)
	plt.savefig(outFile)
	plt.close()

def main():
	df = pd.read_csv(sys.argv[1], sep='\t')
	NMDS_Plot(df)

if __name__ == "__main__":
	main()
