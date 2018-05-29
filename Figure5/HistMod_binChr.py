#!/usr/bin/env python

from __future__ import division
import time
import numpy as np
import pandas as pd
import sys, getopt


def main(argv):
	start_time1 = time.time()

	try:
		opts, args = getopt.getopt(argv,"ho:c:b:",["outFile=","chrom=","binsize="])
	except getopt.GetoptError:
		print 'Usage: ... '
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'Usage: ... '
			sys.exit()
		elif opt in ("-o", "--outFile"):
			outFile = arg
		elif opt in ("-c", "--chrom"):
		    chrom = int(arg)
		elif opt in ("-b", "--binsize"):
			binsize = int(arg)


	if binsize==0:
		print 'Zero bin width, exiting'
		sys.exit(2)

	chrlen_df=pd.DataFrame(columns=['chr','len'])
	chrlen_df.loc[:,'chr'] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y']
	chrlen_df.loc[:,'len'] = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]

	chr_len = int(chrlen_df.loc[chrlen_df.loc[:,'chr']==chrom,'len'])

	n_bins = np.int(np.ceil(chr_len/binsize))

	binnedChr = pd.DataFrame(columns=['chr','start', 'end', 'strand'])
	ends = range(binsize,chr_len,binsize)
	ends.extend([chr_len])
	binnedChr.loc[:,'start'] = range(1,chr_len,binsize)
	binnedChr.loc[:,'end'] = ends
	binnedChr.loc[:,'chr'] = (len(ends))*['chr'+str(chrom)]
	binnedChr.loc[:,'strand'] = (len(ends))*['+']

	binnedChr.to_csv(outFile, sep='\t', header=False, index=False)


	elapsed_time1 = time.time() - start_time1

if __name__ == "__main__":
	import warnings
	warnings.filterwarnings("ignore")
	main(sys.argv[1:])
