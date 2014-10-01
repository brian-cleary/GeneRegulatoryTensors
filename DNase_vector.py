#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np

# $ samtools view /broad/compbio/zhizhuo/projects/enhancer_study/data/2014-04-02/DNaseBAM/E003-DNase.bam -o - | python DNase_vector.py [options]
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'o:c:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-o'):
			outfile = arg
		elif opt in ('-c'):
			chr = arg
	D = np.zeros((4*10**8,))
	max_idx = 0
	for line in sys.stdin:
		ls = line.strip().split()
		if ls[2] == chr:
			idx = int(ls[3])
			D[idx] += 1
			if idx > max_idx:
				max_idx = idx
	D = D[:max_idx]
	np.save(outfile,D/D.sum()*max_idx)
