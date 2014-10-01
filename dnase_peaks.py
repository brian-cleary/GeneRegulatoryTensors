#!/usr/bin/env python

import sys, getopt
import glob, os
import numpy as np
import pywt


if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'o:d:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-o'):
			outfile = arg
		elif opt in ('-d'):
			dnase_file = arg
	D = np.load(dnase_file)
	s = D.shape
	# approximation at 4nt
	ca4 = pywt.downcoef('a',D,'sym4',level=2)
	ca4 = ca4/np.percentile(ca4,99.9)
	f = open(outfile,'w')
	for i in xrange(s[0]):
		x = ca4[int(i/4.)]
		f.write('%d %f\n' % (i,x))
	f.close()
