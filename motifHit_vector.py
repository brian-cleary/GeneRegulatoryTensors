#!/usr/bin/env python

import sys, getopt
import glob, os
import gzip
import numpy as np
from collections import defaultdict
from operator import itemgetter
from scipy.spatial import distance

def ensemble_dict(fp):
	f = open(fp)
	header = f.readline()
	EnsembleIds = {}
	for line in f:
		ls = line.strip().split()
		EnsembleIds[ls[0]] = ls[1]
	f.close()
	return EnsembleIds

def expression_vector(fp,tissueid,EnsembleIds,TFlist):
	f = open(fp)
	header = f.readline().strip().split()
	colidx = header.index(tissueid)
	E = defaultdict(float)
	for line in f:
		ls = line.strip().split()
		gene = EnsembleIds.get(ls[0],None)
		expression = float(ls[colidx])
		if (gene in TFlist) and (expression > 0):
			E[gene] = expression
	f.close()
	E = sorted(E.iteritems(),key=itemgetter(1))
	return dict([(x[0],100.*i/len(E)) for i,x in enumerate(E)])

def expression_vector_newTissue(fp,EnsembleIds,TFlist):
	f = open(fp)
	E = defaultdict(float)
	for line in f:
		ls = line.strip().split('\t')
		# e.g.: ['chr1', 'gencodeV7', 'gene', '157043049', '157044494', '0.001287', '+', '.', 'gene_id "ENSG00000224520.1"; transcript_ids "ENST00000414328.1,"; RPKM1 "0.089518"; RPKM2 "0.029692"; iIDR "0.141";']
		info = ls[8].split(';')
		if 'ENSG' in info[0]:
			gene = EnsembleIds.get(info[0][info[0].index('"')+1:info[0].index('.')],None)
			expression = float(info[2].replace('"','').split()[1])
			if (gene in TFlist) and (expression > 0):
				E[gene] = expression
	f.close()
	E = sorted(E.iteritems(),key=itemgetter(1))
	return dict([(x[0],100.*i/len(E)) for i,x in enumerate(E)])

def tf_list(fp):
	f = open(fp)
	TF = []
	for line in f:
		TF.append(line.strip())
	f.close()
	return TF

def write_motif_matches(fp,op,T,chr,Expression,PeakFiles,corr_thresh=.2,score_thresh=.05,peak_thresh=.75,window_adj=8):
	f = gzip.open(fp)
	g = open(op,'w')
	g.write(' '.join(['motif','chr','start','end','DNase correlation'] + T)+'\n')
	DScoresBlock = []
	idx = -1
	for line in f:
		ls = line.strip().split()
		full_motif = ls[0][:ls[0].index('_8mer')]
		motif = full_motif[:full_motif.index('_')]
		score = float(ls[6])
		bp_score = [binding_potential(score,e.get(motif,0.)) for e in Expression]
		if max(bp_score) > score_thresh:
			start = int(ls[2]) - window_adj
			end = int(ls[3]) + window_adj
			if end > idx+len(DScoresBlock):
				DScoresBlock = DScoresBlock[start-idx:]
				for i in range(len(DScoresBlock),end-start-len(DScoresBlock)):
					DScoresBlock.append(read_dnase_peak_lines(start+i,PeakFiles))
				idx = start
			Dmax = np.max(DScoresBlock[start-idx:end-idx],0)
			#DUMBNESS WARNING: peak_thresh SET IN DUMB WAY DUMB
			if max(Dmax) > peak_thresh:
				c = 1 - distance.correlation(bp_score,Dmax)
				if c > corr_thresh:
					combined_score = combine_scores(Dmax,bp_score,c,corr_thresh)
					g.write(' '.join([full_motif,chr,str(start),str(end),str(c)] + ['%f'%bp for bp in combined_score])+'\n')
	f.close()
	g.close()

def write_motif_matches_newTissue(fp,op,T,chr,Expression,PeakFiles,newTissuePeaks,newTissueExpression,corr_thresh=.2,score_thresh=.05,peak_thresh=.75,window_adj=8):
	PeakFiles.append(newTissuePeaks)
	f = gzip.open(fp)
	g = open(op,'w')
	g.write(' '.join(['motif','chr','start','end','DNase correlation','newTissue'])+'\n')
	DScoresBlock = []
	idx = -1
	for line in f:
		ls = line.strip().split()
		full_motif = ls[0][:ls[0].index('_8mer')]
		motif = full_motif[:full_motif.index('_')]
		score = float(ls[6])
		# WOULD JUST BE GREAT TO NOT RECALCULATE THE CORRELATIONS EVERY TIME...
		bp_score = [binding_potential(score,e.get(motif,0.)) for e in Expression]
		bp_score_newTissue = binding_potential(score,newTissueExpression.get(motif,0.))
		if (max(bp_score) > score_thresh) and (bp_score_newTissue > score_thresh*5):
			start = int(ls[2]) - window_adj
			end = int(ls[3]) + window_adj
			if end > idx+len(DScoresBlock):
				DScoresBlock = DScoresBlock[start-idx:]
				for i in range(len(DScoresBlock),end-start-len(DScoresBlock)):
					DScoresBlock.append(read_dnase_peak_lines(start+i,PeakFiles))
				idx = start
			Dmax = np.max(DScoresBlock[start-idx:end-idx],0)
			#DUMBNESS WARNING: peak_thresh SET IN DUMB WAY DUMB
			if (max(Dmax[:-1]) > peak_thresh) and (Dmax[-1] > peak_thresh):
				c = 1 - distance.correlation(bp_score,Dmax[:-1])
				if c > corr_thresh:
					combined_score = combine_scores(Dmax[-1],bp_score_newTissue,c,corr_thresh)
					g.write(' '.join([full_motif,chr,str(start),str(end),str(c),str(combined_score)])+'\n')
	f.close()
	g.close()

def combine_scores(D,bp_score,c,kd_d=2,kd_c=.5):
	return D**4/(kd_d**4 + D**4)*np.array(bp_score)*(1 + c/(kd_c + c))

def binding_potential(motif_score,expression,q=1,kd=50.):
	if (expression > 0) and (motif_score > 0):
		return expression**q/(kd**q/motif_score + expression**q)
	else:
		return 0.

def tissue_ids(tissue_file):
	f = open(tissue_file)
	T = []
	for line in f:
		T.append(line.strip())
	f.close()
	return T

def read_dnase_peak_lines(motif_position,PeakFiles):
	DScores = [(-1,0) for _ in PeakFiles]
	eof = False
	for i,pf in enumerate(PeakFiles):
		while DScores[i][0] < motif_position:
			ls = pf.readline().strip().split()
			if len(ls) > 0:
				DScores[i] = (int(ls[0]),float(ls[1]))
			else:
				eof = True
				break
		if eof:
			break
	if eof:
		return None
	else:
		return [ds[1] for ds in DScores]

def precision_recall(fp,thresh):
	O = []
	R = []
	B = {}
	f = open(fp)
	for line in f:
		ls = line.strip().split('\t')
		if ls[5] == '.':
			O.append(0)
			R.append(0)
		elif (ls[5],ls[6],ls[7]) in B:
			O.append(1)
			R.append(0)
		else:
			O.append(1)
			R.append(1)
			B[(ls[5],ls[6],ls[7])] = True
		if float(ls[3]) < thresh:
			break
	return np.cumsum(O,dtype=np.float64),np.cumsum(R,dtype=np.float64)

def plot_pr(tf,n,thresh):
	p,r = precision_recall('%s.txt'%tf,thresh)
	pylab.plot(r/n,p/(np.arange(len(p))+1),'.')
	pylab.xlim((0,1))
	pylab.ylim((0,1))
	pylab.xlabel('Recall')
	pylab.ylabel('Precision')
	pylab.title('%s: n=%d, thresh=%.2f' % (tf,len(p),thresh))
	pylab.show()


#expression_file = '/broad/compbio/zhizhuo/projects/enhancer_study/result/2014-03-07/3a/57epigenomes.RPKM.pc'
#ensembleid_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/EnsemblIDToGene.txt'
#tf_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/tflist.txt'
#tissue_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/tissue_list.txt'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'m:e:i:o:f:t:d:c:1:2:',["motifhits=","expressionmatrix=","ensemblid="])
	except:
		print help_message
		sys.exit(2)
	newExpressionFile = None
	newPeakFile = None
	for opt, arg in opts:
		if opt in ('-m','--motifhits'):
			motifhit_prefix = arg
		elif opt in ('-e','--expressionmatrix'):
			expression_file = arg
		elif opt in ('-i','--ensemblid'):
			ensembleid_file = arg
		elif opt in ('-o'):
			outdir = arg
		elif opt in ('-f'):
			tf_file = arg
		elif opt in ('-t'):
			tissue_file = arg
		elif opt in ('-d'):
			dnase_path = arg
		elif opt in ('-c'):
			# FILES ARE SPLIT INTO CHROMOSOMES, CHUNKS OF 5MB
			chr,chr_pos,chr_idx = arg.split(',')
		elif opt in ('-1'):
			newExpressionFile = arg
		elif opt in ('-2'):
			newPeakFile = arg
	motifhit_file = '%s.%s.%s.matches.txt.gz' % (motifhit_prefix,chr,chr_pos)
	outfile = '%s/%s.%s.txt' % (outdir,chr,chr_pos)
	EnsembleIds = ensemble_dict(ensembleid_file)
	TFlist = tf_list(tf_file)
	T = tissue_ids(tissue_file)
	Expression = [expression_vector(expression_file,tissueid,EnsembleIds,TFlist) for tissueid in T]
	PeakFiles = [open('%s/%s.peaks.%s.txt.%s' % (dnase_path,tissueid,chr,chr_idx)) for tissueid in T]
	if (newExpressionFile != None) and (newPeakFile != None):
		newExpression = expression_vector_newTissue(newExpressionFile,EnsembleIds,TFlist)
		write_motif_matches_newTissue(motifhit_file,outfile,T,chr,Expression,PeakFiles,open('%s.%s.txt.%s' % (newPeakFile,chr,chr_idx)),newExpression)
	else:
		write_motif_matches(motifhit_file,outfile,T,chr,Expression,PeakFiles)
	
	