import sys,getopt
import numpy as np
from collections import defaultdict

#$ bedtools window -a ENCODE/tss.unique.bed -b all_regions.bed -w 10000 > ENCODE/tss.unique.enhancerOverlap.bed
def gene_inflation_idx(tss_overlap_file,enhancer_labels,gene_labels):
	EtoG = defaultdict(list)
	f = open(tss_overlap_file)
	for line in f:
		ls = line.strip().split('\t')
		if ls[5] != '.':
			tss_pos = int(ls[1])
			enh_pos = int(ls[6])
			if (ls[5],int(ls[6])) in enhancer_labels:
				if ls[3] == '+':
					genomic_dist = tss_pos - enh_pos
				else:
					genomic_dist = enh_pos - tss_pos
				EtoG[enhancer_labels[(ls[5],int(ls[6]))]].append((gene_labels[ls[4]],genomic_dist))
	f.close()
	return EtoG

def inflate_to_genes(subs,vals,EtoG,thresh=0):
	idx_enhancer = []
	idx_tissue = []
	idx_tf = []
	idx_gene = []
	scores = []
	for i,score in enumerate(vals):
		for g_idx in EtoG[subs[0][i]]:
			new_score = distance_score(score,g_idx[1])
			if new_score > thresh:
				idx_enhancer.append(subs[0][i])
				idx_tf.append(subs[1][i])
				idx_tissue.append(subs[2][i])
				idx_gene.append(g_idx[0])
				scores.append(new_score)
	return (idx_enhancer,idx_tf,idx_tissue,idx_gene),scores

# This could certainly be improved
def distance_score(score,distance):
	if distance > 0:
		return 5*score/np.arcsinh(distance/5)
	else:
		return 5*score/np.arcsinh(distance/.5)