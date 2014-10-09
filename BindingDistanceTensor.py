import sys,getopt
import numpy as np
from collections import defaultdict

def tf_list(fp):
	f = open(fp)
	TF = []
	for line in f:
		TF.append(line.strip())
	f.close()
	return TF

def tissue_ids(tissue_file):
	f = open(tissue_file)
	T = []
	for line in f:
		T.append(line.strip())
	f.close()
	return T

# need to support RNA-Seq only here...quantile normalize within tissues and all that
# chr start end strand gene chr start end state overlap_len
def gene_chromhmm_states(fp):
	GeneState = defaultdict(float)
	current_gene = None
	current_overlap = 0
	start = None
	end = None
	chr = None
	f = open(fp)
	for line in f:
		ls = line.strip().split('\t')
		if ls[4] != current_gene:
			if current_overlap > 0:
				GeneState[current_gene] = GeneState[current_gene]/current_overlap
			current_overlap = 0
			current_gene = ls[4]
			chr = ls[0]
			start = int(ls[1])
		# DUMB MODEL: eg Tx, Prom are low numbers, Quie and Repr are high
		state = int(ls[8][:ls[8].index('_')])
		GeneState[current_gene] += state*int(ls[9])
		current_overlap += int(ls[9])
		end = int(ls[2])
	f.close()
	GeneState[current_gene] = GeneState[current_gene]/current_overlap
	return GeneState

#$ bedtools window -a ENCODE/tss.unique.bed -b all_regions.bed -w 10000 > ENCODE/tss.unique.enhancerOverlap.bed
def gene_inflation_idx(tss_overlap_file,enhancer_labels,gene_labels):
	EtoG = defaultdict(list)
	f = open(tss_overlap_file)
	for line in f:
		ls = line.strip().split('\t')
		if ls[5] != '.':
			tss_pos = int(ls[1])
			enh_pos = (int(ls[6]) + int(ls[7]))/2
			if (ls[5],int(ls[6])) in enhancer_labels:
				if ls[3] == '+':
					genomic_dist = tss_pos - enh_pos
				else:
					genomic_dist = enh_pos - tss_pos
				EtoG[enhancer_labels[(ls[5],int(ls[6]))]].append((gene_labels[ls[4]],genomic_dist))
	f.close()
	return EtoG

def inflate_to_genes(subs,vals,EtoG,thresh=0.1):
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
	if distance > -1:
		return score/np.arcsinh((distance+1.)/500)
	else:
		return score/np.arcsinh(abs(distance)+1.)

def reduce_to_nonzero_indices(original_idx,labels):
	if isinstance(labels,dict):
		L = [None for _ in range(max(labels.values())+1)]
		for k,v in labels.items():
			L[v] = k
		L = np.array(L)
	elif isinstance(labels,list):
		L = np.array(labels)
	nonzeroIdx = list(set(original_idx))
	lookup = dict([(x,i) for i,x in enumerate(nonzeroIdx)])
	new_idx = [lookup[x] for x in original_idx]
	return new_idx,L[nonzeroIdx]

#expression_file = '/broad/compbio/zhizhuo/projects/enhancer_study/result/2014-03-07/3a/57epigenomes.RPKM.pc'
#ensembleid_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/EnsemblIDToGene.txt'
#tf_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/tflist.txt'
#tissue_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/tissue_list.txt'
#bed_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/filtered_motif_hits4/regions.peaks.txt'
#rank = [150,20,15,100]
#tss_file = '/ahg/regevdata/users/bcleary/MotifAnalysis/ENCODE/tss.unique.enhancerPromoterOverlap.bed'
#chrom_path = '/ahg/regevdata/users/bcleary/MotifAnalysis/data/ChromHMM/'
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'f:t:s:c:o:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-t'):
			tissue_file = arg
		elif opt in ('-f'):
			tf_file = arg
		elif opt in ('-s'):
			tss_file = arg
		elif opt in ('-c'):
			chrom_path = arg
		elif opt in ('-o'):
			outfile = arg
	TF = tf_list(tf_file)
	Tissues = tissue_ids(tissue_file)
	GeneTissueStates = [gene_chromhmm_states('%s/%s.tss.bed' % (chrom_path,t)) for t in Tissues]
	gene_labels = dict([(k,i) for i,k in enumerate(GeneTissueStates[0].keys())])
	subs = tuple(np.load(outfile + '.binding.subs.npy'))
	vals = np.load(outfile + '.binding.vals.npy')
	enhancer_labels = np.load(outfile + '.binding.enhancer_labels.npy').item()
	EtoG = gene_inflation_idx(tss_file,enhancer_labels,gene_labels)
	subs,vals = inflate_to_genes(subs,vals,EtoG)
	Labels = [enhancer_labels,TF,Tissues,gene_labels]
	newSubs = []
	newLabels = []
	for i in range(len(subs)):
		ns,l = reduce_to_nonzero_indices(subs[i],Labels[i])
		newSubs.append(ns)
		newLabels.append(l)
	
	np.save(outfile + '.bindingDistance.subs.npy',newSubs)
	np.save(outfile + '.bindingDistance.vals.npy',vals)
	np.save(outfile + '.bindingDistance.enhancer_labels.npy',newLabels[0])
	np.save(outfile + '.bindingDistance.tf_labels.npy',newLabels[1])
	np.save(outfile + '.bindingDistance.tissue_labels.npy',newLabels[2])
	np.save(outfile + '.bindingDistance.gene_labels.npy',newLabels[3])