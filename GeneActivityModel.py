import sys,getopt
import numpy as np
import sktensor as sk
from collections import defaultdict
import itertools
from scipy.sparse.linalg import svds

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

def expand_tfs_polynomial(subs,vals,shape,tf_labels,enh_mode=0,tf_mode=1,tissue_mode=2,gene_mode=3,order=3):
	newSubs = subs
	newVals = vals
	BindingAtGene = {}
	for i in xrange(len(vals)):
			g_idx = subs[gene_mode][i]
			if g_idx not in BindingAtGene:
				BindingAtGene[g_idx] = defaultdict(list)
			BindingAtGene[(g_idx,subs[enh_mode][i],subs[tissue_mode][i])].append((subs[tf_mode][i],vals[i]))
	newTFs = []
	newTF_dict = {}
	num_single_tfs = len(tf_labels)
	for idx,binding_list in BindingAtGene.iteritems():
		for o in range(2,order+1):
			for tf_combo in itertools.combinations(binding_list,o):
				tfs = tuple([t[0] for t in tf_combo])
				scores = [t[1] for t in tf_combo]
				if tfs not in newTF_dict:
					newTFs.append(tfs)
					newTF_dict[tfs] = len(newTF_dict)
				newSubs[enh_mode].append(idx[1])
				newSubs[tf_mode].append(num_single_tfs + newTF_dict[tfs])
				newSubs[tissue_mode].append(idx[2])
				newSubs[gene_mode].append(idx[0])
				newVals.append(np.product(scores))
	return newSubs,newVals,tf_labels+newTFs

def regression_svd(subs,vals,shape,enh_mode=0,tf_mode=1,kmax=10000):
	S = sk.sptensor(subs,vals,shape=shape)
	# constant enhancer weights...not so good maybe
	We = np.ones((1,S.shape[enh_mode]))
	S = S.ttv(We,enh_mode).unfold(tf_mode).tocsr()
	# might need transpose here
	u,s,vt = svds(S,k=min(min(S.shape),kmax))
	return S,u,s

def truncated_inverse(u,s,k):
	s = np.diag(1/s[:k])
	return np.dot(np.dot(u[:,:k],s),u.T[:k])

def truncated_regression_solution(S,Zinv,Y):
	return np.dot(np.dot(Zinv,S.T),Y)
