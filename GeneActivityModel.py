import sys,getopt
sys.path.insert(0,'/ahg/regev/users/bcleary/src/lib/python2.7/site-packages/')
import numpy as np
import sktensor as sk
from collections import defaultdict
import itertools
from scipy.sparse.linalg import svds
from scipy.sparse import csr_matrix

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

def expand_tfs_polynomial(subs,vals,shape,tf_labels,enh_mode=0,tf_mode=1,tissue_mode=2,gene_mode=3,order=2,thresh_pctile=10):
	newSubs = [list(s) for s in subs]
	newVals = list(vals)
	BindingAtGene = defaultdict(list)
	tf_labels = list(tf_labels)
	for i in xrange(len(vals)):
			g_idx = subs[gene_mode][i]
			binding_key = (g_idx,subs[enh_mode][i],subs[tissue_mode][i])
			BindingAtGene[binding_key].append((tf_labels[subs[tf_mode][i]],vals[i]))
	newTFs = []
	newTF_dict = {}
	num_single_tfs = len(tf_labels)
	thresh = np.percentile(vals,thresh_pctile)**order
	for idx,binding_list in BindingAtGene.iteritems():
		for o in range(2,order+1):
			for tf_combo in itertools.combinations_with_replacement(binding_list,o):
				tfs = tuple(sorted([t[0] for t in tf_combo]))
				scores = [t[1] for t in tf_combo]
				newScore = np.product(scores)
				if newScore > thresh:
					if tfs not in newTF_dict:
						newTFs.append(tfs)
						newTF_dict[tfs] = len(newTF_dict)
					newSubs[enh_mode].append(idx[1])
					newSubs[tf_mode].append(num_single_tfs + newTF_dict[tfs])
					newSubs[tissue_mode].append(idx[2])
					newSubs[gene_mode].append(idx[0])
					newVals.append(newScore)
	newTFs = [','.join(tfs) for tfs in newTFs]
	return newSubs,newVals,tf_labels+newTFs

# Should consider a model with nonnegative weights, inverse binding scores, so forth...
def regression_svd(subs,vals,shape,enh_mode=0,tf_mode=1,kmax=10000):
	S = sk.sptensor(subs,vals,shape=shape)
	# constant enhancer weights...not so good maybe
	We = np.ones(S.shape[enh_mode])
	if enh_mode < tf_mode:
		tf_mode -= 1
	S = S.ttv(We,modes=[enh_mode]).unfold(tf_mode).tocsr()
	u,s,vt = svds(S,k=min(min(S.shape),kmax)-1)
	return S,u,s

def truncated_inverse(u,s,k):
	s = np.diag(1/s[-k:])
	return np.dot(np.dot(u[:,-k:],s),u.T[-k:]).T

def truncated_regression_solution(S,Zinv,GeneTissueState,gene_labels):
	gene_labels = list(gene_labels)
	geneDict = dict([(g,i) for i,g in enumerate(gene_labels)])
	Y = np.zeros((1,len(GeneTissueState),len(gene_labels)))
	for i,GT in enumerate(GeneTissueState):
		for gene,state in GT.iteritems():
			if gene in geneDict:
				Y[0,i,geneDict[gene]] = state
	Y = sk.dtensor(Y).unfold(0).transpose()
	# seems to lose sparsity during dot here...
	SZ = S.T.dot(Zinv.T).T
	return np.dot(SZ,Y)

def eliminate_elements(subs,vals,labels,mode,tissue_mode=2,num_tissues=None,minMaxPercentile=0.5):
	if num_tissues == None:
		num_tissues = max(subs[tissue_mode])+1
	num_mode = len(labels)
	TissueElementSums = np.zeros((num_tissues,num_mode))
	for i,v in enumerate(vals):
		TissueElementSums[subs[tissue_mode][i],subs[mode][i]] += v
	TissueElementPercentiles = np.zeros(TissueElementSums.shape)
	for i,t in enumerate(TissueElementSums):
		nnz = np.nonzero(t)[0]
		xs = np.argsort(t[nnz])
		for rank,idx in enumerate(xs):
			TissueElementPercentiles[i,nnz[idx]] = rank/float(len(nnz))
	MaxPercentiles = TissueElementPercentiles.max(0)
	return MaxPercentiles
	newSubs = [[] for _ in range(len(subs))]
	newVals = []
	newLabel_dict = {}
	for i,v in enumerate(vals):
		s_idx = subs[mode][i]
		if MaxPercentiles[s_idx] > minMaxPercentile:
			for j,s in enumerate(subs):
				if j != mode:
					newSubs[j].append(s[i])
			newVals.append(v)
			if s_idx not in newLabel_dict:
				newLabel_dict[s_idx] = len(newLabel_dict)
			newSubs[mode].append(newLabel_dict[s_idx])
	newLabels = sorted([(v,labels[k]) for k,v in newLabel_dict.iteritems()])
	newLabels = [l[1] for l in newLabels]
	return newSubs,newVals,newLabels

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
		elif opt in ('-c'):
			chrom_path = arg
		elif opt in ('-o'):
			outfile = arg
	Tissues = tissue_ids(tissue_file)
	tf_labels = np.load(outfile + '.bindingDistance.tf_labels.npy')
	subs = np.load(outfile + '.bindingDistance.subs.npy')
	vals = np.load(outfile + '.bindingDistance.vals.npy')
	vals[(vals == np.inf)] = vals[(vals < np.inf)].max()
	shape = tuple([max(s)+1 for s in subs])
	newSubs,newVals,newTF_labels = expand_tfs_polynomial(subs,vals,shape,tf_labels)
	np.save(outfile + '.bindingDistance.TFexpanded.subs.npy',newSubs)
	np.save(outfile + '.bindingDistance.TFexpanded.vals.npy',newVals)
	np.save(outfile + '.bindingDistance.TFexpanded.tf_labels.npy',newTF_labels)
	shape = tuple([max(s)+1 for s in newSubs])
	newSubs = tuple(newSubs)
	S,u,s = regression_svd(newSubs,newVals,shape)
	np.save(outfile + '.regressionModel.svd.u.npy',u)
	np.save(outfile + '.regressionModel.svd.s.npy',s)
	GeneTissueStates = [gene_chromhmm_states('%s/%s.tss.bed' % (chrom_path,t)) for t in Tissues]
	gene_labels = np.load(outfile + '.bindingDistance.gene_labels.npy')
	for k in [10,100,1000,5000,10000]:
		StSinv = truncated_inverse(u,s,k)
		Wtf = truncated_regression_solution(S,StSinv,GeneTissueStates,gene_labels)
		np.save(outfile + '.regressionModel.Wtf.%d.npy' % k,Wtf)
