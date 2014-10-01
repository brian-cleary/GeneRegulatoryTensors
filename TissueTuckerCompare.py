import sys,getopt
import numpy as np
import sktensor as sk
from collections import defaultdict
from scipy.spatial import distance

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

# format: chr start end blank blank blank chr start end motif
def parse_bed(fp,TF):
	f = open(fp)
	idx_enhancer = []
	idx_tissue = []
	idx_tf = []
	scores = []
	last_enhancer = None
	e_idx = -1
	last_tf_position = {}
	enhancer_tf_scores = defaultdict(float)
	for line in f:
		ls = line.strip().split('\t')
		enhancer = (ls[0],int(ls[1]))
		if enhancer != last_enhancer:
			for k,v in enhancer_tf_scores.items():
				idx_enhancer.append(e_idx)
				idx_tf.append(k[0])
				scores.append(v)
			enhancer_tf_scores = defaultdict(float)
			last_enhancer = enhancer
			e_idx += 1
		if ls[6] != '.':
			# TODO: ALLOW FILE TO SPECIFY SCORES
			score = 1.
			motif_start = int(ls[7])
			motif_end = int(ls[8])
			for tf_motif in ls[9].replace('"','').split(','):
				if '_' in tf_motif:
					tf = tf_motif[:tf_motif.index('_')]
				else:
					tf = tf_motif
				if tf in TF:
					if (tf,) in last_tf_position:
						last_pos,last_score = last_tf_position[(tf,)]
						if (len(set(range(motif_start,motif_end)) & last_pos) > 0):
							if score > last_score:
								enhancer_tf_scores[(TF.index(tf),)] += score - last_score
								last_tf_position[(tf,)] = (set(range(motif_start,motif_end)),score)
						else:
							enhancer_tf_scores[(TF.index(tf),)] += score
							last_tf_position[(tf,)] = (set(range(motif_start,motif_end)),score)
					else:
						enhancer_tf_scores[(TF.index(tf),)] += score
						last_tf_position[(tf,)] = (set(range(motif_start,motif_end)),score)
	f.close()
	for k,v in enhancer_tf_scores.items():
		idx_enhancer.append(e_idx)
		idx_tf.append(k[0])
		scores.append(v)
	e_idx += 1
	idx_tissue = [0]*len(idx_tf)
	return (idx_enhancer,idx_tf,idx_tissue),scores,e_idx+1

def invert_core(core,mode):
	u,s,v = np.linalg.svd(core.unfold(mode),full_matrices=False)
	return s,np.dot(np.dot(v.T,np.diag(s**-1)),u.T)

def tucker_tissue(subs,vals,shape,InvertedUnfoldedCore,E,TF,Gene):
	S = sk.sptensor(subs,vals,shape=shape,dtype=np.float64)
	A = S.ttm(E.T,0).ttm(TF.T,1).ttm(Gene.T,3).unfold(2)
	return np.dot(InvertedUnfoldedCore.T,A.T)

def tissue_diff_by_mode(core,U,mode=1,tissue_mode=2):
	# U = [E,TF,TissueDiff,Gene]
	# U[tissue_mode].shape should be eg (1,20)
	remaining_modes = np.setdiff1d(range(len(U)),(mode,tissue_mode))
	Dshape = [U[r].shape[0] for r in remaining_modes]
	smallest_mode = (remaining_modes[np.argmin(Dshape)],np.argmin(Dshape))
	remaining_modes = np.setdiff1d(remaining_modes,[smallest_mode[0]])[0]
	idx = [slice(0,s) for s in core.shape]
	idx[mode] = 0
	idx[tissue_mode] = 0
	D = np.zeros(U[mode].shape[0])
	core_tissue = core.ttm(U[tissue_mode],tissue_mode)
	for i in xrange(U[mode].shape[0]):
		X = core_tissue.ttm(U[mode][i:i+1],mode)
		X = X[idx]
		if smallest_mode[0] > remaining_modes:
			X = X.T
		X = np.dot(X,U[remaining_modes].T).T
		for j in xrange(Dshape[smallest_mode[1]]):
			d = np.dot(X,U[smallest_mode[0]][j:j+1].T)
			D[i] += difference_measure(d)
	return D

def difference_measure(d):
	return d.sum()

# $ awk '{print $1,$2,$3}' E003.narrowPeak | awk -v OFS="\t" '$1=$1' | sed 's/$/\tPOU5F1/g' > E003.narrowPeak.labeled
# $ bedtools intersect -a /ahg/regevdata/users/bcleary/MotifAnalysis/data/enhancer_bed/all_regions.bed -b E003.narrowPeak.labeled -loj > E003.narrowPeak.labeled.enhancerOverlap
if __name__ == "__main__":
	try:
		opts, args = getopt.getopt(sys.argv[1:],'b:f:t:i:o:m:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-b'):
			bed_file = arg
		elif opt in ('-f'):
			tf_file = arg
		elif opt in ('-t'):
			tissue_file = arg
		elif opt in ('-i'):
			infiles = arg
		elif opt in ('-m'):
			m = int(arg)
		elif opt in ('-o'):
			outfile = arg
	#TF = tf_list(tf_file)
	#Tissues = tissue_ids(tissue_file)
	#subs,vals,num_enhancers = parse_bed(bed_file,TF)
	#print 'shape of tensor:',(num_enhancers,len(TF),len(Tissues))
	#print 'max indices:',max(subs[0]),max(subs[1]),max(subs[2])
	#print 'number of elements:',len(vals)
	
	core = sk.dtensor(np.load(infiles + '.core.npy'))
	E = np.load(infiles + '.enhancers.npy')
	TF = np.load(infiles + '.tfs.npy')
	Tiss = np.load(infiles + '.tissues.npy')
	Genes = np.load(infiles + '.genes.npy')
	
	TissLabels = np.load(infiles + '.tissue_labels.npy')
	
	#s,core = invert_core(core,m)	
	#T = tucker_tissue(subs,vals,(num_enhancers,len(TF),1),G,E,TF)
	#d = distance.cdist(T.T,Tiss,'cosine')
	
	T = np.zeros((1,Tiss.shape[1]))
	T[0] = Tiss[0]
	Tdiff = T - Tiss[m:m+1,:]
	TFdiff = tissue_diff_by_mode(core,[E,TF,Tdiff,Genes])
	np.save(infiles + '.TFdiff.E003.%s.npy' % TissLabels[m],TFdiff)
	
	#np.save(outfile + '.lmbda.npy',s)
	#np.save(outfile + '.tissue.npy',T)
	#np.save(outfile + '.tissue_distance.npy',d)
	
	#print np.array(Tissues)[np.argsort(d[0])]
	#print d[0][np.argsort(d[0])]



