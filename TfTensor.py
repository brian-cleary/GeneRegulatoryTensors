import sys,getopt
# use scipy-0.14.0
sys.path.insert(0,'/home/unix/bcleary/src/lib/python2.7/site-packages/')
import numpy as np
import sktensor as sk
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

# $ find . -name "chr*.txt" | xargs -n 1 tail -n +2 | awk '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32}' | sort -k1,1 -k2,2n | awk -v OFS="\t" '$1=$1' > all_peaks.txt
# $ bedtools intersect -a all_regions.bed -b all_peaks.txt -loj -sorted > regions.peaks.txt
# format: chr start end blank blank blank chr start end motif correlation scores-->|
def parse_bed(fp,TF,score_thresh=.75):
	f = open(fp)
	enhancer_labels = {}
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
				enhancer_labels[last_enhancer] = e_idx
				idx_enhancer.append(e_idx)
				idx_tissue.append(k[1])
				idx_tf.append(k[0])
				scores.append(v)
			enhancer_tf_scores = defaultdict(float)
			last_enhancer = enhancer
			e_idx += 1
		if ls[6] != '.':
			for i in range(11,len(ls)):
				score = float(ls[i])
				if score > score_thresh:
					tf = ls[9][:ls[9].index('_')]
					motif_start = int(ls[7])
					motif_end = int(ls[8])
					if (tf,i - 11) in last_tf_position:
						last_pos,last_score = last_tf_position[(tf,i - 11)]
						if (len(set(range(motif_start,motif_end)) & last_pos) > 0):
							if score > last_score:
								enhancer_tf_scores[(TF.index(tf),i - 11)] += score - last_score
								last_tf_position[(tf,i - 11)] = (set(range(motif_start,motif_end)),score)
						else:
							enhancer_tf_scores[(TF.index(tf),i - 11)] += score
							last_tf_position[(tf,i - 11)] = (set(range(motif_start,motif_end)),score)
					else:
						enhancer_tf_scores[(TF.index(tf),i - 11)] += score
						last_tf_position[(tf,i - 11)] = (set(range(motif_start,motif_end)),score)
	f.close()
	for k,v in enhancer_tf_scores.items():
		enhancer_labels[last_enhancer] = e_idx
		idx_enhancer.append(e_idx)
		idx_tissue.append(k[1])
		idx_tf.append(k[0])
		scores.append(v)
	return (idx_enhancer,idx_tf,idx_tissue),scores,enhancer_labels

#$ bedtools window -a ENCODE/tss.unique.bed -b all_regions.bed -w 10000 > ENCODE/tss.unique.enhancerOverlap.bed
def gene_inflation_idx(tss_overlap_file,enhancer_labels,gene_labels):
	EtoG = defaultdict(list)
	f = open(tss_overlap_file)
	for line in f:
		ls = line.strip().split('\t')
		if ls[5] != '.':
			if (ls[5],int(ls[6])) in enhancer_labels:
				EtoG[enhancer_labels[(ls[5],int(ls[6]))]].append((gene_labels[ls[4]],ls[4]))
	f.close()
	return EtoG

def inflate_to_genes(subs,vals,EtoG,GeneTissueStates,thresh=0):
	idx_enhancer = []
	idx_tissue = []
	idx_tf = []
	idx_gene = []
	scores = []
	for i,score in enumerate(vals):
		for g_idx in EtoG[subs[0][i]]:
			gene_state_score = GeneTissueStates[subs[2][i]][g_idx[1]]
			new_score = tf_gene_activity_score(score,gene_state_score)
			if new_score > thresh:
				idx_enhancer.append(subs[0][i])
				idx_tf.append(subs[1][i])
				idx_tissue.append(subs[2][i])
				idx_gene.append(g_idx[0])
				scores.append(new_score)
	return (idx_enhancer,idx_tf,idx_tissue,idx_gene),scores

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

def tf_gene_activity_score(binding_score,gene_state_score,kd=7.,q=4):
	if gene_state_score < 10:
		return binding_score*1./gene_state_score**4/(1./kd**q + 1./gene_state_score**q)
	else:
		return 0.
			

# TODO: CLEANUP PATCH TO sktensor THAT ALLOWS FOR MEMORY-EFFICIENT TUCKER (MET) - no intermediate blowups
# MODIFIED FILES: core.py, sptensor.py, tucker.py
def tensor_decomposition(subs,vals,rank,shape,num_proc):
	S = sk.sptensor(subs,vals,shape=shape,dtype=np.float64)
	return sk.tucker.hooi(S,rank,edims=[0,3],cpus=num_proc)

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
		opts, args = getopt.getopt(sys.argv[1:],'b:f:t:r:s:c:o:n:',["outfile="])
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
		elif opt in ('-r'):
			rank = [int(s) for s in arg.split(',')]
		elif opt in ('-s'):
			tss_file = arg
		elif opt in ('-c'):
			chrom_path = arg
		elif opt in ('-n'):
			num_proc = int(arg)
		elif opt in ('-o'):
			outfile = arg
	TF = tf_list(tf_file)
	Tissues = tissue_ids(tissue_file)
	GeneTissueStates = [gene_chromhmm_states('%s/%s.tss.bed' % (chrom_path,t)) for t in Tissues]
	gene_labels = dict([(k,i) for i,k in enumerate(GeneTissueStates[0].keys())])
	subs,vals,enhancer_labels = parse_bed(bed_file,TF)
	EtoG = gene_inflation_idx(tss_file,enhancer_labels,gene_labels)
	subs,vals = inflate_to_genes(subs,vals,EtoG,GeneTissueStates)
	del EtoG
	Labels = [enhancer_labels,TF,Tissues,gene_labels]
	del enhancer_labels
	del gene_labels
	newSubs = []
	newLabels = []
	for i in range(len(subs)):
		ns,l = reduce_to_nonzero_indices(subs[i],Labels[i])
		newSubs.append(ns)
		newLabels.append(l)
	del Labels
	newSubs = tuple(newSubs)
	shape = tuple([len(l) for l in newLabels])
	print 'shape of tensor:',shape
	print 'max indices:',[max(s) for s in newSubs]
	print 'number of elements:',len(vals)
	print 'rank:',rank
	core,U = tensor_decomposition(newSubs,vals,rank,shape,num_proc)
	np.save(outfile + '.core.npy',core)
	np.save(outfile + '.enhancers.npy',U[0])
	np.save(outfile + '.tfs.npy',U[1])
	np.save(outfile + '.tissues.npy',U[2])
	np.save(outfile + '.genes.npy',U[3])
	np.save(outfile + '.enhancer_labels.npy',newLabels[0])
	np.save(outfile + '.tf_labels.npy',newLabels[1])
	np.save(outfile + '.tissue_labels.npy',newLabels[2])
	np.save(outfile + '.gene_labels.npy',newLabels[3])


	