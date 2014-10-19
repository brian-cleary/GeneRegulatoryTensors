import sys,getopt
# use scipy-0.14.0
sys.path.insert(0,'/home/unix/bcleary/src/lib/python2.7/site-packages/')
import numpy as np
import sktensor as sk
from collections import defaultdict
			
def apply_tf_weights(subs,vals,Wtf,tf_mode=1):
	return vals*Wtf[subs[tf_mode]]

def eliminate_elements(subs,vals,labels,mode,tissue_mode=2,num_tissues=None,minMaxPercentile=0.60):
	if num_tissues == None:
		num_tissues = max(subs[tissue_mode])+1
	num_mode = len(labels)
	TissueElementSums = np.zeros((num_tissues,num_mode))
	for i,v in enumerate(vals):
		TissueElementSums[subs[tissue_mode][i],subs[mode][i]] += v
	TissueElementPercentiles = np.zeros(TissueElementSums.shape)
	for i,t in enumerate(TissueElementSums):
		nnz = np.nonzero(t)[0]
		xs = np.argsort(abs(t[nnz]))
		for rank,idx in enumerate(xs):
			TissueElementPercentiles[i,nnz[idx]] = rank/float(len(nnz))
	MaxPercentiles = TissueElementPercentiles.max(0)
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

# TODO: CLEANUP PATCH TO sktensor THAT ALLOWS FOR MEMORY-EFFICIENT TUCKER (MET) - no intermediate blowups
# MODIFIED FILES: core.py, sptensor.py, tucker.py
def tensor_decomposition(subs,vals,rank,shape):
	S = sk.sptensor(subs,vals,shape=shape,dtype=np.float64)
	return sk.tucker.hooi(S,rank,edims=[0,3])

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
		opts, args = getopt.getopt(sys.argv[1:],'r:o:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-r'):
			rank = [int(s) for s in arg.split(',')]
		elif opt in ('-o'):
			outfile = arg
	subs = np.load(outfile + '.bindingDistance.TFexpanded.subs.npy')
	vals = np.load(outfile + '.bindingDistance.TFexpanded.vals.npy')
	Wtf = np.load(outfile + '.regressionModel.Wtf.100.npy')
	Wtf = Wtf[:,0]
	WtfVals = apply_tf_weights(subs,vals,Wtf)
	del vals
	
	enhancer_labels = np.load(outfile + '.bindingDistance.enhancer_labels.npy')
	newSubs,newVals,newEnhLabels = eliminate_elements(subs,WtfVals,enhancer_labels,0,num_tissues=27)
	del WtfVals
	del subs
	
	tf_labels = np.load(outfile + '.bindingDistance.TFexpanded.tf_labels.npy')
	newSubs,newVals,newTFLabels = eliminate_elements(newSubs,newVals,tf_labels,1,num_tissues=27)
	
	gene_labels = np.load(outfile + '.bindingDistance.gene_labels.npy')
	newSubs,newVals,newGeneLabels = eliminate_elements(newSubs,newVals,gene_labels,3,num_tissues=27)
	
	np.save(outfile + '.finalModel.subs.npy',newSubs)
	np.save(outfile + '.finalModel.vals.npy',newVals)
	np.save(outfile + '.finalModel.enhancer_labels.npy',newEnhLabels)
	np.save(outfile + '.finalModel.tf_labels.npy',newTFLabels)
	np.save(outfile + '.finalModel.gene_labels.npy',newGeneLabels)
	
	shape = tuple([max(s)+1 for s in newSubs])
	newSubs = tuple(newSubs)
	print 'shape of tensor:',shape
	print 'number of elements:',len(newVals)
	print 'rank:',rank
	
	core,U = tensor_decomposition(newSubs,newVals,rank,shape)
	
	np.save(outfile + '.finalModel.core.npy',core)
	np.save(outfile + '.finalModel.enhancers.npy',U[0])
	np.save(outfile + '.finalModel.tfs.npy',U[1])
	np.save(outfile + '.finalModel.tissues.npy',U[2])
	np.save(outfile + '.finalModel.genes.npy',U[3])


	