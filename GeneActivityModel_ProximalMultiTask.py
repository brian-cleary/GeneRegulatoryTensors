import sys,getopt
sys.path.insert(0,'/ahg/regev/users/bcleary/src/lib/python2.7/site-packages/')
import numpy as np
from collections import defaultdict
from operator import itemgetter
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import svds
from scipy.linalg import eigh
import sktensor as sk

def tissue_ids(tissue_file):
	f = open(tissue_file)
	T = []
	for line in f:
		T.append(line.strip())
	f.close()
	return T

def ensemble_dict(fp):
	f = open(fp)
	header = f.readline()
	EnsembleIds = {}
	for line in f:
		ls = line.strip().split()
		EnsembleIds[ls[0]] = ls[1]
	f.close()
	return EnsembleIds

def expression_vector(fp,Tissues,EnsembleIds):
	f = open(fp)
	header = f.readline().strip().split()
	Cidx = [header.index(tissueid) for tissueid in Tissues]
	E = [defaultdict(float) for _ in Tissues]
	for line in f:
		ls = line.strip().split()
		gene = EnsembleIds.get(ls[0],None)
		for i,colidx in enumerate(Cidx):
			expression = float(ls[colidx])
			if expression > 0:
				E[i][gene] = expression
	f.close()
	E = [sorted(e.iteritems(),key=itemgetter(1)) for e in E]
	return [dict([(x[0],100.*i/len(e)) for i,x in enumerate(e)]) for e in E]

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

def smooth_expression(e,k=30,q=4):
	return e**q/(k**q + e**q)

def regression_response(ExpressionVector,gene_labels):
	gene_labels = list(gene_labels)
	geneDict = dict([(g,i) for i,g in enumerate(gene_labels)])
	Y = np.zeros((len(ExpressionVector),len(gene_labels)))
	for i,E in enumerate(ExpressionVector):
		for gene,quantile in E.iteritems():
			if gene in geneDict:
				Y[i,geneDict[gene]] = smooth_expression(quantile)
	return Y

def initialize_variables(subs,vals,Y,enh_mode=0,tf_mode=1,tissue_mode=2,gene_mode=3):
	shape = tuple([max(s)+1 for s in subs])
	GeneSubs = [[[],[]] for _ in xrange(shape[gene_mode])]
	GeneVals = [[] for _ in xrange(shape[gene_mode])]
	for i,v in enumerate(vals):
		# tissues x tfs
		g = subs[gene_mode][i]
		GeneSubs[g][0].append(subs[tissue_mode][i])
		GeneSubs[g][1].append(subs[tf_mode][i])
		GeneVals[g].append(v)
	XtX = []
	Steps = []
	XY = []
	for g in xrange(shape[gene_mode]):
		X = coo_matrix((np.array(GeneVals[g]),(np.array(GeneSubs[g][0]),np.array(GeneSubs[g][1]))),shape=(shape[tissue_mode],shape[tf_mode]))
		X = X.tocsr()
		# densifies with the centering here...could probably just work on the nnz though
		#X = X - X.mean(0)
		xtx = X.T.dot(X)
		y = Y[:,g]
		#y = y - np.average(y)
		XY.append(X.T.dot(y))
		XtX.append(xtx)
		u,s,vt = svds(xtx,k=1)
		Steps.append(1./s[0])
	return XtX,XY,np.array(Steps)

def grad_descend(XtX,XY,W,lda1,lda2,alpha,Oinv,Steps,n,tol=10**-4,maxItr=250):
	d = W.shape[0]
	t = W.shape[1]
	coReg = lda1*(1-alpha)*np.eye(t)+lda2*Oinv
	change = 1
	itr = 0
	while (change > tol) and (itr < maxItr):
		w_grad = np.zeros((d,t))
		for i in xrange(t):
			w_grad[:,i] = 2*(XtX[i].dot(W[:,i]) - XY[i])/n
		w_grad += W.dot(coReg)
		W = W - Steps*w_grad
		change = np.linalg.norm(Steps*w_grad)
		itr += 1
	return W

def grad_single(X,Y,tol=10**-5,maxItr=500):
	change = 1
	itr = 0
	W = np.zeros(X.shape[1])
	XtX = X.T.dot(X)
	XY = X.T.dot(Y)
	u,s,vt = svds(XtX,k=1)
	step = 1./s[0]
	while (change > tol) and (itr < maxItr):
		w_grad = 2*(XtX.dot(W)-XY)/Y.shape[0]
		W = W - step*w_grad
		change = np.linalg.norm(step*w_grad)
		print change
		itr += 1
	return W

def proximal_grad(W,delta):
	return (W - delta*W/np.sign(W+10**-10))*(abs(W) > delta)

def update_Oinv(W):
	WtW = W.T.dot(W)
	w,v = eigh(WtW)
	return (w**.5).sum()*v.dot(np.diag(w**-.5)).dot(v.T)

def iterate_prox_relations(subs,vals,Y,lda1,lda2,alpha,tol=10**-4,maxItr=100):
	XtX,XY,Steps = initialize_variables(subs,vals,Y)
	W = np.zeros((XtX[0].shape[0],Y.shape[1]))
	Oinv = np.eye(Y.shape[1])/Y.shape[1]
	change = [1,1]
	itr = 0
	while (np.average(change) > tol) and (itr < maxItr):
		oldW = W
		oldOinv = Oinv
		W = grad_descend(XtX,XY,W,lda1,lda2,alpha,Oinv,Steps,Y.shape[0])
		W = proximal_grad(W,Steps*lda1*alpha)
		Oinv = update_Oinv(W)
		change[0] = np.linalg.norm(W-oldW)
		change[1] = np.linalg.norm(Oinv-oldOinv)
		itr += 1
		print itr,change
		sys.stdout.flush()
	return W,np.linalg.inv(Oinv)

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
		opts, args = getopt.getopt(sys.argv[1:],'e:t:s:c:f:o:',["outfile="])
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
		elif opt in ('-f'):
			ensembleid_file = arg
		elif opt in ('-e'):
			expression_file = arg
	Tissues = tissue_ids(tissue_file)
	EnsembleIds = ensemble_dict(ensembleid_file)
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
	ExpressionValues = expression_vector(expression_file,Tissues,EnsembleIds)
	gene_labels = np.load(outfile + '.bindingDistance.gene_labels.npy')
	Y = regression_response(ExpressionValues,gene_labels)
	W,O = iterate_prox_relations(newSubs,newVals,Y,.2,.01,.75)
	np.save(outfile + '.finalModel.W.npy',W)
	np.save(outfile + '.finalModel.O.npy',O)
	