import sys,getopt
from collections import defaultdict
import numpy as np

def tf_list(fp):
	f = open(fp)
	TF = []
	for line in f:
		TF.append(line.strip())
	f.close()
	return TF

# $ find . -name "chr*.txt" | xargs -n 1 tail -n +2 | awk '{print $2,$3,$4,$1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32}' | sort -k1,1 -k2,2n | awk -v OFS="\t" '$1=$1' > all_peaks.txt
# $ bedtools intersect -a all_regions.bed -b all_peaks.txt -loj -sorted > regions.peaks.txt
# format: chr start end chr start end motif correlation scores-->|
def parse_bed(fp,TF,score_thresh=.65):
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
		if ls[3] != '.':
			for i in range(8,len(ls)):
				score = float(ls[i])
				if score > score_thresh:
					tf = ls[6][:ls[6].index('_')]
					motif_start = int(ls[4])
					motif_end = int(ls[5])
					if (tf,i - 8) in last_tf_position:
						last_pos,last_score = last_tf_position[(tf,i - 8)]
						if (len(set(range(motif_start,motif_end)) & last_pos) > 0):
							if score > last_score:
								enhancer_tf_scores[(TF.index(tf),i - 8)] += score - last_score
								last_tf_position[(tf,i - 8)] = (set(range(motif_start,motif_end)),score)
						else:
							enhancer_tf_scores[(TF.index(tf),i - 8)] += score
							last_tf_position[(tf,i - 8)] = (set(range(motif_start,motif_end)),score)
					else:
						enhancer_tf_scores[(TF.index(tf),i - 8)] += score
						last_tf_position[(tf,i - 8)] = (set(range(motif_start,motif_end)),score)
	f.close()
	for k,v in enhancer_tf_scores.items():
		enhancer_labels[last_enhancer] = e_idx
		idx_enhancer.append(e_idx)
		idx_tissue.append(k[1])
		idx_tf.append(k[0])
		scores.append(v)
	return (idx_enhancer,idx_tf,idx_tissue),scores,enhancer_labels

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
		opts, args = getopt.getopt(sys.argv[1:],'b:f:o:',["outfile="])
	except:
		print help_message
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-b'):
			bed_file = arg
		elif opt in ('-f'):
			tf_file = arg
		elif opt in ('-o'):
			outfile = arg
	TF = tf_list(tf_file)
	subs,vals,enhancer_labels = parse_bed(bed_file,TF)
	np.save(outfile + '.binding.subs.npy',subs)
	np.save(outfile + '.binding.vals.npy',vals)
	np.save(outfile + '.binding.enhancer_labels.npy',enhancer_labels)