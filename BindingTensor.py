import sys,getopt
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