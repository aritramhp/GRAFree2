import argparse
import os
import copy
from os import listdir
from os.path import isfile, join
import time
import math
import Bio
from Bio import SeqIO

import FeatureVector


parser = argparse.ArgumentParser(description='GRAFree2 for numerical representation of nucleotide sequence')
parser.add_argument('--dataset', type=str, default='sample_execution/data', required=True, help="dataset path")
parser.add_argument('--outdir', type=str, default='sample_execution/grafree2', required=True, help="output location")

#--------------------------------------------------------------------------
# Largest sequence
#--------------------------------------------------------------------------
def FindMaxFrag(input_seqfile):
	max_f = 0

	for species in SeqIO.parse(input_seqfile, "fasta"):
		sequence=species.seq
		seq_len = len(sequence)
		frag_size = round(2*math.log(seq_len, 4)+.5)
		no_of_f = round(seq_len/frag_size + 0.5)
		max_f= max(max_f, no_of_f)

	return round(max_f+0.5)



def main(input_seqfile, gfp_outdir, dft_outdir, feat_outdir):
	CASES = ['1','2','3']
	dir_str = '_'.join(CASES)
	
	REDO = True
	WRITE = True
	
	max_no_of_frag = FindMaxFrag(input_seqfile)

	list_of_frag = [max_no_of_frag]

	# Derive the feature vector for each species
	for fragment in list_of_frag:
		feature_vector = FeatureVector.Derive(input_seqfile, CASES, fragment, gfp_outdir, dft_outdir, feat_outdir, WRITE)
			


# ---------------------------------------------------------------------------------------
if __name__ == '__main__':
	args = parser.parse_args()
	DATASET_DIR = args.dataset
	OUT_DIR = args.outdir

	dataset_list = os.listdir(DATASET_DIR)
	print('dataset_list:', int(len(dataset_list)))

	for seqfile in dataset_list:
		# GFP directory
		gfp_outdir = OUT_DIR + '/GFP/' + str(seqfile)
		if os.path.isdir(gfp_outdir) == False:
			mkdr_cmd = 'mkdir -p ' + gfp_outdir
			os.system(mkdr_cmd)
		# Drift directory
		drift_outdir = OUT_DIR + '/Drift/' + str(seqfile)
		if os.path.isdir(drift_outdir) == False:
			mkdr_cmd = 'mkdir -p ' + drift_outdir
			os.system(mkdr_cmd)
		# FeatVect directory
		featvec_outdir = OUT_DIR + '/FeatVect/' + str(seqfile)
		if os.path.isdir(featvec_outdir) == False:
			mkdr_cmd = 'mkdir -p ' + featvec_outdir
			os.system(mkdr_cmd)

		print('#seqfile: ', seqfile)
		input_seqfile = DATASET_DIR + '/' + str(seqfile)
		main(input_seqfile, gfp_outdir, drift_outdir, featvec_outdir)
