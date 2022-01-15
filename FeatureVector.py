import pickle
import os
import numpy
import math
import json
import Bio
from Bio import SeqIO

import Drift
import DriftAnalyses
import GFP

# ------------------------------------------------------------------
# Derive the feature vector for all the species
# ------------------------------------------------------------------
def Derive(seqfile, CASES, fragment, gfp_outdir, dft_outdir, feat_outdir, WRITE):
	# print('This execution discards the previous stored data!!!')
	# Derive GFP
	print('Deriving GFPs...')
	GFP_dict = DeriveGFP(seqfile, CASES)
	if WRITE:
		WriteData(GFP_dict, gfp_outdir)
			
	# Derive Drift
	print('Deriving Drifts...')
	Drift_dict = DeriveDrift(GFP_dict, CASES)
	if WRITE:
		WriteData(Drift_dict, dft_outdir)
			
	# Derive feature vector
	feature_vector_dict = dict()
	for species in Drift_dict.keys():
		feature_vector = []
		for case in range(len(CASES)):
			drift = Drift_dict[species][case]
			featvect = DriftAnalyses.Analysis(drift, fragment)
			feature_vector.append(featvect)
		try:
			feature_vector_dict.update({species: feature_vector})
		except:
			print("Error in deriving feature vector. Species: ", species)
	if WRITE:
		WriteData(feature_vector_dict, feat_outdir)
	
	# return feature_vector_dict


# ------------------------------------------------------------------
# Derives the GFPs for all the species and stores in a directory
# ------------------------------------------------------------------
def DeriveGFP(seqfile, CASES):
	GFP_dict = dict()
	for record in SeqIO.parse(seqfile, "fasta"):
		GFP_xy = []
		species = record.id
		sequence = record.seq
		for case in CASES:
			gfp = GFP.GraphicalFootPrint(sequence, case)
			GFP_xy.append(gfp)
		try:
			GFP_dict.update({species: GFP_xy})
		except:
			print('Error in computing GFP. Species: ', species)

	return GFP_dict


# ------------------------------------------------------------------
# Derives the Drifts for all the species and stores in a directory
# ------------------------------------------------------------------
def DeriveDrift(GFP_dict, CASES):
	Drift_dict = dict()
	for species in GFP_dict.keys():
		drift_xy = []
		for case in range(len(CASES)):
			GFP_xy = GFP_dict[species][case]
			block_size = ComputeBlockLength(GFP_xy)
			# block_size = 18
			dft = Drift.Drift(GFP_xy, block_size)
			drift_xy.append(dft)
		try:
			Drift_dict.update({species: drift_xy})
		except:
			print('Error in computing Drift. Species: ', species)
	
	return Drift_dict


# ------------------------------------------------------------------
# Store the data (including data structure)
# ------------------------------------------------------------------
def WriteData(data, dirname):
	for k in data.keys():
		with open(dirname + '/' + k + '.json', 'w') as fp:
			json.dump(data[k], fp)




#------------------------------------------------------------------
# Compute block length 
#------------------------------------------------------------------
def ComputeBlockLength(GFP_xy):
	low = 5
	high = 600
	block_size = low
	step = .15
	window = 100
	gap = 10
	flag = False
	while True:
		# print('Inter block: ', block_size)
		drift_xy = Drift.Drift(GFP_xy, block_size)
		ffp_0 = ComputeFFP(drift_xy)
		
		drift_xy = Drift.Drift(GFP_xy, block_size+gap)
		ffp_1 = ComputeFFP(drift_xy)

		change = (ffp_1-ffp_0)*100./ffp_0
		if change >= step:
			if flag:
				window = max(int(window/2), 1)
			block_size = block_size + window
			if block_size > high:
				step *= 2
				gap /= 2
				block_size = low
				window = 100
			flag = False
		elif change < step:
			# print('case-2')
			if window == 1:
				break
			else:
				window = int(window/2)
				block_size = block_size - window
			if block_size < low:
				step /= 2
				gap *= 2
				print(step, ffp_1, ffp_0)
				block_size = low
				window = 100
			flag = True

	return block_size

#------------------------------------------------------------------
# Compute block length 
#------------------------------------------------------------------
def ComputeBlockLength(GFP_xy):
	low = 5
	high = 600
	block_size = low
	step = .15
	window = 100
	gap = 10
	flag = False
	while True:
		# print('Inter block: ', block_size)
		drift_xy = Drift.Drift(GFP_xy, block_size)
		ffp_0 = ComputeFFP(drift_xy)
		
		drift_xy = Drift.Drift(GFP_xy, block_size+gap)
		ffp_1 = ComputeFFP(drift_xy)

		change = (ffp_1-ffp_0)*100./ffp_0
		if change <= 0:
			if block_size+gap >= high:
				break
			else:
				gap = gap*2
		else:
			if change >= step:
				if flag:
					window = max(int(window/2), 1)
				block_size = block_size + window
				if block_size > high:
					step *= 2
					gap /= 2
					block_size = low
					window = 100
				flag = False
			elif change < step:
				# print('case-2')
				if window == 1:
					break
				else:
					window = int(window/2)
					block_size = block_size - window
				if block_size < low:
					step /= 2
					gap *= 2
					block_size = low
					window = 100
				flag = True

	return block_size

#--------------------------------------------------------------------------------------
# Compute FFP from the drift coordinates
#--------------------------------------------------------------------------------------
def ComputeFFP(drift_xy):
	freq_dict = dict()
	sum_count = 0
	for i in range(len(drift_xy)):
		if str(drift_xy[i]) in freq_dict.keys():
			freq_dict[str(drift_xy[i])] += 1
			sum_count += 1
		else:
			freq_dict.update({str(drift_xy[i]) : 1})
			sum_count += 1
	
	FFP_dict=freq_dict.copy()
	for k in FFP_dict.keys():
		FFP_dict[k] = FFP_dict[k]*1.0/sum_count
	
	#print FFP_dict
	ffp = 0
	for k in FFP_dict.keys():
		ffp = ffp + FFP_dict[k]*math.log(FFP_dict[k])
	ffp = ffp *(-1)
	
	return ffp

