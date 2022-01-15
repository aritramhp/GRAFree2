#!/usr/bin/env python
import math
import pickle
import Bio
from Bio import SeqIO


# -----------------------------------------------
# Count the number of occurrence of nucleotides
# -----------------------------------------------
def CountNucleotide(sequence):
	nu = ['A', 'T', 'G', 'C', 'N']
	nu_count = dict()
	for i in nu:
		try:
			count_i = sequence.count(i)
			nu_count.update({i: count_i})
		except KeyError:
			print('Error occurred during nucleotide count.')
	return nu_count


# Count the number of occurrence of each di nucleotides
def CountDinucleotide(sequence):
	di_nu = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC', 'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
	dnu_count = dict()
	for i in di_nu:
		try:
			count_i = sequence.count(i)
			dnu_count.update({i: count_i})
		except KeyError:
			print('Error occurred during di-nucleotides count.')
	return dnu_count


#-----------------------------------------------
# Count the number of occurrence of each codons
#-----------------------------------------------
def CountCodon(sequence):
	codon = ['AAA', 'AAT', 'AAG', 'AAC', 'ATA', 'ATT', 'ATG', 'ATC', 'AGA', 'AGT', 'AGG', 'AGC',
	         'ACA', 'ACT', 'ACG', 'ACC', 'TAA', 'TAT', 'TAG', 'TAC', 'TTA', 'TTT', 'TTG', 'TTC',
	         'TGA', 'TGT', 'TGG', 'TGC', 'TCA', 'TCT', 'TCG', 'TCC', 'GAA', 'GAT', 'GAG', 'GAC',
	         'GTA', 'GTT', 'GTG', 'GTC', 'GGA', 'GGT', 'GGG', 'GGC', 'GCA', 'GCT', 'GCG', 'GCC',
	         'CAA', 'CAT', 'CAG', 'CAC', 'CTA', 'CTT', 'CTG', 'CTC', 'CGA', 'CGT', 'CGG', 'CGC',
	         'CCA', 'CCT', 'CCG', 'CCC']
	codon_count = dict()
	for i in codon:
		try:
			count_i = sequence.count(i)
			codon_count.update({i: count_i})
		except KeyError:
			print('Error occurred during codon count.')
	return codon_count


# -----------------------------------------------
# Compute entropy
#------------------------------------------------
def ComputeEntropy(sequence):
	dinuc_count = CountDinucleotide(sequence)
	codon_count = CountCodon(sequence)

	entropy = []
	for i in range(len(sequence)):
		if i >= 2:
			cdn = sequence[i - 2:i + 1]
			din = sequence[i - 1:i + 1]
			e = math.log2(codon_count[cdn] / dinuc_count[din]) * (-1)
		else:
			e = 1
		entropy.append(e)
	return entropy

#---------------------------------------------------------------------------------------
'''
	GRAPHICAL-FOOT-PRINT OF DNA
	This function (proposed by JM) read the sequence and
	the pointer move RIGHT -- if it read 'G'
	the pointer move LEFT -- if it read 'A'
	the pointer move UP -- if it read 'C'
	the pointer move DOWN -- if it read 'T'
'''
#---------------------------------------------------------------------------------------
def GraphicalFootPrint(sequence, case):	
	sequence = str(sequence).replace('N','').replace('Y','').replace('K','').replace('R','').replace('W','').replace('M','').replace('S','').replace('D','').replace('H','')
	seq_len = len(sequence)
	gfp_xn = []
	gfp_yn = []
	xn=[0]
	yn=[0]
			
	# We proposed Entropy + weighted vector (Modified case-1)
	if case == '1':
		nu_count = CountNucleotide(sequence)
		A_ratio = nu_count['A']*1.0 / seq_len
		C_ratio = nu_count['C']*1.0 / seq_len
		T_ratio = nu_count['T']*1.0 / seq_len
		G_ratio = nu_count['G']*1.0 / seq_len
		entropy = ComputeEntropy(sequence)			
		for i in range(seq_len):				
			if sequence[i] == 'G':
				xn.append(xn[-1] + entropy[i])
				yn.append(yn[-1] + G_ratio)
			elif sequence[i] == 'A':
				xn.append(xn[-1] - entropy[i])
				yn.append(yn[-1] - A_ratio)
			elif sequence[i] == 'C':
				xn.append(xn[-1] - C_ratio)
				yn.append(yn[-1] + entropy[i])
			elif sequence[i] == 'T':
				xn.append(xn[-1] + T_ratio)
				yn.append(yn[-1] - entropy[i])
			
	# We proposed Entropy + weighted vector (Modified case-2)
	elif case == '2':
		nu_count = CountNucleotide(sequence)
		A_ratio = nu_count['A']*1.0 / seq_len
		C_ratio = nu_count['C']*1.0 / seq_len
		T_ratio = nu_count['T']*1.0 / seq_len
		G_ratio = nu_count['G']*1.0 / seq_len
		entropy = ComputeEntropy(sequence)
		for i in range(seq_len):
			if sequence[i] == 'C':
				xn.append(xn[-1] + entropy[i])
				yn.append(yn[-1] + C_ratio)
			elif sequence[i] == 'G':
				xn.append(xn[-1] - entropy[i])
				yn.append(yn[-1] - G_ratio)
			elif sequence[i] == 'T':
				xn.append(xn[-1] - T_ratio)
				yn.append(yn[-1] + entropy[i])
			elif sequence[i] == 'A':
				xn.append(xn[-1] + A_ratio)
				yn.append(yn[-1] - entropy[i])
			
	# We proposed Entropy + weighted vector (Modified case-3)
	elif case == '3':
		nu_count = CountNucleotide(sequence)
		A_ratio = nu_count['A']*1.0 / seq_len
		C_ratio = nu_count['C']*1.0 / seq_len
		T_ratio = nu_count['T']*1.0 / seq_len
		G_ratio = nu_count['G']*1.0 / seq_len
		entropy = ComputeEntropy(sequence)
		for i in range(seq_len):
			if sequence[i] == 'A':
				xn.append(xn[-1] + entropy[i])
				yn.append(yn[-1] + A_ratio)
			elif sequence[i] == 'C':
				xn.append(xn[-1] - entropy[i])
				yn.append(yn[-1] - C_ratio)
			elif sequence[i] == 'T':
				xn.append(xn[-1] - T_ratio)
				yn.append(yn[-1] + entropy[i])
			elif sequence[i] == 'G':
				xn.append(xn[-1] + G_ratio)
				yn.append(yn[-1] - entropy[i])
	
	gfp = [[x,y] for x, y in zip(xn,yn)]		
	return gfp