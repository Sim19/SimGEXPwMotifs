# coding=utf-8

'''
Script to load motif data for data 01_simulation
'''

#### PACKAGES #####
# NUMERIC
import numpy as np

# DATA FRAMES
import pandas as pd

# SYSTEM/OS
import sys
import os

cwd = os.getcwd()
sys.path.append(cwd + '/code/data_prep/')
import prepMotif4Py as pm
#####################


#### FUNCTIONS #####
def getMotif(Genes=None):
	'''
		motif-file in right format: source script from code/data_prep/prepMotif4Py
		motif: nparray G x C from real data
	'''
	# load necessary files
	if Genes is None or Genes == 978:
		motiffile = cwd + '/data/l1000/motif_lincs_400_100.txt'
		condfile = cwd + '/data/l1000/mcf7_landmark_level5.gctx'
		# applt. script to get motif into right format
		motif = pm.main(motif_file=motiffile, cond_file=condfile)
	else:
		if Genes.__class__ is not int:
			Genes = int(Genes)
		motiffile = cwd + '/data/motif/motif_mostVariable' + str(Genes) + '_400_100.txt'
		motif = pm.read_file(motiffile, header=4)
		motif_index = []
		for idx in range(motif.shape[0]):
			motif_index.append(motif[motif.columns[0]][idx].split(' ')[1])
		# set new index with only Gene Names by Ensembl
		motif.index = motif_index
		# drop column that only contained name and chromosome of gene
		motif = motif.drop(motif.columns[0], 1)
		
	# get motifs for genes under consideration - SHUFFLE!
	motif_min = motif.min(axis=0)
	motif_max = motif.max(axis=0)
	motif -= motif_min
	motif /= (motif_max - motif_min + 1e-4)
	
	return motif


def getTG(motif):
	
	G,T = motif.shape
	
	return (T,G)


def getM(motif):
	
	# compute M by taking the dot-product \sum_t m_tg * m_t`g
	M = np.dot(motif, motif.T)
	M_df = pd.DataFrame(data=M, index=motif.index, columns=motif.index)
	
	return M_df

