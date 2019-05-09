# coding=utf-8

'''
Script to generate V for data in 01_simulation
'''

#### PACKAGES ####

import os, sys

# ANALYSIS
import numpy as np
from numpy import random

# LOAD OWN PACKAGES
cwd = os.getcwd()
sys.path.append(cwd + '/code/01_simulation/datageneration/')

import corrMatrix as cm
###################

#### FUNCTIONS ####
# get parameters C
def createV(C):
	'''
	Create correlation matrices
	:param C (int):	dimension of matrices
	:return: 		list of generated correlation matrices
	'''

	# identical
	V_id = cm.corr_diag(dim=C, const=1)

	# block matrix
	bs = map( lambda x: int(round(x, 0)), np.array([C*.2, C*.5, C*.2, C*.1]))
	i=0
	while sum(bs) != C:
		i += 1
		if sum(bs) < C:
			bs[i] += 1
		elif sum(bs) > C:
			bs[i] -= 1
		
	# lowrank matrices
	# 	put offset such that cov is invertible
	offset = 1e-04*np.eye(C)
	rank = [2, C * 1/2]
	#rank = map(lambda x: x * C * 1/4, [2,3])
	if 2 not in rank:
		rank.append(2)
		rank.sort()

	random.seed(1)		
	V_lr_1 = cm.corr_lowrank(C, rank[0], block=True) + offset
	V_lr_2 = cm.corr_lowrank(C, rank[1], block=True) + offset

	# freeform
	V_ff = np.eye(C) + np.reshape(random.rand(C*C), (C, C)) 	# diagonal matrix + random noise
	V_freeform = np.dot(V_ff, V_ff.T)
	
	dict_lowrank = {
		'lowrank_' + str(rank[0]): V_lr_1,
		'lowrank_' + str(rank[1]): V_lr_2,
		'id': V_id,
		'freeform': V_freeform
		}
	
	return dict_lowrank
#####################
