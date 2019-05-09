# coding=utf-8

### PACKAGES ####
import numpy as np
#################

#### FUNCTIONS ####
def covar_block(dim, block_dims, id=False):
	
	'''
	Function to create Covariance Matrix with several block structures
	:param dim: dimension of entire cov-matrix
	:param block_dims: dimensions of blocks in array in right order
	:param id: boolen to whether add an identity matrix
	:return: block-matrix
	'''
	start = []
	end = []
	dims = []
	start.append(0)
	numBlocks = block_dims.__len__()
	for i in range(numBlocks):
		dims.append(block_dims[i])
		end.append(start[i] + dims[i])
		Cov_i = covar_block_single(dim=dim, start=start[i], end=end[i])
		if i == 0:
			Cov = Cov_i
		else:
			Cov += Cov_i
		
		if numBlocks > 1 and i > 0:
			for j in range(i):
				Cov_i_edge = covar_block_side(dim,
											  start1=start[i],
											  end1 = end[i],
											  start2 = start[i-j-1],
											  end2 = end[i-j-1])
				Cov += Cov_i_edge

		start.append(start[i] + dims[i])
		
	if id:
		Cov += np.eye(dim)
	
	return Cov
	
	
	
def covar_block_single(dim, start, end, id=False):
	#initialize matrix - set all elements to zero
	cov = np.zeros((dim, dim))
	
	# get dimension of block matrix
	cov[start:end, start:end] = np.ones(end - start)
	if id:
		cov[start:end, start:end] += np.eye(end - start)
	
	return cov


def covar_block_side(dim, start1, end1, start2, end2):
	# initialize matrix - set all elements to zero
	cov = np.zeros((dim, dim))
	
	#get dimensions of block matrix
	cov[start1:end1, start2:end2] = -1e-2 * np.ones((end1 - start1, end2 - start2))
	#make symmetric
	cov[start2:end2, start1:end1] = -1e-2 * np.ones((end2 - start2, end1 - start1))

	return cov
##################