# coding=utf-8
# ----------------------------------------------------------------
# 	corrMatrix:
#
# 	generates correlation matrices of:
# 		-diagonal (with samen and different diagonal elements)
# 		-block (constant correlation matrix)
# 			different blocks on diagonal with different covariance values,
# 			but the same within a block
# 		-lowrank:
# 			created via vector multiplication of independent vectors.
# ----------------------------------------------------------------

#### PACKAGES ####
import os
import sys

cwd = os.getcwd()

## NUMERICS
import numpy as np
## RANDOM
from numpy import random

## BLOCK-STRUCTURED COVARIANCE
sys.path.append(cwd + '/code/01_simulation/datageneration')
import lowrank_block as lb
###################

#### FUNCTIONS #####
# diagonal matrix
def corr_diag(dim=None, const=None):
	
	assert dim is not None, 'corr_diag: size of matrix needs to be given'
	
	ident = np.eye(dim)

	if const is None:
		print('corr_diag: no constant for diagonal values given. Set to 1.')
		random.seed(dim)
		const = random.rand(dim) * 10
	
	const = np.array(const, ndmin=1)
	if const.shape[0] is not dim:
		const = np.tile(const, dim)
		const = const[:dim]
	
	diag = ident * const
	
	return diag


#block matrix
def corr_block(dim=None, blocksize=None, rho=[0.9], const=None):
	
	assert dim is not None, 'corr_block: size of matrix needs to be given'
	if blocksize is None:
		blocksize = np.array([dim/2, dim/2], ndmin=1)
	else:
		blocksize = np.array(blocksize, ndmin=1)
		assert sum(blocksize) == dim, 'corr_block: blocksize needs to be a positives integers summing up to the dimension dim'
		
	rho = np.array(rho, ndmin=1)
	
	assert all(rho >= 0) and all(rho <= 1), 'corr_block: rho needs to be between zero and one'
	
	if rho.shape[0] is not dim:
		rho = np.repeat(rho, dim)
		
		# in case rho is not of length one, the vector needs to be cut
		rho = rho[:dim]
	
	# matrix to be filled with block
	block = np.zeros((dim, dim))
	
	zeros = np.where(blocksize == 0)
	for zero in zeros:
		blocksize = np.delete(blocksize, zero)
	
	nblock = blocksize.shape[0]
	
	#
	if const is None:
		print 'corr_block: no constant for blocks given. Set to 1'
		const = np.ones(nblock)
	else:
		const = np.array(const, ndmin=1)
		
	if const.shape[0] is not nblock:
		if const.shape[0] < nblock:
			rep = sum(blocksize)
			const = np.repeat(const, rep)
		else:
			rep = blocksize
			const = np.repeat(const, rep)
		
		# of const is not of length one, the vector needs to be cut
		const = const[:nblock]


	ind1 = 0
	for n in range(nblock):
		
		bs = blocksize[n]
		
		mat = np.ones((bs, bs))
		mat *= rho[n]
		diag = (1 - rho[n]) * np.eye(bs)
		mat += diag
		mat *= const[n]
		
		# write block into matrix
		ind2 = sum(blocksize[:n+1])
		block[ind1:ind2, ind1:ind2] = mat
		
		ind1 = ind2
	
	return block
	

# create random vector
def random_vec(dim, const=10):
	vec = np.floor(random.rand(dim)*2*const - const)
	
	return vec


# create vector independent to vecs
def indep_vec(dim, vecs=None):
	
	# set boolean for whether all vectors are independent
	indep_bool = False
	
	while not indep_bool:
		#create a new vector:
		vec_new = random_vec(dim)
		
		# add new vecotr to other vectors
		if vecs is not None:
		
			# add new vector to other vectord
			vecs_new = np.c_[vecs, vec_new]
	
			# check independence of m vectors
			# deteminant of mxm matrix must be non-zero
			if np.abs( np.linalg.det( np.dot(vecs_new.T, vecs_new) ) ) >  1e-2:
			
				#return true if independent
				indep_bool = True
			
		else:
			vecs_new = vec_new
			indep_bool = True
	
	return vecs_new



# lowrank matrix
def corr_lowrank(dim, rank=1, block=False):
	
	min_eig = -1
	step_all=0
	mat = np.eye(dim)
	
	np.random.seed(rank*dim)
	while min_eig < 0 and step_all < 1000:
		
		if block:
			blocksize = np.ones(rank)
			for i in range(rank-1):
				blocks_len = sum(blocksize)
				size = np.random.randint(1, dim - blocks_len + 1)
				blocksize[i] = size
				
			blocksize[rank-1] = dim - (blocks_len - 1) - size + 1
			blocksize = [int(elem) for elem in blocksize]
			assert sum(blocksize) == dim, 'blocksize has to sum up to %s' %dim
			mat = lb.covar_block(dim=dim, block_dims=blocksize, id=False)
			
		else:
			for r in range(rank):
				if r == 0:
					vecs = indep_vec(dim=dim, vecs=None)
				else:
					vecs = indep_vec(dim=dim, vecs=vecs)
		
			mat = np.dot(vecs, vecs.T)
		
		# make sure the eigenvalues being zero are not negative due to floating poi	nt numbers introducing truncation errors
		min_eig = np.min(np.real(np.linalg.eigvals(mat)))
		if min_eig < 0:
			mat -= 10*min_eig * np.eye(*mat.shape)
			min_eig = np.min(np.real(np.linalg.eigvals(mat)))
		
		step_all +=1
		
	assert np.linalg.matrix_rank(mat) == rank, 'rank is not %d' %rank
	
	return mat
######################
	
