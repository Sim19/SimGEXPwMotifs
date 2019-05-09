# coding=utf-8

'''
This script contains the function beta_rand and ycondbeta_rand
to randomly generate data from \beta and y|\beta distribution.
Additionally added is the function ycondbeta_seed to set a seed in ycondbeta_rand.

ycondbeta_rand makes use of the fact that the variation among genes is independent
of genes and only dependent on variation induced by different conditions
'''

#### PACKAGES ####
import numpy as np
##################


#### FUNCTIONS ####
def beta_rand(mean, cov, size=1):
	
	# draw from a multivariate normal size-times
	beta = np.random.multivariate_normal(mean=mean, cov=cov, size=size)
	
	if size is None:
		beta = beta.reshape((beta.shape[0], 1))
		
	return beta

# @profile
def ycondbeta_rand(coeff, beta, cov, mu_G=None, mu_C=None):
	'''
	draw random data from y|beta with E(y|beta) = mu_G + mu_C + coeff x beta
	:param coeff:
	:param beta:
	:param cov:
	:param mu_G:
	:param mu_C:
	:return: randomly generated y|beta
	'''
	
	# compute the mean
	# 	coeff x beta
	coeff_beta = np.dot(coeff, beta)

	# 	control shape
	shape_coeff_beta = coeff_beta.shape
	
	if hasattr(shape_coeff_beta, '__len__'):
		if shape_coeff_beta.__len__() == 2:
			G, C = shape_coeff_beta
			# print 'G: %d, C: %d' %(G,C)
		elif shape_coeff_beta.__len__() == 1:
			shape_beta = beta.shape
			shape_coeff = coeff.shape
			if hasattr(shape_beta, '__len__'):
				if shape_beta.__len__() == 1:
					T_beta = shape_beta[0]
				elif shape_beta.__len__() == 2:
					T_beta, C_beta = shape_beta
		
			if hasattr(shape_coeff, '__len__'):
				if shape_coeff.__len__() == 1:
					T_coeff = shape_coeff[0]
					G_coeff = 1
				elif shape_coeff.__len__() == 2:
					G_coeff,T_coeff = shape_coeff
			# print 'T_coeff: %d, G_coeff: %d' %(T_coeff,G_coeff)
			if not T_beta == T_coeff:
				raise ValueError("coefficients and beta have weird dimensions.")
			G, C = G_coeff, C_beta
	else:
		G = C = 1
	
	# 	mu
	mu = np.zeros((G, C))
	if mu_G is not None:
		if not mu_G.shape == (G, C):
			if hasattr(mu_G.shape, '__len__'):
				if mu_G.shape.__len__() == 1 and mu_G.shape[0] == G:
					mu_G = np.kron( mu_G, np.ones((C, 1))).T
				else:
					raise ValueError("mu_G must have dimension G nor GxC")
			else:
				mu_G = mu_C * np.ones(G,C)
			
		# mu = mu + mu_G
		
	if mu_C is not None:
		if not mu_C.shape == (G, C):
			if hasattr(mu_C.shape, '__len__'):
				if mu_C.shape.__len__() == 1 and mu_C.shape[0] == C:
					mu_C = np.kron( np.ones((G, 1)), mu_C)
				else:
					raise ValueError("mu_C must have dimension C or GxC")
			else:
				mu_C = mu_C * np.ones(G,C)
				
		mu = mu + mu_C
		
	mean = mu + coeff_beta
	
	# y - constant variation in genes assumed
	if mean.shape[1] is not 1:
		mean = mean.reshape(-1, order='F')
	
	if cov.shape[0] == mean.shape[0] and cov.shape[1] == mean.shape[0]:
		y = np.random.multivariate_normal(mean=mean, cov=cov)
		
	else:
		raise ValueError('mean is of dimension %d x %d and'
						 'cov is of dimension %d x %d but '
						 'must be of dimension GC x GC'
						 % (mean.shape[0], mean.shape[1],
						 cov.shape[0], cov.shape[1]))
	
	return y

def ycondbeta_seed(seed, coeff, beta, cov, mu_G=None, mu_C=None):
	
	'''
	Function to call  ycondbeta_rand with a seed
	Aim of this function is to be called repetitively to generate large amount of
	generated data with repetition
	:param seed: seed
	:param coeff:
	:param beta:
	:param cov:
	:param mu_G:
	:param mu_C:
	:return: generated data from ycondbeta_rand with seed set to seed
	'''
	
	np.random.seed(seed)
	
	y = ycondbeta_rand(coeff=coeff, beta=beta, cov=cov, mu_G=mu_G, mu_C=mu_C)
	
	return y
####################
