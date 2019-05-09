# coding=utf-8
'''
functions for initializing delta and sigmasquared for datageneration
'''


#### PACKAGES ####
import numpy as np
import numpy.linalg as la
##################


#### FUNCTIONS ####
def eigVal(matrix):
	eig = la.eigvals(matrix)
	return eig


def tweakEigval(V_eig, Sigma_eig, M_eig, delta, frac, method=None):
	
	# 	compute eigenvalue of M \otimes V
	M_V_eig = np.array([ m*v for m in M_eig for v in V_eig])

	# compute constant factor
	if frac == 0:
		return 0
	elif frac == 1:
		return np.sum(Sigma_eig) * delta
	else:
		const = frac/(1-frac) * delta
	
	# run over different methods:
	if method is None or method == 'mean_V_mean_M':
		eigs =  np.mean(Sigma_eig) / (np.mean(V_eig) * np.mean(M_eig))
	elif method == 'mean_V_M':
		eigs =  np.mean(Sigma_eig) / (np.mean(M_V_eig))
	elif method == 'median_V_M':
		eigs = np.median(Sigma_eig) / (np.median(M_V_eig))
	elif method == 'median_V_median_M':
		eigs = np.median(Sigma_eig) / (np.median(V_eig) * np.median(M_eig))
	elif method == 'trace':
		eigs = np.sum(Sigma_eig) / np.sum(V_eig)
	elif not isinstance(method, basestring):
		raise ValueError('method is not a string')
	else:
		raise ValueError('method %s does not exist' %(method))
	
	# compute sigma_beta2 as product of constant and eigenvalues
	sigma_beta2 = const * eigs
	
	return sigma_beta2.real


def checkRange(rsq, frac=0.2):
	
	if (frac-0.01) < rsq < (frac+0.01):
		return True
	else:
		return False

#########################
