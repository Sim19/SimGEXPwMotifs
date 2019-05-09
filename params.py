# coding=utf-8
'''
This script generates data necessary for viz_tools
'''

#### PACKAGES ####
# SYSTEM, OS
import os
import sys

# DATA ANALYSIS
from numpy import random
import numpy as np

# IMPORT/EXPORT DATA
import cPickle as pickle

cwd=os.getcwd()

# OWN PACKAGES
sys.path.append(cwd + '/code/01_simulation/datageneration/')

### DATA-GENERATION
import cond_covar as cc
import motif as mf
import params_dict as pd
import tweak_param as tp
####################


#### FUNCTIONS ####
def mu_G(G, C):
	muG = np.zeros(G)
	for i in range(G):
		muG[i] = (2 * random.randn(C)).mean()
	
	return muG


def mu_C(G, C):
	muC = np.zeros(C)
	for j in range(C):
		muC[j] = (2 * random.randn(G)).mean()
	
	return muC


# 	NOISE
def genSigma(C):
	''' Sigma - condition-condition noise matrix
		assumption: low-rank
		1rst step: draw each element from a standard normal distribution
	
	'''
	Sigma0 = random.randn(C, C)  # cond-cond noise matrix
	Sigma = np.dot(Sigma0.T, Sigma0)
	
	return Sigma


def genDelta():
	return random.rand(1)


def main(C, Genes=978, Sig='random', R2method=None, frac=None, TF=None, fracNoise=0.7):
	# conditional covariances
	V = cc.createV(C=C)
	
	# motif data
	motif = mf.getMotif(Genes)
	if TF is not None:
		random.seed(0)
	
		## subset on TF-motifs
		idx = random.randint(0, 623, TF)
		motif = motif.iloc[:, idx]
		
	M = mf.getM(motif)
	
	# dimension of data
	T, G = mf.getTG(motif)
	assert Genes == G, ValueError('the loaded motif file does not correspond to the number of genes indicated')
	# genetic covariance
	I_G = np.eye(G)
	
	# mean
	muG = mu_G(G, C)
	muC = mu_C(G, C)
	
	# Sigma
	if 'random' in Sig:
		Sigma = genSigma(C)
	elif 'identity' in Sig:
		Sigma = np.eye(C)
	else:
		raise ValueError('Noise-form has to contain random or identity')
	
	if frac is None:
		frac=0.2
		
	Sigma = {i: Sigma for i in V.keys()}
	
	
	### Prepare for tweikEigenvalue
	
	# compute eigenvalues
	M_eig = tp.eigVal(M)
	Sigma_eig = [tp.eigVal(Sigma[i]) for i in Sigma.keys()]
	V_eig = [tp.eigVal(V[i]) for i in V.keys()]
	
	# noiseSigma: if Sigma dependent on V, add V with frac=0.7
	if 'V' in Sig:
		#compute factor to multiply V with for
		#	 Noise=gamma*V + delta*Sigma
		gamma = [tp.tweakEigval(V_eig=V_eig[i], Sigma_eig=Sigma_eig[i], M_eig=M_eig, delta=1, frac=fracNoise, method='trace') for i in range(V_eig.__len__())]
		
		gamma_dict = {}
		for n, v in zip(V.keys(), gamma):
			gamma_dict.setdefault(n, []).append(v)
		
		
		Sigma = {i: gamma_dict[i] * V[i] + Sigma[i] for i in Sigma.keys()}
	
	V = {i: V[i]/V[i].trace() for i in V.keys()}
	Sigma = {i: Sigma[i]/Sigma[i].trace() for i in Sigma.keys()}
	
	#### RECOMPUTE TO FIND FACTOR FOR V : Noise = (frac):(1-frac)
	#  delta
	normSigma = [ np.linalg.norm(Sigma[i], ord=2) for i in Sigma.keys() ]
	delta = [ (1 - frac) * 1 / normSigma[i] for i in range(normSigma.__len__())]
	
	# sigma_beta
	sigbeta2 = [tp.tweakEigval(V_eig=V_eig[i], Sigma_eig=Sigma_eig[i], M_eig=M_eig, delta=delta[i], frac=frac, method=R2method) for i in range(V_eig.__len__())]
	
	sigma_beta = {}
	for n, v in zip(V.keys(), sigbeta2):
		sigma_beta.setdefault(n, []).append(v)
	delta_dict = {}
	for n, v in zip(Sigma.keys(), delta):
		delta_dict.setdefault(n, []).append(v)
		
	return pd.params_dict(C=C, G=G, T=T, V=V, Sigma=Sigma, motif=motif, delta=delta_dict, sigma_beta=sigma_beta, I_G=I_G, mu_G=muG, mu_C=muC, frac=frac, TF=TF, fracNoise=fracNoise)



### TEST PROCEDURE - THIS SCRIPT IS NOT CALLED INDIVIDUALLY
if __name__ == '__main__':
	Cond = 4
	Sigma = 'identityV'
	param = main(cwd, C=Cond, Sig=Sigma)
	print('C: %s' %(param.Cond))
	print('Sigma: %s' %param.Sigma)
	FILENAME = cwd + '/data/viz_tools/params_generated_C_' + str(Cond) + '_Sigma_' + str(
		Sigma) + '.pkl'
	with open(FILENAME, 'w') as f:
		pickle.dump(param, f)
#########################