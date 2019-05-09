#----------------------
# 	from \beta:
#----------------------
''' MODEL APPLICATFION
	draw vec(E) data by drawing from beta and drawing Y from y|beta
	
	1) vec(beta) \sim N ( 0, vec( sigma_beta * I_TF \otimes V))
	2) vec(y)|vec(beta) \sim N ( vec(mu_G + mu_C + M^TF beta), \delta I_G \otimes \Sigma )
'''


#### PACKAGES ####
import os
cwd=os.getcwd()

import numpy as np
import pandas as pd
import pickle
import optparse

from multiprocessing import Pool, cpu_count

from functools import partial

import sys
sys.path.append(cwd + '/code/01_simulation/datageneration/')
sys.path.append(cwd + '/code/funcs/')

## DATA-class
import Ybetaparam as Ybp

## PARAMETER TWEAKING FOR DATA GENERATION
import var_explained as var
import tweak_param as tp
import params as pm

## TOOLS
import grepFiles as gf
##################


#### FUNCTIONS ####
def genData(r, Cond, Genes, TF, coeff, cov_motif, cov_noise, offset_beta=None):
	'''
	Generate motif influence \beta and gene expression Y according to model
	:param r (int):	numer of repetititons
	:param Cond (int):	samples/condition size
	:param Genes (int): gene size
	:param TF (int):	beta size
	:param coeff (np.array):	coefficient matrix (here: motif scores)
	:param cov_motif (np.array):	covariance matrix of beta
	:param cov_noise (np.array): 	covariance matrix of noise
	:param offset_beta (int or pd.DataFrame): constant variation parameter of Transcription Factors (sigma_beta^2)
	:return: dictionary with gene expression data Y and motif influence beta
	'''
	
	# set seed
	np.random.seed(r)
	
	if isinstance(offset_beta, pd.core.frame.DataFrame):
		cov_motif = offset_beta[r] * cov_motif
	
	# generate beta for repitition r: draw TF times from C-dim distribution
	beta_r = np.random.multivariate_normal(mean=np.zeros(Cond),
											   cov=cov_motif,
											   size=TF)
	# Y | beta
	#	assumption: constant variation in Genes
	# set seed  different to seed before
	np.random.seed(r + 100)
	cov_noise = cov_noise

	noise_r = np.random.multivariate_normal(mean=np.zeros(Cond),
												cov=cov_noise, size=Genes)
	Y_r = np.dot(coeff, beta_r) + noise_r

	return {'Y_r': Y_r, 'beta_r': beta_r}


def generateY(rep=None, Genes=None, Cond=None, TF=None, coeff=None,
			  offset_beta=None, covar_cond=None, covar_data=None,
			  offset_data = None, frac=None, fracNoise=None):
	'''
	Function to generate randomly generated responses Y based on generated coefficients beta
	:param rep: 		number of repetitions to draw Y (Genes x Cond)
	:param Genes: 		first dimension of Y 	- number of genes
	:param Cond: 		second dimension of Y 	- number of conditions
	:param TF: 			parameter dimension 	- number of Transcription Factors
	:param coeff:  		coefficients in model (Genes x TF) 	- motifs
	:param offset_beta: constant variation parameter  		- sigma_beta^2
	:param covar_cond: 	covariance of conditions (Cond x Cond)	- V_Cond
	:param covar_data: 	covariance of data	(Cond x Cond)		- Sigma_C
	:param offset_data: constant variation parameter of data 	- delta
	:return: 			Y | beta - multivariate normal distributed with N( 0, sigma_beta^2 * V \otimes coeff'coeff.T + delta * Sigma \otimes I_G)
	'''
	
	# precompute some values
	cov_noise = offset_data * covar_data
	cov_motif = offset_beta * covar_cond

	if isinstance(offset_beta, list):
		if offset_beta.__len__() == 1:
			offset_beta = offset_beta[0]
	
	
	#POOLING
	Y_Pool = Pool(max(cpu_count()//2, 1)).map(partial(genData,
										 TF=TF, Genes=Genes,
										 Cond=Cond, coeff=coeff,
										 cov_motif=cov_motif,
										 cov_noise=cov_noise), range(rep))
	
	# initialize dataframes
	Y = pd.DataFrame(columns=range(rep))
	beta = pd.DataFrame(columns=range(rep))
	offset_beta_df = pd.DataFrame(columns=range(rep))
	
	# run pover all repetitions
	for rand in range(Y_Pool.__len__()):
		Ybeta = Y_Pool[rand]
		keys = Ybeta.keys()
		if keys.__len__() != 2:
			raise('Error: Pooling error')
		
		Y[rand] 	= pd.Series(Ybeta['Y_r'].reshape(-1, order='F'))
		beta[rand] 	= Ybeta['beta_r'].reshape(-1, order='F')
	
		# bisection:
		rsqrd = var.varCoeff(Y=Y[rand], A=coeff, b=Ybeta['beta_r'])
		offset_beta_old = offset_beta
		lim_low = 0
		lim_up = 1e+4 * offset_beta_old
		idx = 0
		abs_diff = False
		
		# if no fraction betweee nsignal and noise is set, set to 0.2
		if frac is None:
			frac=0.2
		
		# find right parameters such that signal is $frac
		while not tp.checkRange(rsq=rsqrd, frac=frac) and idx < 1000:
			
			if rsqrd < (frac-0.01):
				if abs_diff:
					lim_up = 1e+1 * lim_up
					
				lim_low = offset_beta_old
				offset_beta_new = (lim_up - offset_beta_old) /2
			else:
				if abs_diff:
					lim_low = 0
				lim_up = offset_beta_old
				offset_beta_new = (offset_beta_old - lim_low)/2
			
			rsqrd_old = rsqrd
			cov_motif_new = offset_beta_new * covar_cond
			Ybeta = genData(rand, Cond=Cond, Genes=Genes, TF=TF, coeff=coeff, 								cov_motif=cov_motif_new, cov_noise=cov_noise)
			
			rsqrd = var.varCoeff(Y=Ybeta['Y_r'].reshape(-1, order='F'), A=coeff, b=Ybeta['beta_r'])
			
			offset_beta_old = offset_beta_new
			abs_diff = abs(rsqrd-rsqrd_old) < 0.01
			idx += 1
			
		# update parameter
		offset_beta_df[rand] = [offset_beta_old]
		
		# reshape updated gene expression and motif influence
		Y[rand] = Ybeta['Y_r'].reshape(-1, order='F')
		beta[rand] = Ybeta['beta_r'].reshape(-1, order='F')
		
	# create dictionary of input parameters
	params = {'G': Genes, 'C': Cond, 'TF': TF, 'sigma_beta2': offset_beta_df,
			  'delta': offset_data, 'V': covar_cond, 'Sigma': covar_data,
			  'motif_TG': coeff, 'frac': frac, 'fracNoise': fracNoise}

	# retrun Ybp object
	return Ybp.Ybetaparam(Y=Y, beta=beta, parameter=params)


def main(params, path, rep, Sig = None, R2method=None, frac=None):
	
	# save dictionaries
	#  make sure directory exists
	dir_name = '/data/simulation/'
	directory = os.path.dirname(path + dir_name)
	if not os.path.exists(directory):
		try:
			os.makedirs(directory)
		except OSError as e:
			if e.errno != errno.EEXISTF:
				raise

	# for key in ['diag']:
	for key in params.V.keys():
		
		print(str(key))
		
		# generate data
		Y = generateY(rep=rep, Genes=params.Genes, Cond=params.Cond, TF=params.TF, coeff=params.motif, offset_beta=params.sigma_beta[key], offset_data=params.delta[key], covar_cond=params.V[key], covar_data=params.Sigma[key], frac=params.frac, fracNoise=params.fracNoise)

		covar_cond_name = key

		# 	Y_gen
		filenameYbeta = gf.grepFiles(Cond=params.Cond, Genes=params.Genes, TF=params.TF, fracNoise=params.fracNoise, genMethod=covar_cond_name, rep=rep, Sigma=Sig, R2method=R2method, frac=frac, data='Ybetaparam')
		
		
		with open(filenameYbeta, 'wb') as f:
			pickle.dump(Y, f)



if __name__ == '__main__':
	
	parser = optparse.OptionParser()
	parser.add_option("-f", dest="file")
	parser.add_option("-r", type="str", dest="rep")
	options, args = parser.parse_args()
	
	FILE = str(options.file)
	FILENAME = cwd + FILE
	
	print(FILENAME)
	
	with open(FILENAME, 'rb') as f:
		param = pickle.load(f)
	
	print('Done reading file')
	
	rep = int(options.rep)
	
	main(params=param, rep=rep, path=cwd)
##################
