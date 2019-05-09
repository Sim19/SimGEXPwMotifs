# coding=utf-8
'''
generate gene expression data with various covariance matrix shapes,
determinded strength of signal-to-noise ration
'''

# 	====== PACKAGES ======
# LIN ALG
import os
import sys

from numpy import set_printoptions
# RANDOMNESS
from numpy.random import RandomState

cwd = os.getcwd()
print(cwd)
# DATA STORAGE
import optparse

# monitor the time for execution of __main__ script
import time
# 	======================


# 	==== OWN PACKAGES ====
sys.path.append(cwd + '/code/01_simulation/datageneration/')
sys.path.append(cwd + '/code/funcs/')

### DATA-GENERATION ###
import params as pms
import data_generation as dg
# 	======================

# 	==== SET OPTIONS =====
''' set printoptions
	determine displt.y of floating point numbers,
	arrays and other numpy object
	first arg: precision: number of digits of precision per
	floating point output (default 8)
'''
set_printoptions(4)

''' set random state
	make results reproducible *Mersenne Twister pseudo-random generator
	first arg: seed - initialize pseudo-random number generator
'''
random = RandomState(1)
# 	======================


def main(Cond, rep, Genes, Sigma, path, R2method, frac, TF=None, fracNoise=None):
	
	print('generate_data: Number of repetitions: ' + str(rep))
	print('generate_data: read parameters')
	params = pms.main(C=Cond, Genes=Genes, Sig=Sigma, R2method=R2method, frac=frac, TF=TF, fracNoise=fracNoise)
	print('generate_data: done!')
	

	# ----------------------------------------------------------------
	# 	GENERATE DATA
	# ----------------------------------------------------------------
	
	print('generate_data: generate Y and beta:')
	dg.main(params=params, path=path, rep=rep, Sig=Sigma,
			R2method=R2method, frac=frac)
	print('generate_data: done with generation of Y and beta')


if __name__ == '__main__':
	
	cwd = os.getcwd()
	
	''' parser options:
		c: number of conditions
		r: number of repititions
	'''
	parser = optparse.OptionParser()
	parser.add_option("-c", type="int", dest="cond")
	parser.add_option("-G", default=978, type="int", dest="Genes")
	parser.add_option('-r', type='int', dest='rep')
	parser.add_option('-f', default=None, type='float', dest='fraction')
	parser.add_option('--fN', default=None, type='float', dest='fracNoise')
	parser.add_option('-S', default=None, type='str', dest='Sigma')
	parser.add_option('-R', default=None, type='str', dest='R2method')
	parser.add_option('-T', default=None, type='int', dest='TF')
	options, args = parser.parse_args()
	# 	======================
	# ----------------------------------------------------------------

	# ----------------------------------------------------------------
	# 	PARAMETER SETTING
	# ----------------------------------------------------------------
	# DEFINE PHENOTYPE
	# 	condition size
	# C0 = 4
	C = int(options.cond)
	Genes = int(options.Genes)
	

	# amount of repetitions to be drawn for beta
	rep = int(options.rep)
	
	# fraction explained by variation in data (not noise)
	if options.fraction is not None:
		frac = float(options.fraction)
	else:
		frac=None
	
	# Form of Sigma
	if options.Sigma is None:
		Sigma = 'random'
	else:
		Sigma = str(options.Sigma)
	
	# R2-method applied
	if options.R2method is not None:
		R2method = str(options.R2method)
	else:
		R2method = None
		
	if options.TF is not None:
		TF = int(options.TF)
	else:
		TF = None
		
	if options.fracNoise is not None:
		fracNoise = float(options.fracNoise)
	else:
		if 'V' in Sigma:
			Warning('fracNoise is not set even though noise is structured. FracNoise is set to 0.7')
			fracNoise = 0.7
		else:
			fracNoise = None
			
	print('generate_data: \n'
		  '	Number of conditions: %s \n'
		  '	Number of genes: %s \n'
	  	  ' fraction of Signal: %s \n'
		  ' fraction in noise: Signal2Noise %s \n'
		  ' Sigma: %s ' %(str(C), str(Genes), str(frac), str(fracNoise), str(Sigma))
	  )
	
	#monitor time of execution of main():
	startTime = time.time()
	main(path=cwd, Cond=C, Genes=Genes, rep=rep, Sigma=Sigma, R2method=R2method, frac=frac, TF=TF, fracNoise=fracNoise)
	print('generate_data took %d seconds' %( time.time() - startTime))
