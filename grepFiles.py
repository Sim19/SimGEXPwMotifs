# coding=utf-8
'''
get and create Filename from information on data
'''

#### PACKAGES ####
import os, sys
cwd = os.getcwd()
sys.path.append(cwd + '/code/funcs/')
import namingFile as nf
##################


#### FUNCTIONS ####
def grepFiles(rep, Cond=None, Genes=None, Sigma=None, noise=None, method=None, cwd=None, path=None, data='beta', R2method=None, frac=None, genMethod=None, estimMethod=None, initial=False, TF=None, fracNoise=None, corr_type=None):
	'''
	Function to get pattern of filename given details.
	:param rep: 	#repitions
	:param Cond: 	#conditions/col-dimension of Y
	:param Genes: 	#gens/peaks/row-dimension of Y
	:param Sigma: 	#noise of generated data
	:param noise: 	#noise of fitted data
	:param method: 	#method - mostly None
	:param cwd: 	#cwd - current working directory
	:param data: 	#data - type of data - driving filename
	:param path: 	#path for data storage
	:param R2method: #method used to create Signal-Noise-ratio
	:param frac: 	#ratio of Signal-noise of generated data
	:param genMethod: #type of correlation used for data generation
	:param estimMethod: #type of correlation used for data fit
	:param initial:	bool if initialing values used as initialization
	:return: string-filename
	'''
	NAMES=nf.nameFile(Cond=Cond, Genes=Genes, rep=rep, Sigma=Sigma, noise=noise, method=method, cwd=cwd, path=path, data=data, R2method=R2method, frac=frac, genMethod=genMethod, estimMethod=estimMethod, initial=initial, TF=TF, fracNoise=fracNoise, corr_type=corr_type)
	
	filenames = grepFilesFromNames(NAMES)
	
	return filenames


def grepFilesFromNames(NAMES):
	METH = NAMES.method
	data = NAMES.data
	if data[0] == '_':
		data = data[1::]
	COND=NAMES.Cond
	GENES=NAMES.Genes
	REP=NAMES.rep
	SIGMA=NAMES.Sigma
	NOISE=NAMES.noise
	R2METHOD=NAMES.R2method
	FRAC=NAMES.frac
	GENMETH=NAMES.genMethod
	ESTIMMETH=NAMES.estimMethod
	CWD=NAMES.cwd
	PATH=NAMES.path
	INITIAL=NAMES.initial
	TF=NAMES.TF
	FRACNOISE=NAMES.fracNoise
	CORRELATION=NAMES.corr_type
	
	# saving data correctly
	if ESTIMMETH == '_limixV_*':
		ESTIMMETH_save=''
	else:
		ESTIMMETH_save=ESTIMMETH
	
	if CWD is None:
		cwd = os.getcwd()
	else:
		cwd = CWD

	if PATH is None:
		path = '/data/simulation/'
	else:
		path = PATH

	if FRAC == '':
		FRAC_save = '_frac_*'
	else:
		FRAC_save = FRAC
	
	if SIGMA == '':
		SIGMA_save = '_Sigma_*'
	else:
		SIGMA_save = SIGMA
	
	
	if data == 'beta':
		
		FILES = 'beta_res' + METH + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA + FRACNOISE + NOISE + FRACNOISE + R2METHOD + FRAC + INITIAL + '.pkl'
	
	elif data == 'beta_both':
		
		FILES = 'beta_*' + METH + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA + FRACNOISE + NOISE + R2METHOD + FRAC + INITIAL + '.pkl'

	elif data == 'beta_posterior':
		
		FILES = 'beta_posterior' + METH + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA + FRACNOISE + NOISE + R2METHOD + FRAC + INITIAL + '.pkl'
	
	elif data == 'beta_ridge':
		
		FILES = 'beta_ridge' + METH + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA + FRACNOISE + NOISE + R2METHOD + FRAC + INITIAL + '.pkl'
		
	elif data == 'abs_beta':
		
		FILES = 'abs_mean' + '_beta' + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA + FRACNOISE + NOISE + METH + R2METHOD + FRAC + INITIAL + '.pkl'
	
	
	elif data == 'Ybetaparam':
		FILES = 'Ybetaparams_generated' + COND + GENES + TF + GENMETH + REP + SIGMA + FRACNOISE + R2METHOD + FRAC + '.pkl'
	
		
	elif data == 'Rsquared':
		
		FILES = 'Rsquared' + COND + GENES + TF + REP + SIGMA_save + FRACNOISE + NOISE + GENMETH + ESTIMMETH_save + R2METHOD + FRAC_save + INITIAL + '.pkl'
	
	
	elif data == 'VSigmaKiY':
		
		FILES = 'VSigmaKiY' + COND + GENES + TF + REP + SIGMA_save + FRACNOISE + NOISE + GENMETH + ESTIMMETH + R2METHOD + FRAC_save + INITIAL +  '.pkl'
	
	
	elif data == 'verbose':
		
		FILES = 'verbose' + COND + GENES + TF + REP + SIGMA_save + FRACNOISE + NOISE + GENMETH + ESTIMMETH + R2METHOD + FRAC + INITIAL + '.pkl'
		
	
	elif data == 'KLdivV':
		
		FILES = 'KLdivV' + COND + GENES + TF + REP + SIGMA_save + FRACNOISE + NOISE + GENMETH + ESTIMMETH_save + R2METHOD + FRAC_save + INITIAL + '.pkl'
	
	
	elif data == 'KLdivSigma':
		
		FILES = 'KLdivSigma' + COND + GENES + TF + REP + SIGMA_save + FRACNOISE + NOISE + GENMETH + ESTIMMETH_save + R2METHOD + FRAC_save + INITIAL + '.pkl'
	
	
	elif data == 'corr':
		
		FILES = 'Correlation' + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA_save + FRACNOISE + NOISE + R2METHOD + FRAC + INITIAL + '.pkl'
		FILES = 'Correlation' + '_beta_both' + METH + COND + GENES + TF + GENMETH + ESTIMMETH_save + REP + SIGMA_save + FRACNOISE + NOISE + R2METHOD + FRAC + INITIAL + CORRELATION + '.pkl'
		
	
	
	else:
		
		raise ValueError('No such datatype implemented as %s' %data)
	
	FILENAMES = cwd + path + FILES
	
	return FILENAMES

##################
