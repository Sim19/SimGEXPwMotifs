# coding=utf-8
'''
FILETRAIT class - store all information about genereted data
	- create filename given the FILETRAIT-object
	- get FILETRAIT object given filename
'''

#### PACKAGES ####
import os
##################



#### FILETRAIT - CLASS ####
class FILETRAIT(object):
	__slots__ = ['Cond', 'Genes', 'rep', 'Sigma', 'noise', 'method', 'R2method', 'data', 'cwd', 'path', 'frac', 'genMethod', 'estimMethod', 'initial', 'TF', 'fracNoise', 'corr_type']
	
	def __init__(self, Cond, Genes, rep, Sigma, noise, method, R2method, data, cwd, path, frac, genMethod, estimMethod, initial, TF, fracNoise, corr_type):
		self.Cond=Cond
		self.Genes=Genes
		self.rep=rep
		self.Sigma=Sigma
		self.noise=noise
		self.method=method
		self.R2method=R2method
		self.data=data
		self.cwd=cwd
		self.path=path
		self.frac=frac
		self.genMethod=genMethod
		self.estimMethod=estimMethod
		self.initial=initial
		self.TF=TF
		self.fracNoise=fracNoise
		self.corr_type=corr_type
		
	def __getstate__(self):
		return {slot: getattr(self, slot) for slot in self.__slots__}
	
	def __setstate__(self, d):
		for slot in d:
			setattr(self, slot, d[slot])
##############################


#### FUNCTIONS ####
def nameFile(Cond=None, Genes=None, rep=None, Sigma=None, noise=None, method=None, R2method=None, data=None, cwd=None, path=None, frac=None, genMethod=None, estimMethod=None, initial=None, TF=None, fracNoise=None, corr_type=None):
	
	FT = FILETRAIT
	
	if Cond is None:
		COND = '_C_*'
	else:
		COND = '_C_' + str(Cond)
	FT.Cond=COND
	
	if Genes is None:
		GENES = '_G_*'
	else:
		GENES = '_G_' + str(Genes)
	FT.Genes=GENES
	
	if rep is None:
		REP = '_rep_100'
	else:
		REP = '_rep_' + str(rep)
	FT.rep=REP
	
	if Sigma is None:
		SIGMA = ''
	else:
		SIGMA = '_Sigma_' + Sigma
	FT.Sigma=SIGMA
	
	if noise is None:
		NOISE = ''
	else:
		NOISE = '_limixE_' + noise
	FT.noise = NOISE
	
	if method is None:
		METH = ''
	else:
		METH = '_' + method
	FT.method=METH
	
	if R2method is None:
		R2METH = ''
	elif R2method == 'meanVmeanM':
		R2METH = '_R2method_mean_V_mean_M'
	else:
		R2METH = '_R2method_' + R2method
	FT.R2method=R2METH
	
	if data is None:
		DATA = 'beta'
	else:
		DATA = '_' + data
	FT.data=DATA
	
	if cwd is None:
		cwd = os.getcwd()
	FT.cwd = cwd
	
	if path is None:
		path = '/data/simulation/'
	FT.path = path
	
	if frac is None:
		frac = ''
	else:
		frac = '_frac_' + str(int(frac*100))
	FT.frac=frac
	
	if genMethod is None:
		GENMETHOD = '_V_*'
	else:
		GENMETHOD = '_V_' + genMethod
	FT.genMethod=GENMETHOD
	
	if estimMethod is None:
		ESTIMMETHOD = '_limixV_*'
	else:
		ESTIMMETHOD = '_limixV_' + estimMethod
	FT.estimMethod=ESTIMMETHOD
	
	if initial is False or initial is None:
		INITIAL = ''
	else:
		INITIAL = '_initOriginal'
	FT.initial=INITIAL
	
	if TF is None:
		TF = ''
	else:
		TF = '_TF_' + str(TF)
	FT.TF = TF
	
	if fracNoise is not None and 'V' in Sigma:
		FRACNOISE = '_fracS_' + str(int(fracNoise*100))
	else:
		FRACNOISE = ''
	FT.fracNoise = FRACNOISE
	
	if corr_type is not None:
		CORRELATION = '_corr_' + str(corr_type)
	else:
		CORRELATION = ''
	FT.corr_type = CORRELATION
	
	return FT



def getIndex(string, list):
	'''
	given a string, get index of next object in list
		- this should give the parameter value
	'''
	try:
		element = list[ list.index(string) + 1 ]
		
		if string == 'R2method':
			elements = list[ (list.index('R2method') + 1): list.index('frac')]
			el = ''
			for elem in elements:
				el += elem
			element = el
			
		if string == 'V':
			if 'limixV' in list:
				elements = list[ (list.index('V') + 1):list.index('limixV')]
			else:
				elements = list[ (list.index('V') + 1):list.index('rep')]
			if elements.__len__() > 1:
				element = elements[0] + '_' + elements[1]
			else:
				element = elements[0]
			
		if string == 'limixV':
			rep_index = list.index('rep')
			r2meth_index = list.index('R2method')
			limixV_index = list.index('limixV')
			if (limixV_index < rep_index < r2meth_index):
				elements = list[ (limixV_index + 1):rep_index]
			else:
				elements = list[ (limixV_index + 1):r2meth_index]
				
			if elements.__len__() > 1:
				element = elements[0] + '_' + elements[1]
			else:
				element = elements[0]
	
	except:
		element = None
	
	return element
	



def getTraitsFromFile(filename):
	
	# create class FILETRAIT
	FT = FILETRAIT
	
	# split filename into pieces
	file = filename.split('/')[-1].split('.pkl')[0]
	parts = file.split('_')
	
	# put pieces together and save in FILETRAIT
	FT.data = parts[0]
	
	FT.Cond	 = int( getIndex('C', list=parts) )
	FT.Genes = int( getIndex('G', list=parts) )
	FT.Sigma = getIndex('Sigma', list=parts)
	FT.rep	 = int( getIndex('rep', list=parts) )
	FT.frac		= float( getIndex('frac', list=parts))/100
	FT.genMethod	= getIndex('V', list=parts)
	FT.estimMethod	= getIndex('limixV', list=parts)
	FT.noise	= getIndex('limixE', list=parts)
	FT.R2method = getIndex('R2method', list=parts)
	FT.method = None
	FT.cwd = None
	FT.path = None
	if 'initOriginal' in filename:
		FT.initial = True
	else:
		FT.initial = False
		
	tfs = getIndex('TF', list=parts)
	if tfs is None:
		tfs=0
	FT.TF=int(tfs)
	
	if 'V' in FT.Sigma:
		FT.fracNoise=float(getIndex('fracS', list=parts))/100
	else:
		FT.fracNoise=None
	
	if '_corr_' in filename:
		FT.corr_type = getIndex('corr', list=parts)
	else:
		FT.corr_type='pearson'
		
		
	return FT



def makeStringFileTrait(FILETRAIT):
	'''
	turn FILETRAIT-object into FILETRAIT object with strings for
	file-namig
	:param FILETRAIT:
	:return: FILETRAIT-obejct filled with naming strings
	'''
	FT = nameFile(Cond=FILETRAIT.Cond, Genes=FILETRAIT.Genes, rep=FILETRAIT.rep, Sigma=FILETRAIT.Sigma, noise=FILETRAIT.noise, method=FILETRAIT.method, R2method=FILETRAIT.R2method, data=FILETRAIT.data, cwd=FILETRAIT.cwd, path=FILETRAIT.path, frac=FILETRAIT.frac, genMethod=FILETRAIT.genMethod, estimMethod=FILETRAIT.estimMethod, initial=FILETRAIT.initial, TF=FILETRAIT.TF, fracNoise=FILETRAIT.fracNoise, corr_type=FILETRAIT.corr_type)
	
	return FT
##############################