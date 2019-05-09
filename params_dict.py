# coding=utf-8
# define class to return values in that class
class params_dict(object):
	__slots__ = ['Cond', 'Genes', 'TF', 'V', 'Sigma', 'motif', 'delta', 'sigma_beta', 'I_G', 'mu_C', 'mu_G', 'frac', 'TF', 'fracNoise']
	
	def __init__(self, C, G, T, V, Sigma, motif, delta, sigma_beta, I_G, mu_G, mu_C, frac, TF, fracNoise):
		self.Cond = C
		self.Genes = G
		self.TF = T
		self.V = V
		self.Sigma = Sigma
		self.motif = motif
		self.delta = delta
		self.sigma_beta = sigma_beta
		self.I_G = I_G
		self.mu_G = mu_G
		self.mu_C = mu_C
		self.frac = frac
		self.TF = TF
		self.fracNoise = fracNoise
	
	def __getstate__(self):
		return {slot: getattr(self, slot) for slot in self.__slots__}
	
	def __setstate__(self, d):
		for slot in d:
			setattr(self, slot, d[slot])