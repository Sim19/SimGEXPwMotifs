# coding=utf-8
'''
define class to return values in that class
'''

class Ybetaparam(object):
	__slots__ = ['Y', 'beta', 'parameter']
	
	def __init__(self, Y, beta, parameter):
		self.Y = Y
		self.beta = beta
		self.parameter = parameter
		
	def __getstate__(self):
		return {slot: getattr(self, slot) for slot in self.__slots__}
	
	def __setstate__(self, d):
		for slot in d:
			setattr(self, slot, d[slot])