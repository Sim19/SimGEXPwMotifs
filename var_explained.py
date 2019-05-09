# coding=utf-8
'''
Functions to compute R2 value in linear regression
of
$$ Y = M \beta + noise $$
'''

#### PACKAGES ####
import numpy as np
##################


#### FUNCTIONS #####
def var_explained(Y, A, b):
	
	var = np.cov(Y, np.dot(A, b).reshape(-1, order='F'))[0,1]
	
	return var


def varCoeff(Y, A, b):
	
	res = Y - np.dot(A, b).reshape(-1, order='F')
	total = Y - Y.mean()
	
	SS_res = np.dot(res.T, res)
	SS_tot = np.dot(total.T, total)
	
	Rsquared = 1 - SS_res/SS_tot
	
	return Rsquared


def varCoeff2(Y, Y_hat):
	
	res = (Y - Y_hat).reshape(-1, order='F')
	total = (Y - Y.mean()).reshape(-1, order='F')
	
	SS_res = np.dot(res.T, res)
	SS_tot = np.dot(total.T, total)
	
	Rsquared = 1 - SS_res/SS_tot
	
	return Rsquared


def correctParam(param_comp, value, param_estim):

	diff = value/param_comp
	
	if diff < 0.8 or diff > 1.2:
		
		param_update = diff * param_estim
		
	return param_update^2
##################



