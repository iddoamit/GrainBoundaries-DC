import numpy as np
import matplotlib.pyplot as plt
import math

def frequency_gen(data):
	new_data = np.copy(data)
	new_data.sort()
	freq = np.arange(1, len(new_data)+1)
	freq = (freq-0.3)/(len(new_data)+0.4)
	return new_data, freq

def Rsquare(y, ycalc):
	num = np.sqrt(np.sum((y - ycalc)**2))
	den = np.sqrt(np.sum((y - np.mean(y))**2))
	return 1-num/den

def search_delta(data, freq, res=100):
	delta_vec = np.linspace(0, data[0], res+1)
	delta_vec = delta_vec[:-1]
	R2 = np.zeros_like(delta_vec)
	for k in range(len(delta_vec)):
		P = np.polyfit(np.log(data-delta_vec[k]), np.log(-np.log(1-freq)), 1)
		y = np.log(-np.log(1-freq))
		ycalc = np.polyval(P, np.log(data-delta_vec[k]))
		R2[k] = Rsquare(y, ycalc)
	delta = delta_vec[R2==max(R2)]
#	print('Delta is '+str(delta))
	return delta

def Weibull_parameters(data, freq, delta):
	P = np.polyfit(np.log(data-delta), np.log(-np.log(1-freq)), 1)
	beta = P[0]
	theta = np.exp(-P[1]/P[0])
#	print('beta is '+str(beta))
#	print('theta is '+str(theta))
	return beta, theta

def CDFPDF(data, beta, theta, delta):
	x = np.linspace(data[0], data[-1], 1000)
	y = 1 - np.exp(-(((x-delta)/theta)**beta))
	z = np.gradient(y, x)
	return x, y, z

def WeibMean(beta, theta, delta):
	mean = delta + theta*math.gamma(1+(1/beta))
	return mean

def WeibMedian(beta, theta, delta):
	median = delta + theta*((np.log(2))**(1/beta))
	return median

def WeibMode(beta, theta, delta):
	if beta > 1:
		mode = delta + theta*(((beta-1)/beta)**(1/beta))
	else:
		mode = 0
	return mode

def WeibVar(beta, theta, delta):
	A = math.gamma(1+(2/beta))
	B = (math.gamma(1+(1/beta)))**2
	var = (theta**2)*(A-B)
	return var

def running_Weibull(data):
	data_use = np.copy(data)
	S = data_use.shape
	running_mean = np.zeros(S[0])
	running_std = np.zeros(S[0])
	Weib_params = np.zeros((S[0], 3))

	for k in range(S[0]):
		sub_data = data_use[k, :]
		sub_data, freq = frequency_gen(sub_data)
		delta = search_delta(sub_data, freq)
		beta, theta = Weibull_parameters(sub_data, freq, delta)
		Weib_params[k, 0] = beta
		Weib_params[k, 1] = theta
		Weib_params[k, 2] = delta
		running_mean[k] = WeibMean(beta, theta, delta)
		running_std[k] = np.sqrt(WeibVar(beta, theta, delta))
	return running_mean, running_std, Weib_params
