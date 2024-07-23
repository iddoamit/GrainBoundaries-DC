import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

here = os.path.dirname(os.path.realpath(__file__))+'/'
datapath = here+'Netlists/'

sys.path.append(here)
import Weibull as Weibull

grains = np.arange(2, 10+1)


A = 1
radius = 3e-5

BETA = np.zeros((10001, len(grains)))
MEAN = np.zeros((10001, len(grains)))
STD = np.zeros((10001, len(grains)))

for g in range(len(grains)):
	print('Processing circuit '+str(grains[g]))
	filename = 'results'+str(grains[g])+'.csv'
	data = pd.read_csv(datapath+filename).to_numpy()
	bias = data[:, 1]
	curr = data[:, 2:]
	resistance = np.zeros_like(curr)
	for i in range(resistance.shape[1]):
		resistance[:, i] = bias/curr[:, i]
	L = grains[g]*2*radius
	resistivity = resistance*A/L
	running_mean, running_std, Weib_params = Weibull.running_Weibull(resistivity)
	BETA[:, g] = Weib_params[:, 0]
	MEAN[:, g] = running_mean
	STD[:, g] = running_std
	print('Processed circuit '+str(grains[g]))
