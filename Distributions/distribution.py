import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import pandas as pd
import os
from datetime import datetime

start = datetime.now()

here = os.path.dirname(os.path.realpath(__file__))+'/'

dimension = 3
ND = 1e16
L = 1e-3
volume = L**dimension
radius = 4.5e-5
num_doping = np.round((ND**(dimension/3))*volume)

nd_mat = np.random.uniform(0, L, (int(num_doping), dimension))

boundary = 5e-5
num_crystallites = 20**dimension
pos = np.random.uniform(0+boundary, L-boundary, (num_crystallites, dimension))

nd_incl = np.zeros(pos.shape[0])
for k in range(len(nd_incl)):
	dist = np.sqrt(np.sum((nd_mat-pos[k, :])**2, 1))
	nd_incl[k] = dist[dist<=radius].shape[0]
	if k%80 == 0:
		print('Finished '+str(k//80)+'%')

data = {'4.5em5': nd_incl}
DF = pd.DataFrame(data=data)
DF.to_csv(here+'ND_distribution_4_5em5.csv', index=None)

end = datetime.now()
delta = end-start
print('Simulation ended at '+end.time().strftime('%H:%M:%S'))
print('Simulation took '+f'{delta.total_seconds():.3f}'+' seconds')
print('All done!')