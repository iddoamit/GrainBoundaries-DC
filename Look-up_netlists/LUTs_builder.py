import numpy as np
import pandas as pd
import sympy as sp
from datetime import datetime

start = datetime.now()

path = '/LUTS/LUTsDep/'

q = 1.6e-19
e0 = 8.85e-14
er = 11.7
ee = e0*er
k = 8.6e-5
T = 300
kT = k*T
meff = 0.86
NC = 2.5e19*(meff**1.5)
A = 1
AA = 120*meff

def Vbi_calculator(ND, NT):
	x1, x2 = sp.symbols('x1 x2')
	eqn1 = sp.Eq(q*ND[0]*(x1**2)/(2*ee), q*ND[1]*(x2**2)/(2*ee)+kT*np.log(ND[0]/ND[1]))
	eqn2 = sp.Eq(x1*ND[0]+x2*ND[1], NT)
	S = sp.solve((eqn1, eqn2), (x1, x2), dict=True)
	if S[0][x1] > 0 and S[0][x2] > 0:
		x = [float(S[0][x1]), float(S[0][x2])]
	else:
		x = [float(S[1][x1]), float(S[1][x2])]
	Vbi = [q*ND[0]*(x[0]**2)/(2*ee), q*ND[1]*(x[1]**2)/(2*ee)]
	return Vbi

def depletion_calculator(ND, NT, Vbi, length=200):
	Vmax = 0.75*(Vbi[0]+Vbi[1])
	VA_vec = np.linspace(-Vmax, Vmax, length)
	x1, x2, VA = sp.symbols('x1 x2 VA')
	eqn1 = sp.Eq(q*ND[0]*(x1**2)/(2*ee), q*ND[1]*(x2**2)/(2*ee)+kT*np.log(ND[0]/ND[1])-VA)
	eqn2 = sp.Eq(x1*ND[0]+x2*ND[1], NT)
	Xvalues = np.zeros((length, 2), dtype=float)
	for k in range(length):
		eqn3 = sp.Eq(VA, VA_vec[k])
		S = sp.solve((eqn1, eqn2, eqn3), (x1, x2, VA), dict=True)
		if S[0][x1] > 0 and S[0][x2] > 0:
			Xvalues[k, :] = [float(S[0][x1]), float(S[0][x2])]
		else:
			Xvalues[k, :] = [float(S[1][x1]), float(S[1][x2])]
	return VA_vec, Xvalues
	

def flux_calculator(ND, Xvalues):
	phi = np.zeros_like(Xvalues)
	for k in range(len(Xvalues)):
		phi[k, 0] = q*ND[0]*(Xvalues[k, 0]**2)/(2*ee) + kT*np.log(NC/ND[0])
		phi[k, 1] = q*ND[1]*(Xvalues[k, 1]**2)/(2*ee) + kT*np.log(NC/ND[1])
	flux = np.exp(-phi[:, 0]/kT) - np.exp(-phi[:, 1]/kT)
	return flux

def resistivity_calculator(VA_vec, flux):
	res = VA_vec/(flux*A*AA*(T**2))
	return res

def LUT_writer(ND, VA_vec, res, Xvalues):
	filename = 'LUT_'+f'{ND[0]/1e13:03.0f}'+'_'f'{ND[1]/1e13:03.0f}'+'_1e13.csv'
	data = {'Bias': VA_vec, 'Resistance': res, 'X_Left': Xvalues[:, 0], 'X_Right': Xvalues[:, 1]}
	DF = pd.DataFrame(data=data)
	DF.to_csv(path+filename, index=None)


ND_vec = np.arange(86, 114+1, 0.5)*1e14
NT = 4e11

num_files = (len(ND_vec)**2)
for i in range(len(ND_vec)):
	for j in range(len(ND_vec)):
		Vbi = Vbi_calculator([ND_vec[i], ND_vec[j]], NT)
		VA_vec, xvalues = depletion_calculator([ND_vec[i], ND_vec[j]], NT, Vbi)
		flux = flux_calculator([ND_vec[i], ND_vec[j]], xvalues)
		res = resistivity_calculator(VA_vec, flux)
		LUT_writer([ND_vec[i], ND_vec[j]], VA_vec, res, xvalues)
		print('Completed file '+str((i)*len(ND_vec)+j+1)+' ouf of '+str(num_files))

end = datetime.now()
delta = end-start
print('Simulation ended at '+end.time().strftime('%H:%M:%S'))
print('Simulation took '+f'{delta.total_seconds():.3f}'+' seconds')
print('All done!')