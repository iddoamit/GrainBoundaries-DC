import numpy as np
import pandas as pd
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
mu = 1400
L = 1e-4
radius = 3e-5
volume = 4/3*np.pi*(radius**3)

def check_profile_boundaries(profile):
	check = np.zeros_like(profile, dtype=bool)
	for k in range(len(profile)):
		if profile[k] >= 860 and profile[k] <= 1145:
			check[k] = True
	double_check = sum(check) == len(profile)
	return ~double_check

def create_profile(grains=6, mu=1130.97/volume, sigma=33.713/volume):
	repeat = True
	while repeat:
		doping = np.random.normal(mu, sigma, grains)
		profile = []
		for k in range(len(doping)):
			profile.append(int(doping[k]//5e13*5))
		repeat = check_profile_boundaries(profile)

	return profile

def LUT_formatter(sub_profile):
	filename = 'LUT_'+f'{sub_profile[0]:02.0f}'+'_'+f'{sub_profile[1]:02.0f}'+'_1e13.csv'
	data = pd.read_csv(path+filename, header=0).to_numpy()
	string = ''
	for k in range(len(data)):
		string = string+f'{data[k][0]:.4f}'+','
		string = string+f'{data[k][1]:.4f}'+','
	string = string[:-1]
	string = string+')'
	return string


def grain_resistance(nominal_doping):
	ND = float(nominal_doping*1e13)
	sigma = q*mu*ND
	return radius/(sigma*A)

def grain_resistance_LUT(sub_profile, side):
	Side_flag = 0
	if side == 'Right':
		Side_flag = 1
	filename = 'LUT_'+f'{sub_profile[0]:02.0f}'+'_'+f'{sub_profile[1]:02.0f}'+'_1e13.csv'
	data = pd.read_csv(path+filename, header=0).to_numpy()
	ND = float(sub_profile[Side_flag]*1e13)
	sigma = q*mu*ND
	string = ''
	for k in range(len(data)):
		string = string+f'{data[k][0]:.4f}'+','
		string = string+f'{data[k][Side_flag+2]:.4e}'+','
	string = string[:-1]
	string = string+')'
	return string

def netlist_writer(profile, identifier):
	filename = 'netlist_LUT6_'+f'{identifier:04.0f}'+'.net'
	ff = open(path+'Netlists/'+filename, 'w')
	ff.write('* Simulation '+str(identifier)+', with grain profiles ')
	for k in range(len(profile)):
		ff.write(f'{profile[k]*1e13:.2e}')
		ff.write(', ')
	ff.write('\n')

	ff.write('V1 n001 0 1\n')

	for k in range(len(profile)):
		if k == 0:
			res1 = grain_resistance(profile[0])
			res2 = grain_resistance_LUT([profile[0], profile[1]], 'Left')
			ff.write('Rg0L'+' '+'n001 n002'+' ')
			ff.write(f'{res1:.4e}'+'\n')
			ff.write('Rg0R'+' '+'n002 n003'+' ')
			ff.write('R=table(V(n003,n004),')
			ff.write(res2)
			ff.write('\n')
		elif k == len(profile)-1:
			res1 = grain_resistance_LUT([profile[-2], profile[-1]], 'Right')
			res2 = grain_resistance(profile[-1])
			ff.write('Rg'+str(len(profile))+'L'+' ')
			ff.write('n'+f'{3*(k+1)-2:03.0f}'+' '+'n'+f'{3*(k+1)-1:03.0f}'+' ')
			ff.write('R=table(V(n'+f'{3*(k+1)-3:03.0f}'+',n'+f'{3*(k+1)-2:03.0f}'+'),')
			ff.write(res1)
			ff.write('\n')
			ff.write('Rg'+str(len(profile))+'R'+' ')
			ff.write('n'+f'{3*(k+1)-1:03.0f}'+' '+'0'+' ')
			ff.write(f'{res2:.4e}'+'\n')
		else:
			res1 = grain_resistance_LUT([profile[k-1], profile[k]], 'Right')
			res2 = grain_resistance_LUT([profile[k], profile[k+1]], 'Left')
			ff.write('Rg'+str(k+1)+'L'+' ')
			ff.write('n'+f'{3*(k+1)-2:03.0f}'+' '+'n'+f'{3*(k+1)-1:03.0f}'+' ')
			ff.write('R=table(V(n'+f'{3*(k+1)-3:03.0f}'+',n'+f'{3*(k+1)-2:03.0f}'+'),')
			ff.write(res1)
			ff.write('\n')
			ff.write('Rg'+str(k+1)+'R'+' ')
			ff.write('n'+f'{3*(k+1)-1:03.0f}'+' '+'n'+f'{3*(k+1):03.0f}'+' ')
			ff.write('R=table(V(n'+f'{3*(k+1):03.0f}'+',n'+f'{3*(k+1)+1:03.0f}'+'),')
			ff.write(res2)
			ff.write('\n')
		
	for k in range(len(profile)-1):
		ff.write('Rb'+str(k+1)+' ')
		ff.write('n'+f'{3*(k+1):03.0f}'+' ')
		ff.write('n'+f'{3*(k+1)+1:03.0f}'+' ')
		ff.write('R=table(V(n'+f'{3*(k+1):03.0f}'+',n'+f'{3*(k+1)+1:03.0f}'+'),')
		ff.write(LUT_formatter([profile[k], profile[k+1]]))
		ff.write('\n')

	ff.write('.dc V1 -0.5 0.5 0.0001'+'\n')
	ff.write('.backanno'+'\n')
	ff.write('.end'+'\n')
	ff.close()

for i in range(1000):
	profile = create_profile()
	netlist_writer(profile, i)

end = datetime.now()
delta = end-start
print('Simulation ended at '+end.time().strftime('%H:%M:%S'))
print('Simulation took '+f'{delta.total_seconds():.3f}'+' seconds')
print('All done!')

