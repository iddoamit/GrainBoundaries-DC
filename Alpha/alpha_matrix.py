import numpy as np
import sympy as sp
from datetime import datetime
import pandas as pd
import os

here = os.path.dirname(os.path.realpath(__file__))+'/'

start = datetime.now()
print('Simulation starting at '+start.time().strftime('%H:%M:%S'))

q = 1.6e-19
e0 = 8.85e-14
er = 11.7
ee = e0*er
k = 8.6e-5
T = 300
kT = k*T
NT = 4e11

VAvec = np.linspace(-0.5, 0.5, 200+1)
ind = np.ones_like(VAvec, dtype=bool)
ind[len(VAvec)//2] = 0


ND = 1e16*np.arange(0.80, 1.21, 0.04)

x1, x2, VA = sp.symbols('x1 x2 VA')

x1vec = np.zeros((len(ND), len(ND), len(VAvec)), dtype=float)
x2vec = np.zeros_like(x1vec)

for left in range(len(ND)):
       for right in range(len(ND)):
              for v in range(len(VAvec)):
                     eqn1 = sp.Eq(q*ND[left]*(x1**2)/(2*ee), q*ND[right]*(x2**2)/(2*ee) + kT*np.log(ND[left]/ND[right]) - VA)
                     eqn2 = sp.Eq(x1*ND[left] + x2*ND[right], NT)
                     eqn3 = sp.Eq(VA, VAvec[v])
                     S = sp.solve((eqn1, eqn2, eqn3), (x1, x2, VA))
                     if S[0][0] >= 0 and S[0][1] >= 0:
                            x1vec[left, right, v] = float(S[0][0])
                            x2vec[left, right, v] = float(S[0][1])
                     else:
                            x1vec[left, right, v] = float(S[1][0])
                            x2vec[left, right, v] = float(S[1][1])
       print('Finished doping level '+str(left+1)+' out of '+str(len(ND)))
       lap = datetime.now()
       print('So far, '+f'{(lap-start).total_seconds():.3f}'+' have passed')

VbiL = np.empty_like(x1vec)
VbiR = np.empty_like(x1vec)
for n in range(len(ND)):
       VbiL[n, :, :] = q*ND[n]*(x1vec[n, :, :]**2)/(2*ee)
       VbiR[:, n, :] = q*ND[n]*(x2vec[:, n, :]**2)/(2*ee)


DiffL = np.zeros_like(VbiL)
DiffR = np.zeros_like(VbiR)
for v in range(len(VAvec)):
       DiffL[:, :, v] = VbiL[:, :, v] - VbiL[:, :, len(VAvec)//2]
       DiffR[:, :, v] = VbiR[:, :, v] - VbiR[:, :, len(VAvec)//2]

alpha = np.zeros_like(DiffL)
for v in range(len(VAvec)):
       alpha[:, :, v] = -DiffL[:, :, v]/VAvec[v]
       

end = datetime.now()
delta = end-start
print('Simulation ended at '+end.time().strftime('%H:%M:%S'))
print('Simulation took '+f'{delta.total_seconds():.3f}'+' seconds')

for n in range(len(ND)):
       data = {'VA': VAvec}
       for n2 in range(len(ND)):
              data.update({f'{ND[n2]:.2e}': alpha[n, n2, :]})
       DF = pd.DataFrame(data=data)
       filename = 'Alpha_'+f'{ND[n]:.2e}'+'.csv'
       DF.to_csv(here+filename, index=None)
print('All done!')