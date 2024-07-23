import sys
sys.path.append('/Users/xxxx/.local/pipx/venvs/ltspice/lib/python3.12/site-packages/')
import ltspice
import os
import numpy as np
import time
import subprocess
import pandas as pd
from datetime import datetime

start = datetime.now()

path = '/LUTS/LUTsDep/Netlists/'
files = os.listdir(path)

i = 0
for file in files:
	subprocess.run(['/Applications/LTspice.app/Contents/MacOS/LTspice', '-b', path+file])
	#os.popen('/Applications/LTspice.app/Contents/MacOS/LTspice -b '+path+file)
	#time.sleep(5)
	try:
		LT = ltspice.Ltspice(path+file[:-4]+'.raw')
		LT.parse()
		if i == 0:
			i_data = LT.get_data('V(n001)')
		i_data = np.vstack((i_data, -LT.get_data('I(V1)')))
		i += 1
		print('Completed iteration '+str(i))
	except:
		print('Processing file '+file+' did not succeed.')
	


DF = pd.DataFrame(data=i_data.T)
DF.to_csv(path+'results.csv')

end = datetime.now()
delta = end-start
print('Simulation ended at '+end.time().strftime('%H:%M:%S'))
print('Simulation took '+f'{delta.total_seconds():.3f}'+' seconds')
print('All done!')