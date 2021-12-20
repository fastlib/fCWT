import numpy as np
import pywt
from scipy import signal

import timeit
import time

print('Generate data for CWT benchmark...')

#Generate data
fs = 64
n = 100000
x = np.arange(n)
y = np.sin(2*np.pi*((1+(7*x)/n)*(x/fs)))

#runs
runs = 100;


#======== BENCHMARK PYWAVELET ============#
print('Start PyWavelet CWT benchmark...')

#Initialize timing array
pyw_times = np.array([]);

#Perform CWT via PyWavelet
for i in range(runs):
    st = timeit.default_timer()
    coef, freqs=pywt.cwt(y,np.logspace(1,6,500*6, base=2),'cmor1.5-1.0',method='fft')
    pyw_times = np.append(pyw_times,[timeit.default_timer() - st]);
    print(pyw_times[i]);
    time.sleep(10);

#Output
print('PyWavelet time M: ', np.mean(pyw_times)) 
print('PyWavelet time SD: ', np.std(pyw_times)) 


#======== BENCHMARK SCIPY ============#
print('Start SciPy CWT benchmark...')

#Initialize timing array
scp_times = np.array([]);

#Perform CWT via SciPy
for i in range(runs):
    st = timeit.default_timer()
    cwtmatr = signal.cwt(y, signal.morlet2, np.logspace(1,6,500*6, base=2),w=np.pi*2)
    scp_times = np.append(scp_times,[timeit.default_timer() - st]);
    print(scp_times[i]);
    time.sleep(10);

#Output
print('SciPy time M: ', np.mean(scp_times)) 
print('SciPy time SD: ', np.std(scp_times)) 