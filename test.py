import fcwt
import numpy as np
import matplotlib.pyplot as plt

#Initialize
fs = 1000
n = fs*100 #100 seconds
ts = np.arange(n)

#Generate linear chirp
signal = np.sin(2*np.pi*((0.1+(2*ts)/n)*(ts/fs)))

f0 = 0.1 #lowest frequency
f1 = 5 #highest frequency
fn = 200 #number of frequencies

#Calculate CWT without plotting...
freqs, out = fcwt.cwt(signal, fs, f0, f1, fn)

#... or calculate and plot CWT
fcwt.plot(signal, fs, f0=f0, f1=f1, fn=fn)