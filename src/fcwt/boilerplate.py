from .fcwt import Morlet, Scales, FCWT, FCWT_LINSCALES, FCWT_LOGSCALES, FCWT_LINFREQS
import numpy as np
import matplotlib.pyplot as plt

def cwt(input, fs, f0, f1, fn, nthreads=1, scaling="lin", fast=False, norm=True):

    #check if input is array and not matrix
    if input.ndim > 1:
        raise ValueError("Input must be a vector")

    #check if input is single precision and change to single precision if not
    if input.dtype != 'single':
        input = input.astype('single')

    morl = Morlet(2.0) #use Morlet wavelet with a wavelet-parameter of 2.0

    #Generate scales

    if scaling == "lin":
        scales = Scales(morl,FCWT_LINFREQS,fs,f0,f1,fn)
    elif scaling == "log":
        scales = Scales(morl,FCWT_LOGSCALES,fs,f0,f1,fn)
    else:
        scales = Scales(morl,FCWT_LOGSCALES,fs,f0,f1,fn)

    fcwt = FCWT(morl, nthreads, fast, norm)

    output = np.zeros((fn,input.size), dtype='csingle')
    freqs = np.zeros((fn), dtype='single')

    fcwt.cwt(input,scales,output)
    scales.getFrequencies(freqs)

    return freqs, output

def plot(input, fs, f0=0, f1=0, fn=0, nthreads=1, scaling="lin", fast=False, norm=True):
    
    f0 = f0 if f0 else 1/(input.size/fs)
    f1 = f1 if f1 else fs/2
    fn = fn if fn else 100
    freqs, output = cwt(input, fs, f0, f1, fn, nthreads=nthreads, scaling=scaling, fast=fast, norm=norm)

    #create two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    #plot input signal
    ax1.plot(input)
    ax1.set_title('Input signal')

    #plot spectrogram
    ax2.imshow(np.abs(output),aspect='auto')

    #set time on x-axis every 10s and frequency on y-axis for 10 ticks
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Frequency (Hz)')
    ax2.set_title('CWT')
    ax2.set_xticks(np.arange(0,input.size,fs*10),np.arange(0,input.size/fs,10))
    ax2.set_yticks(np.arange(0,fn,fn/10),np.round(freqs[::int(fn/10)]))


    plt.show()
    
