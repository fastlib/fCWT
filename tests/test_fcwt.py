from fcwt import *
import os
import numpy as np
import matplotlib.pyplot as plt
#import timeit

def cwt2(x,fs,f0,f1,fn):
    
    #Initialize
    morl = Morlet(2.0)

    out = np.zeros((fn,len(x)), dtype='csingle')
    freqs = np.zeros((fn), dtype='single')
    scales = Scales(morl,FCWT_LOGSCALES,fs,f0,f1,fn)
    scales.getFrequencies(freqs)
    fcwt = FCWT(morl, 8, True, True)
    fcwt.cwt(x, scales, out)

    return out, freqs

def test_odd():
    #Generate data
    fs = 1000
    n = 1000 
    x = np.arange(n)
    hz = 5
    input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
    f0 = 1
    f1 = 20
    fn = 1000

    morl = Morlet(16.0)

    out = np.zeros((fn,len(x)), dtype='csingle')
    freqs = np.zeros((fn), dtype='single')
    scales = Scales(morl,FCWT_LINFREQS,fs,f0,f1,fn)
    scales.getFrequencies(freqs)
    fcwt = FCWT(morl, 1, True, True)
    fcwt.cwt(input, scales, out)

    max = np.argmax(np.mean(np.abs(out),axis=1))
    outhz = freqs[max]
    maxes = np.mean(np.abs(out),axis=1)

    for i,max in enumerate(maxes):
        print(freqs[i],max)
    
    assert np.round(hz*2)/2 == np.round(outhz*2)/2

def test_shape():
    #Generate data
    fs = 100
    n = 2000000
    x = np.arange(n)
    input = np.sin(2*np.pi*((1+(7*x)/n)*(x/fs)), dtype='single')
    
    f0 = 1
    f1 = 32
    fn = 30

    out, scs = cwt2(input,fs,f0,f1,fn)
    print(np.max(np.abs(out)))
    assert out.shape[0]==fn and out.shape[1]==n


def test_complex_input():
    #Generate data
    fs = 100
    n = 2000000
    x = np.arange(n)

    #generate complex sine wave
    input = np.sin(2*np.pi*((1+(7*x)/n)*(x/fs)), dtype='single')
    input = input + 1j*input
    
    f0 = 1
    f1 = 32
    fn = 30

    morl = Morlet(2.0)

    out = np.zeros((fn,len(x)), dtype='csingle')
    freqs = np.zeros((fn), dtype='single')
    scales = Scales(morl,FCWT_LOGSCALES,fs,f0,f1,fn)
    fcwt = FCWT(morl, 8, True, True)
    fcwt.ccwt(input, scales, out)

    assert out.shape[0]==fn and out.shape[1]==n

def test_wavelet():
    #Generate data
    morl = Morlet(2.0)
    sup = morl.getSupport(2.0)

    wavelet = np.zeros(sup*2+1,dtype='csingle')
    morl.getWavelet(2.0,wavelet)
    
    assert np.abs(np.sum(wavelet)) != 0

def test_morlet():
    #Generate data
    fs = 1000
    n = 10000
    x = np.arange(n)
    hz = 5
    input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
    f0 = 1
    f1 = 20
    fn = 1000

    morl = Morlet(16.0)

    out = np.zeros((fn,len(x)), dtype='csingle')
    freqs = np.zeros((fn), dtype='single')
    scales = Scales(morl,FCWT_LINFREQS,fs,f0,f1,fn)
    scales.getFrequencies(freqs)
    fcwt = FCWT(morl, 8, True, True)
    fcwt.cwt(input, scales, out)

    max = np.argmax(np.mean(np.abs(out),axis=1))
    outhz = freqs[max]
    print(np.mean(np.abs(out),axis=1))
    maxes = np.mean(np.abs(out),axis=1)
    for i,max in enumerate(maxes):
        print(freqs[i],max)
    
    assert np.round(hz*10)/10 == np.round(outhz*10)/10


def test_quick():
    #Generate data
    fs = 1000
    n = 10000
    x = np.arange(n)
    hz = 5
    input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
    f0 = 1
    f1 = 20
    fn = 1000

    freqs, out = cwt(input, fs, f0, f1, fn)

    max = np.argmax(np.mean(np.abs(out),axis=1))
    outhz = freqs[max]
    print(np.mean(np.abs(out),axis=1))
    maxes = np.mean(np.abs(out),axis=1)
    for i,max in enumerate(maxes):
        print(freqs[i],max)
    
    assert np.round(hz*10)/10 == np.round(outhz*10)/10


# def test_dog():
#     #Generate data
#     fs = 64
#     n = 10000
#     x = np.arange(n)
#     hz = 1+np.random.rand(1)*10
#     input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
#     f0 = 1
#     f1 = 32
#     fn = 300

#     morl = DOG(2)

#     out = np.zeros((fn,len(x)), dtype='csingle')
#     freqs = np.zeros((fn), dtype='single')
#     scales = Scales(morl,FCWT_LINFREQS,fs,f0,f1,fn)
#     scales.getFrequencies(freqs)
#     fcwt = FCWT(morl, 8, True)
#     fcwt.cwt(input, scales, out)

#     max = np.argmax(np.mean(np.abs(out),axis=1))
#     outhz = freqs[max]
#     print(freqs)
#     print(hz,outhz)
    
#     assert np.round(hz*10)/10 == np.round(outhz*10)/10


# def test_paul():
#     #Generate data
#     fs = 64
#     n = 10000
#     x = np.arange(n)
#     hz = 1+np.random.rand(1)*10
#     input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
#     f0 = 1
#     f1 = 32
#     fn = 300

#     morl = Paul(4)

#     out = np.zeros((fn,len(x)), dtype='csingle')
#     freqs = np.zeros((fn), dtype='single')
#     scales = Scales(morl,FCWT_LINFREQS,fs,f0,f1,fn)
#     scales.getFrequencies(freqs)
#     fcwt = FCWT(morl, 8, True)
#     fcwt.cwt(input, scales, out)

#     max = np.argmax(np.mean(np.abs(out),axis=1))
#     outhz = freqs[max]
#     print(freqs)
#     print(hz,outhz)
    
#     assert np.round(hz*10)/10 == np.round(outhz*10)/10


def test_lowfrequency():
    #Generate data
    fs = 100
    n = 60*fs
    x = np.arange(n)
    hz = 1+np.random.rand(1)*31
    input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
    f0 = 0.01
    f1 = 10
    fn = 300

    out, freqs = cwt2(input,fs,f0,f1,fn)

def test_plan():
    #Generate data
    fs = 64
    n = 10000
    x = np.arange(n)
    hz = 1+np.random.rand(1)*10
    input = np.sin(2*np.pi*hz*(x/fs), dtype='single')
    
    f0 = 1
    f1 = 32
    fn = 300

    morl = Morlet(2.0)
    fcwt = FCWT(morl, 1, True, False)
    fcwt.create_FFT_optimization_plan(2048,0)
    fname = "n2048_t1.wis"

    assert os.path.isfile(fname) 

    fcwt = FCWT(morl, 8, True, False)
    fcwt.create_FFT_optimization_plan(2048,0)
    fname = "n2048_t8.wis"

    assert os.path.isfile(fname) 