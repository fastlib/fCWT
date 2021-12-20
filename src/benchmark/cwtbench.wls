#!/usr/bin/env wolframscript

n = 100000;
fs = 64;
signal = Table[Sin[2*3.14159265*((1+(7*i)/n)*(i/fs))],{i,n}]

runs = 100;
For[k = 1, k <= runs, k++,
    {t,b} = AbsoluteTiming[ContinuousWaveletTransform[signal,MorletWavelet[],{6,500}]];
    Pause[10];
    Print[t]
]