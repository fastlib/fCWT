% fCWT.m Help file for fCWT MEX file
%
% fCWT performs a fast Continuous Wavelet Transformation (CWT) of an input
% signal. 
% 
% Usage:
% [tfm,f] = fCWT(signal, c0, fs, f0, f1, fn, nthreads)
%
% Inputs:
%   signal      - A row vector of singles containing the input signal
%   c0          - The Morlet Wavelet parameter used to tune the time-frequency
%                   tradeoff. A higher value means a higher frequency accuracy 
%                   but lower time accuracy and vica versa.
%   fs          - Signal sample frequency.
%   f0          - Starting frequency.
%   f1          - Ending frequency.
%   fn          - Number of frequencies.
%   nthreads    - Number of threads to use. By default, the number of threads
%
% Output:
%   tfm         - Complex-valued time-frequency matrix. By default, the matrix
%                   is transposed such that rows correspond to the time axis and columns to
%                   the frequency axis.
%   f           - Array containing the frequencies that correspond to the matrix
%                   columns
%
