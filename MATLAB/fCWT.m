% fCWT.m Help file for fCWT MEX file
%
% fCWT performs a fast Continuous Wavelet Transformation (CWT) of an input
% signal. 
% 
% Usage:
% [tfm,f] = fCWT(signal, c0, fs, s0, s1, nds)
%
% Inputs:
%   signal - A row vector of singles containing the input signal
%   c0     - The Morlet Wavelet parameter used to tune the time-frequency
%            tradeoff. A higher value means a higher frequency accuracy 
%            but lower time accuracy and vica versa.
%   fs     - Signal sample frequency
%   s0     - Starting octave of the scale range. Scales are defined
%            dyadically (i.e., 2^(octave+nds/ds)). Consequently, octaves are defined
%            as integer powers of 2. A starting octave of 2 means scales start at
%            s=2^2=4. 
%   s1     - Ending octave of the scale range.
%   nds    - Number of octave subdivisions.
%
% Output:
%   tfm    - Complex-valued time-frequency matrix. By default, the matrix
%            is transposed such that rows correspond to the time axis and columns to
%            the frequency axis.
%   f      - Array containing the frequencies that correspond to the matrix
%            columns
%
