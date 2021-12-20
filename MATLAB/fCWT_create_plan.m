% fCWT_create_plan.m Help file for fCWT MEX file
%
% fCWT performs a fast Continuous Wavelet Transformation (CWT) of an input
% signal. This function creates an optimization plan that optimizes fCWT for
% all signallengths up to N on hardware setup X. 
% 
% Usage:
% fCWT_create_plan(N, nthreads, optimization)
%
% Inputs:
%   N            - Signal length up to which optimization plans are being created 
%   nthreads     - Plans are optimized for specific thread use. The number
%                  of threads should be equal to the number of threads used 
%                  in the actual application.
%   optimization - Should be 'estimate', 'measure', 'patient' or
%                  'exhaustive'. Types tradeoff speed and quality. Estimate is the
%                  quickest to calculate, but often does not find the optimal optimization plan. 
%                  On the other hand, 'exhaustive' definitely finds the best optimization
%                  plan, but takes a lot of time due to the exhaustive tree search it
%                  performs. 
%
%                  See FFTW's documentation for more details: 
%                  https://www.fftw.org/fftw3_doc/Planner-Flags.html
%
% Output:
%   none 
%
