%% psd
% 
%  Power Spectral Density
%
%% Syntax
%
%    [pxx,f]=psdensity(x,Fs,len,flg)
%
%% Arguments
%
%    Input: 
%   
%    x           Signal
%    Fs          Sampling rate
%   'len         number of frequencies
%    flg         flg='true' -> normalized
%
%
%    Output:
%
%    pxx          Power Spectral Densities
%    f            frequencies
%
%
%
%% Description
%
%   This function computes the power spectral density
%
%
%   Autor: Ronaldo Nunes (ronaldovnunes@gmail.com)
%

function [pxx,f]=psdensity(x,Fs,len,flg)

% 4 is the default nw
[pxx_puro,f] = pmtm(x,4,2*len,Fs);


    if flg==true
        pxx_max=max(pxx_puro);
        pxx=pxx_puro/pxx_max;
    else
        pxx=pxx_puro;
    end
    
end