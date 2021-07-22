%% filtro
% 
%  Filter LFP data using a low-pass filter
%
%% Syntax
%
%    filtro(dados,freq_limit,samplingRate)
%
%% Argumentos
%
%    Input: 
%   
%    dados            LFP time-series
%    freq_limit       frequency limit of low-pass filter     
%    samplingRate     Sampling rate of LFP data 
%    
%    Output 
%    
%    LFP           Filtered LFP signal
%
%% Descricao
%
%  This function filter LFP signals using a low-pass filter.
%  The filter is a butterworth with order 6
% 
%
% Autor: Ronaldo Nunes (ronaldovnunes@gmail.com)
% 

function dados=filtro(dados,freq_limit,samplingRate)
    [b,a] = butter(6,freq_limit/(samplingRate/2));
    dados=filtfilt(b,a,dados);
end