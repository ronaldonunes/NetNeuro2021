%% setLFP
% 
%  Pre-preocessing LFP data
%
%% Syntax
%
%    LFP=setLFP(LFP,freq_inicial,freq_final)
%
%% Argumentos
%
%    Input: 
%   
%    LFP            LFP time-series
%    freq_inicial   Sampling rate before downsample (orignal sampling rate)     
%    freq_final     Sampling rate after downsample 
%    
%    Output 
%    
%    LFP           Pre-processed LFP data
%
%% Descricao
%
%  This function perform remove dc component, perform downsample, filtering
%  (to avoid aliasing) and detrend in LFP data
% 
%
% Autor: Ronaldo Nunes (ronaldovnunes@gmail.com)
% 
%% References
%
% Ref: http://www.mathworks.com/help/signal/ug/...
%      filtering-before-downsampling.html?requestedDomain=www.mathworks.com
%
%
% 

function LFP=setLFP(LFP,freq_inicial,freq_final)

nChannels=size(LFP,2);

% mean (remove DC component)
for i=1:nChannels
    LFP(:,i)=LFP(:,i)-mean(LFP(:,i));
end

% Low-pass filter in freq_final
LFP=filtro(LFP,freq_final,freq_inicial);

% Downsample
LFP=downsample(LFP,freq_inicial/freq_final);

% detrend
LFP=detrend(LFP);

end
