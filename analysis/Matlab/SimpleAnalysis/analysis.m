addpath(genpath('npy-matlab-master/'))

seed=1;

timeIni=0; % seconds
timeEnd=3; % seconds
window=50*10^(-3); % seconds
dt=0.05*10^(-3); % seconds
N=2000; % Total number of neurons
Ne=1600; % number of excitatory neurons
v1=5;
s1=0;
% Plot Raster Plot
frEVector=zeros(19,30,5,60000);
frIVector=zeros(19,30,5,60000);

idx=0;
for s1=0:5:145
    idx=idx+1;
    for i=1:19
        for seed=1:5
%             b = readNPY(strcat('/home/ronaldo/Documentos/UvA/MouseKennedySNN/Data/Seed',...
%                 num2str(seed),'/V1_',num2str(v1),'/S1_',num2str(s1),'/Spikes_',num2str(i),'.npy'));


            b = readNPY(strcat('/home/ronaldo/MouseKennedySNN/Data/Seed',...
                num2str(seed),'/V1_',num2str(v1),'/S1_',num2str(s1),'/Spikes_',num2str(i),'.npy'));


            [frE,frI]=firingRate(b,N,Ne,timeIni,timeEnd,dt,window,i);

            frEVector(i,idx,seed,:)=frE(3,:);
            frIVector(i,idx,seed,:)=frI(3,:);  
        end
    end
end


teste=squeeze(mean(frEVector,3));
teste=teste(:,:,41000:49000);
teste=squeeze(mean(teste,3));
