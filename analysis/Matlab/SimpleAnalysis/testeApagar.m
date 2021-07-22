timeIni=0; % seconds
timeEnd=31; % seconds
window=100*10^(-3); % seconds
dt=0.1*10^(-3); % seconds
N=2000; % Total number of neurons
Ne=1600; % number of excitatory neurons
s1=0;
v1=1;
seed=1;
transiente=10000;
AreaList={'GU','SSp-un','VISC','MOp','PL','SSs','SSp-bfd','ACAd',...
    'RL','AL','DP','AUDpo','AM','V1','MM','LM','PM','RSPd','P'};

% Plot Raster Plot
figure(1)
for i=1:19
    b = readNPY(strcat('/home/ronaldo/Dados/ProjectUfabc/MouseData/SimulationData/Seed1/Spikes_',num2str(i),'.npy'));
    subplot(5,4,i)
    plot(b(:,1)/(10^3),b(:,2),'b.','MarkerSize',0.1)
    xlim([20 22])
    title(AreaList(i))
    
    if i==17
        xlabel('Time(ms)')
        ylabel('Neuron Index')
    end
end

% % Plot Firing Rate
% figure(2)
% for i=1:19
%     b = readNPY(strcat('/home/ronaldo/Dados/ProjectUfabc/MouseData/SimulationData/Seed1/Spikes_',num2str(i),'.npy'));
%     [frE,frI]=firingRate(b,N,Ne,timeIni,timeEnd,dt,window,i);
%     subplot(5,4,i)
%     plot(frE(1,:),frE(3,:),'b');
%     xlim([20 22])
%     hold on
%     plot(frI(1,:),frI(3,:),'r');
%     xlim([20 22])
%     
% end

% LFP Tratado
figure(4)
for i=1:19
   b = readNPY(strcat('/home/ronaldo/Dados/ProjectUfabc/MouseData/SimulationData/Seed1/LFP_',num2str(i),'.npy'));
    b=setLFP(b(transiente:end,1),1/dt,1000);
    subplot(5,4,i)
    plot(linspace(0,30,length(b(:,1))),b(:,1))
    xlim([20 22])
    title(AreaList(i))
    
       if i==17
        xlabel('Time(ms)')
        ylabel('LFP (mV)')
    end
end


% PSD
figure(5)
for i=1:19
   b = readNPY(strcat('/home/ronaldo/Dados/ProjectUfabc/MouseData/SimulationData/Seed1/LFP_',num2str(i),'.npy'));
    
    subplot(5,4,i)
    b=setLFP(b(transiente:end,1),1/dt,1000);
    [pxxVector,f]=psdensity(b,1000,50000,false);
    plot(f,pxxVector);
    xlim([0 100])
    title(AreaList(i))
    
   if i==17
    xlabel('Frequency (Hz)')
    ylabel('PSD')
   end
end
