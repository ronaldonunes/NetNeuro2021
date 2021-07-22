timeIni=0; % seconds
timeEnd=6; % seconds
window=100*10^(-3); % seconds
dt=0.05*10^(-3); % seconds
N=2000; % Total number of neurons
Ne=1600; % number of excitatory neurons
s1=0;
v1=3;


% Plot Raster Plot
figure(3)
for i=14
    b = readNPY(strcat('/home/ronaldo/Desktop/SNN/deterministic/V1_',num2str(v1),'/S1_',num2str(s1),'/Spikes_',num2str(i),'.npy'));
    plot(b(:,1)/1000,b(:,2),'b.','MarkerSize',0.1)
    xlim([1 6])
end

% Plot Firing Rate
figure(4)
for i=14
    [frE,frI]=firingRate(b,N,Ne,timeIni,timeEnd,dt,window,i);
   
    plot(dt:dt:timeEnd,frE(3,:),'b');
    hold on
    plot(dt:dt:timeEnd,frI(3,:),'r');
    xlim([1 6])
end

figure(5)
for i=14
    c = readNPY(strcat('/home/ronaldo/Desktop/SNN/deterministic/V1_',num2str(v1),'/S1_',num2str(s1),'/stimulus.npy'));
    plot(dt:dt:timeEnd,c./max(c),'k','LineWidth',3)
    ylim([0 1])
    xlim([1 6])
end


