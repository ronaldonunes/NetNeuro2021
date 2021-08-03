timeIni=0; % seconds
timeEnd=30; % seconds
window=100*10^(-3); % seconds
dt=0.1*10^(-3); % seconds
N=2000; % Total number of neurons
Ne=1600; % number of excitatory neurons


mydir  = pwd;
idcs   = strfind(mydir,'NetNeuro2021/');
newdir = mydir(1:idcs(end)-1);
%% Plot Average PSD over 5 trials
pxxVector=zeros(10,10001);
    
for j=1:10
    b = readNPY(strcat(newdir,'NetNeuro2021/data/Seed',num2str(j),'/lfpDownsampled.npy'));
    LFP=b(4,1001:2000);
    [pxxVector(j,:),f]=psdensity(LFP,1000,10000,false);
end

boundedline(f,mean(pxxVector,1),std(pxxVector,'',1),'k');
xlim([0 80])
 

meanPSD=mean(pxxVector,1);
stdPSD=std(pxxVector,'',1);
freq=f;


save(strcat(newdir,'NetNeuro2021/data/psd.mat'),'meanPSD','stdPSD','freq')

%% Plot Average Firing Rate over 5 trials
figure(2)    
frEVector=zeros(10,300000);
frIVector=zeros(10,300000);
i=4;
for j=1:10
    b = readNPY(strcat(newdir,'NetNeuro2021/data/Seed',num2str(j),'/Spikes_',num2str(i),'.npy'));
    [frE,frI]=firingRate(b,N,Ne,timeIni,timeEnd,dt,window,i);

    frEVector(j,:)=frE(3,:);
    frIVector(j,:)=frI(3,:);

end

boundedline(dt:dt:timeEnd,mean(frEVector,1),std(frEVector,'',1),'b');
hold on
boundedline(dt:dt:timeEnd,mean(frIVector,1),std(frIVector,'',1),'r');

time=(dt:dt:timeEnd)*(10^3);
frEmean=mean(frEVector,1);
frEstd=std(frEVector,'',1);
frImean=mean(frIVector,1);
frIstd=std(frIVector,'',1);

save(strcat(newdir,'NetNeuro2021/data/firingRate.mat'),'time','frEmean','frEstd','frImean','frIstd')