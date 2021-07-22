timeIni=0; % seconds
timeEnd=3; % seconds
window=100*10^(-3); % seconds
dt=0.1*10^(-3); % seconds
N=2000; % Total number of neurons
Ne=1600; % number of excitatory neurons
s1=0;
v1=1;
seed=1;

path='/home/ronaldo/Dados/ProjectUfabc/MouseData/SimulationData/Seed1/';

% Plot LFP for 1 trial
figure(1)
for i=1:19
    b = readNPY(strcat(path,'LFP_',num2str(i),'.npy'));
    LFP=setLFP(b(20000:end),1/dt,200);
    subplot(5,4,i)
    plot(linspace(1,3,length(LFP)),LFP)
    
    if (i<16)
        set(gca,'XTick',[])
    end
    
end

% Plot Raster Plot for 1 trials
figure(2)
for i=1:19
    b = readNPY(strcat(path,'Spikes_',num2str(i),'.npy'));
    subplot(5,4,i)
    plot(b(:,1)/1000,b(:,2),'b.','MarkerSize',0.1)
    xlim([1 3])
    
     if (i<16)
        set(gca,'XTick',[])
    end
    
end

% Plot Average Firing Rate over 5 trials
figure(3)
for i=1:19
    
    frEVector=zeros(5,30000);
    frIVector=zeros(5,30000);
    
    for j=1:5
        path=strcat('/home/ronaldo/Desktop/SNN/deterministic/Seed',num2str(j),'/rate4/V1_0.5/');
        b = readNPY(strcat(path,'/S1_0.0/Spikes_',num2str(i),'.npy'));
        [frE,frI]=firingRate(b,N,Ne,timeIni,timeEnd,dt,window,i);
    
        frEVector(j,:)=frE(3,:);
        frIVector(j,:)=frI(3,:);
        
    end
    
    subplot(5,4,i)
    boundedline(dt:dt:timeEnd,mean(frEVector,1),std(frEVector,'',1),'b');
    hold on
    boundedline(dt:dt:timeEnd,mean(frIVector,1),std(frIVector,'',1),'r');
    xlim([1 3])
    
     if (i<16)
        set(gca,'XTick',[])
    end
    
end

% Plot AVerage PSD over 5 trials
figure(4)
for i=1:19
    
    pxxVector=zeros(5,50001);
    
    for j=1:5
        path=strcat('/home/ronaldo/Desktop/SNN/deterministic/Seed',num2str(j),'/rate4/V1_0.5/');
        b = readNPY(strcat(path,'/S1_0.0/LFP_',num2str(i),'.npy'));
        LFP=setLFP(b(20000:end),1/dt,200);
        [pxxVector(j,:),f]=psdensity(LFP,200,50000,false);
    end
    
    subplot(5,4,i)
    boundedline(f,mean(pxxVector,1),std(pxxVector,'',1),'k');
    
     if (i<16)
        set(gca,'XTick',[])
    end
    
    
end