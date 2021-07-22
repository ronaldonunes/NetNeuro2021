% Firing Rate


function [frE,frI]=firingRate(spk,N,Ne,timeIni,timeEnd,dt,window,i)

    % Convert spike times to second
    spk(:,1)=spk(:,1)*10^(-3);
    % Excitatory spikes
    spkEx=spk(spk(:,2)<Ne+(N*(i-1)),:);
    % Inhibitory spikes
    spkIn=spk(spk(:,2)>=Ne+(N*(i-1)),:);
    
    frE=firingRateCompute(spkEx,Ne,timeIni,timeEnd,dt,window);
    frI=firingRateCompute(spkIn,N-Ne,timeIni,timeEnd,dt,window);

end

function fr=firingRateCompute(spk,n,timeIni,timeEnd,dt,window)

    % Array to firing rate 
    firing_rate=zeros(3,(timeEnd-timeIni)/dt);

    % Store times 
    firing_rate(1,:)=linspace(timeIni,timeEnd,timeEnd/dt);

    % Number of spikes
    spkTimes=unique(spk(:,1));
    Ncount = histc(spk(:,1), spkTimes);

    % I tried several times using vectorization but I did not get it
    % It is not important now (keep going with this Gambiarra).
    %% Rever esse trecho
    k=1;
    for i=1:length(firing_rate) 

        if firing_rate(1,i)>spkTimes(k)
            % (https://neuronaldynamics.epfl.ch/online/Ch12.S2.html)
            firing_rate(2,i-1)=Ncount(k)*(1/dt/n);
            k=k+1;
        end
        %  after last spikeTime get out from the loop
        if k>length(spkTimes)  
            break;
        end
        
    end

    firing_rate(3,:)=slidingWindow(firing_rate(2,:),window,dt);
    fr=firing_rate;
    
end


function sw=slidingWindow(serie,width,dt)

    width_dt = (width / 2 / dt)*2 + 1;
    window = ones(1,width_dt);

    sw=conv(serie, window * 1. / sum(window), 'same');

end

