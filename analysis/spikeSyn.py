#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Spike Sychronization

@author: ronaldo
"""

import pyspike as spk
import os


def synchronization(path,neurons,spike,timeRange):
    
    listTimes=[]
    
    for j in neurons:
        t=spike[spike[:,1]==j,0]    
        listTimes.append(t)
    
    # Save file with spikes
    with open(path+'TempInter.txt', 'w') as f:
        for item in listTimes:
            f.write("%s\n" %  ' '.join(map(str, item)))
            
    spike_trains = spk.load_spike_trains_from_txt(path+'TempInter.txt',edges=timeRange)        
    
    # Compute average spike synchronization        
    avrg_spike_sync_profile = spk.spike_sync_profile(spike_trains)
    avrSyn = avrg_spike_sync_profile.avrg(interval=timeRange)
       
    # Delete TempSpikes
    if os.path.isfile(path+'TempInter.txt'):
        os.remove(path+'TempInter.txt')

    return avrSyn
