#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 10:59:12 2020

@author: ronaldo
"""
from brian2 import *
from partialCorrelation import *
from plv import *
from mne.time_frequency import psd_array_multitaper
import numpy as np
import scipy as sp
import gpdcpy 
import mmse
import spikeSyn
import time

def setLFP(lfp,transiente,freq_inicial,freq_final):
    # Details about filter filtfilt 
    # Check: https://dsp.stackexchange.com/questions/11466/
                # differences-between-python-and-matlab-filtfilt-function
    
    if lfp.ndim==1:
        lfp=np.reshape(lfp, (lfp.size,1))
    
    
    lfp=lfp[transiente:,:]-mean(lfp[transiente:,:],axis=0)
    
    # filtro 4a ordem
    b,a=sp.signal.butter(4,freq_final/(freq_inicial/2))
    lfp=sp.signal.filtfilt(b,a,lfp,axis=0,padlen=3*(max(len(b),len(a))-1))
    
    # downsample
    ### Arrumar esse slicing
    lfp=lfp[::int(freq_inicial/freq_final),:]
    
    return np.squeeze(lfp)

def computeFiringRate(spikeMonitor,group,time_ini,time_end,dt):       
    # Compared to PopulationRateMonitor
    
    
    # 2 -> time and firingRate Background
    firing_rate=np.zeros((2,int(((time_end-time_ini)/second)/(dt/second))))
    
    # vetor com tempos
    firing_rate[0][:]=np.arange(int(time_ini/ms),int(time_end/ms),dt/ms)
           
    tempList=[]
    
    for idx in range(0,len(group)):
        
        tempList.append(np.where(spikeMonitor[:,1]==group[idx])[0].tolist()[:])
    
    tempList2 = [val for sublist in tempList for val in sublist] #### ENTENDER
    tempArray=np.copy(np.asarray(tempList2))

    tempArray=np.copy(spikeMonitor[tempArray])
    times=np.sort(tempArray[:,0])  

    # tempos dos disparos e quantidade de neuronios por tempo
    frTemp=np.unique(times,return_counts=True)

    # indices dos tempos dos disparos
    frTemp2=(np.round(frTemp[0]/(dt/ms))).astype(int)

    # Firing Rate
    firing_rate[1][frTemp2-(int(time_ini/dt))]=frTemp[1]* (1.0/dt/second)/shape(group)[0]

    # Average Firing Rate in the interval 
    averageFR=sum(frTemp[1])/(group.shape[0]*((time_end-time_ini)/second))

    return firing_rate, averageFR


def sliding_window(serie,width,dt):
    # width in ms   
    width_dt = int(width / 2 / dt)*2 + 1
    used_width = width_dt * dt
    window = np.ones(width_dt)
    
    return np.convolve(serie, window * 1. / sum(window), mode='same')

def gpdc(lfp,minIP,maxIP,alg,criterion,fs_new,nfreq):

    order,pf,A,_,_,ef,_,_,_=gpdcpy.mvar(lfp,minIP,maxIP,alg,criterion)
    if order==minIP:
        order,pf,A,_,_,ef,_,_,_=gpdcpy.mvar(lfp,1,maxIP,alg,criterion)
        
    _,gpdcMatrix,f=gpdcpy.pdc_function(A,pf,nfreq,fs_new) 

    return gpdcMatrix,f,order



def firingRate(t0,tf,transS,transMs,dt,path,Nareas,N,Ne):

    
    frMatrix=np.zeros((Nareas,3)) # Ex, In, Total

    for i in range(Nareas):
        #i=13 ######
        spike=np.load(path+'Spikes' +'_'+str(i+1)+'.npy')
        
        # Spikes after transient
        spike=spike[spike[:,0]>transMs,:]
        
        groupEx=np.arange(i*N,i*N+Ne)
        groupIn=np.arange(i*N+Ne,(i+1)*N)
        groupTotal=np.arange(i*N,(i+1)*N)
        
        # Excitatory
        teste,avrFrEx=computeFiringRate(spike,groupEx,t0+transS,tf+transS,dt)
        # Inhibitory
        _,avrFrIn=computeFiringRate(spike,groupIn,t0+transS,tf+transS,dt)
        # Total
        _,avrFrTotal=computeFiringRate(spike,groupTotal,t0+transS,tf+transS,dt)
    
        frMatrix[i,:]=np.array([avrFrEx,avrFrIn,avrFrTotal])

    return frMatrix

def synchronization(seed,transMs,path,Nareas,N,Ne,timeRange):

    synMatrix=np.zeros((Nareas,3)) # Ex, In, Tota
    
    for i in range(Nareas):

        spike=np.load(path+'Spikes' +'_'+str(i+1)+'.npy')
        
        # Spikes after transient
        spike=spike[spike[:,0]>transMs,:]

        groupEx=np.arange(i*N,i*N+Ne)
        groupIn=np.arange(i*N+Ne,(i+1)*N)
        groupTotal=np.arange(i*N,(i+1)*N)

        # Neurons
        neuronsTotal=np.unique(spike[:,1])
        neuronsEx=np.intersect1d(neuronsTotal,groupEx)
        neuronsIn=np.intersect1d(neuronsTotal,groupIn)
        
        spkSynEx=spikeSyn.synchronization(seed,neuronsEx,spike,timeRange)
        spkSynIn=spikeSyn.synchronization(seed,neuronsIn,spike,timeRange)
        spkSynTotal=spikeSyn.synchronization(seed,neuronsTotal,spike,timeRange)
        
        synMatrix[i,:]=np.array([spkSynEx,spkSynIn,spkSynTotal])

    return synMatrix
       


