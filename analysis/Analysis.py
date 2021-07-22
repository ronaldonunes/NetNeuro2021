#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 09:31:38 2020

Apply measures to and save data 

@author: ronaldo
"""
from brian2 import *
import numpy as np
import sys
import methods
import analysisParameters
import time

def analysis(seed):
    
    # Parameters for analysis
    p=analysisParameters.parameters(seed)
    
    # Arrays to store data
    lfpMatrix=np.zeros((p['Nareas'],round((p['tf']/p['dt'])/(p['fs']/p['fs_new']))))
    gpdcMatrix=np.zeros((p['Nareas'],p['Nareas'],p['nfreqGPDC']))
    gpdcPMatrix=np.zeros((p['Nareas'],p['Nareas'],p['nfreqGPDC']))
    frMatrix=np.zeros((p['Nareas'],3)) # Ex, In, Total
    synMatrix=np.zeros((p['Nareas'],3)) # Ex, In, Total
    
    
    # Load LFP signals 
    for i in range(p['Nareas']):
        lfp=np.load(p['path']+'LFP' +'_'+str(i+1)+'.npy')
        # Pre processing
        lfpMatrix[i,:]=methods.setLFP(lfp,p['transSteps'],p['fs'],p['fs_new']) 
    
    # Save processed LFP 	
    np.save(p['path']+'lfpDownsampled.npy',lfpMatrix)	

    # GPDC 
    print('GPDC')    
    gpdcMatrix,fGpdc,order=methods.gpdc(lfpMatrix,p['minIP'],p['maxIP'],p['alg'],p['criterion'],p['fs_new'],p['nfreqGPDC'])
    np.savez(p['path']+'gpdc.npz',fGpdc,gpdcMatrix)
    
    # Save AIC order 	
    np.save(p['path']+'orderAIC.npy',order)

    # Firing Rate 
    print('Firing Rate')
    frMatrix=methods.firingRate(p['t0'],p['tf'],p['transS'],p['transMs'],p['dt'],p['path'],p['Nareas'],p['N'],p['Ne'])
    np.savez(p['path']+'firingRate.npz',frMatrix)
    
    # # Spike Synchronization
    print('Spike Synchronization')
    synMatrix=methods.synchronization(seed,p['transMs'],p['path'],p['Nareas'],p['N'],p['Ne'],p['timeRange'])
    np.savez(p['path']+'synchronization.npz',synMatrix)


# Inputs
seedIni=int(sys.argv[1])
seedEnd=int(sys.argv[2])

for seed in range(seedIni, seedEnd+1):
    t0 = time.time()
    analysis(seed)
    t1 = time.time()
    print(t1-t0)                
