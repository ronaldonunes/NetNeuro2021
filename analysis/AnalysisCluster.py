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
from scipy.io import loadmat

# Inputs
seedIni=int(sys.argv[1])
seedEnd=int(sys.argv[2])

# Clusters
frontoparietal=np.array([1, 6, 5,4, 0, 2, 3]) 
visual=np.array([18, 16, 15, 13, 9, 12, 8])     

def computeGPDCCluster(areas,clusterName,seed):
    
    # Parameters for analysis
    p=analysisParameters.parameters(seed)
    
    # Fln
    conn = loadmat(p["Flnpath"])
    conn=conn['Fln']
    
    # Arrays to store data
    lfpMatrix=np.zeros((p['Nareas'],round((p['tf']/p['dt'])/(p['fs']/p['fs_new']))))
    
    # Load LFP signals 
    for i in range(p['Nareas']):
        lfp=np.load(p['path']+'LFP' +'_'+str(i+1)+'.npy')
        # Pre processing
        lfpMatrix[i,:]=methods.setLFP(lfp,p['transSteps'],p['fs'],p['fs_new'])  

    # new Conn
    connTemp=conn[areas,:]
    connTemp=connTemp[:,areas]
    
    # Select LFP signals
    lfpMatrixTemp=lfpMatrix[areas,:]
    
    # GPDC 
    print('GPDC')    
    gpdcMatrix,fGpdc,_=methods.gpdc(lfpMatrixTemp,1,p['maxIP'],p['alg'],p['criterion'],p['fs_new'],p['nfreqGPDC'])
    
    # Save Data
    np.savez(p['path']+'gpdc_'+clusterName+'.npz',np.max(gpdcMatrix,axis=2).T,connTemp,areas)
        
    return 0

for seed in range(seedIni, seedEnd+1):
    
    computeGPDCCluster(frontoparietal,'frontoparietal',seed)
    computeGPDCCluster(visual,'visual',seed)
    

         
