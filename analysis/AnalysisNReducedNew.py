#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 17:57:47 2020

@author: ronaldo
"""


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
from itertools import combinations 

# Inputs
seedIni=int(sys.argv[1])
seedEnd=int(sys.argv[2])
nAreas=int(sys.argv[3])


# Select all combination of 19 in nAreas  
comb = combinations(np.arange(19), nAreas)
# Convert to array
arrayAreas = np.array(list(comb))
# Select 150 of these combinations
trials=np.random.choice(len(arrayAreas),150,replace=False)


for k,trial in enumerate(trials):

    areas=arrayAreas[trial,:]
    
    for seed in range(seedIni, seedEnd+1):
    
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
        gpdcMatrix,fGpdc,order=methods.gpdc(lfpMatrixTemp,p['minIP'],p['maxIP'],p['alg'],p['criterion'],p['fs_new'],p['nfreqGPDC'])
    
        # Save Data
        np.savez(p['path']+'NReducedGpdc_'+str(nAreas)+'areas_trial'+str(k)+'.npz',np.max(gpdcMatrix,axis=2).T,connTemp,areas,order)
    

         
