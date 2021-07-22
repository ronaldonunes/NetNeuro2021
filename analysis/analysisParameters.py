#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 15:41:52 2020

@author: ronaldo
"""
from brian2 import *
import numpy as np
import os

def parameters(seed):
    
    p={"N":2000, # Number of neurons in each network
      "Ne":1600, # Number of excitatory neurons in each network
      "Nareas":19, # Number of areas
      "t0":0*second, # Initial time
      "tf":30*second, # Final time
      "dt":0.1*ms, # integration step
      "fs":1.0/((0.1*ms)/second), # sampling rate before preprocessing
      "fs_new":1000, # sampling rate after preprocessing
      "transS":1*second, # transient time
      "transMs":1000, # transient in miliseconds
      "transSteps":10000, # transient in miliseconds	
      "timeRange":(1000, 31000), # 1s de transiente
      "maxIP":50, # maximum number of point in the past for GPDC
      "minIP":1, # minimum order to start the searching
      "alg":2, # Algorithm for MVAR (OLS)
      "criterion":1, # Information criteria for MVAR model (AIC) 
      "nfreqGPDC":1000, # Number of the frequencies in GPDC
      ## Path to FLN
      "Flnpath":'/home/ronaldo/github/ProjectUfabc/MouseKennedySNN/Data/Connectome/Kennedy.mat',
      'path': '/home/ronaldo/github/ProjectUfabc/MouseKennedySNN/Data/SimulationData/Seed'+str(seed)+'/', #pathFolder                                  
      
      }

    return p       
       
