#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 17:11:12 2020

@author: ronaldo
"""

from brian2 import *
import os


def setParameters():

    par={
      ## Simulation Parameters
      "runs":62, ##############
      "tSim":0.5*second,
      "trans":10000,
      "dt" : 0.1*ms, ############ 
      ## Network parameters  
      "N":2000,
      "NAreas":19,
      "Ne":0.8,
      "probIntra":0.1,
      "probInter":0.05,
      ## Model Parameters
      "gl":0.025*usiemens,
      "g_na" : 12.5*usiemens,
      "g_kd" : 4.74*usiemens,
      "El" : -65*mV,
      "EK" : -80*mV,
      "ENa" : 40*mV,
      "Ee" : 0*mV,
      "Ei" : -70*mV,
      "tE" : 2*ms,           
      "tI" : 8*ms,
      "CmExcMean" : 0.5*nfarad,
      "CmExcStd" : 0.0*nfarad, 
      "CmInhMean" : 0.25*nfarad, #0.25
      "CmInhStd" : 0.0*nfarad,
      "ratePoisson" : 7300*Hz, #10000
      "rateStim": 8000,
     ## Connection parameters
      "wEEMean":2.5, # Excitatory to Excitatory
      "wEEStd":1.0,  
      "wEIMean":2.5, # Excitatory to Inhibitory # mudei aqui 1 is ok
      "wEIStd":1.0,  
      "wIEMean":240, # Inhibitory to Excitatory
      "wIEStd":10,
      "wIIMean":240, # Inhibitory to Inhibitory
      "wIIStd":10,
      "wInputMean":3.2, # Poisson to Networks 
      "wInputStd":1.0,  
      "muE":50, 
      "muI":25, 
      "wStimV1":0, # Stimulus to V1
      "wStimS1":0, # Stimulus to S1
      "delayMeanRR":1,
      "delayStdRR":0,
      ## Initialization
      "VIniMean": -65*mV,
      "VIniStd": 2*mV,
      "gEIni": 0*nS,
      "gIIni": 0*nS,
      "gExtIni": 0*nS,
       ## Path to FLN
      "Flnpath":os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/data/Connectome/Kennedy.mat',
       ## Path to distance
      "Dpath":os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/data/Connectome/KennedyDistance.mat',
       ## Path to save files
      "pathFiles":os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/data/',
      # Axonal velocity
      "speed": 3.5
        }

    return par


def setEquations():
    
    eqs = Equations('''
    dv/dt = (-gl*(v-El)+Iampa+Igaba+Iext+Iext2-
             g_na*(m*m*m)*h*(v-ENa)-
             g_kd*(n*n*n*n)*(v-EK))/Cm : volt
    
    
    Iampa=ge*(Ee-v) :ampere
    Iext=gext*(Ee-v):ampere
    Iext2=gext2*(Ee-v):ampere
    Igaba=gi*(Ei-v) :ampere
    
    Cm:farad
    
    dge/dt=-ge*(1./taue) :siemens
    dgext/dt=-gext*(1./taue) :siemens
    dgext2/dt=-gext2*(1./taue) :siemens
    dgi/dt=-gi*(1./taui) :siemens
    
    m= alpha_m/(alpha_m+beta_m)    :1
    dh/dt = alpha_h*(1-h)-beta_h*h : 1
    dn/dt = alpha_n*(1-n)-beta_n*n : 1
    
    taue=tE:second
    taui=tI:second
    
    alpha_m = 0.1*(mV**-1)*(v+16*mV)/
             (1.0-exp(-(v+16*mV)/(10*mV)))/ms : Hz
    
    beta_m = 4.0*(exp(-(v+41*mV)/(18*mV)))/ms : Hz
    
    alpha_h = 0.07*(exp(-(v+30*mV)/(20*mV)))/ms : Hz
    
    beta_h = 1./(1.0+exp((-v)/(10*mV)))/ms : Hz
    
    
    alpha_n = 0.01*(mV**-1)*(v+20*mV)/
             (1.0-exp(-(v+20*mV)/(10*mV)))/ms : Hz
             
    beta_n = 0.125*(exp(-(v+30*mV)/(80*mV)))/ms : Hz
    
    ''')
        
    return eqs
