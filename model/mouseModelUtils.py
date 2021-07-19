#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"1;4205;0c""
#Created on Wed Feb  5 07:50:57 2020

#@author: ronaldo
#"""
from scipy.io import loadmat
import numpy as np
from brian2 import *
import os
import scipy.stats as stats
# Model SNN utils
import time

prefs.codegen.target = 'cython'  # use the Python fallback



def setConnections(par):
    
    np.random.seed(par["seed"])
    seed(par["seed"])
    
    # N > Number of neurons in each area
    # NAreas > Number of areas
    # probRR > probability for recurrent connections
    # probLR > probability for long range connections
    # path > path for FLN mat file
           
    sources=[] # to store the sources
    targets=[] # to store the targets
    wExc=[]    # to store the excitatory synaptic weights
    wInh=[]    # to store the inhibitory synaptic weights 
    delays=[]  # to store the delays
    
    # all neurons 
    index=np.arange(par["N"]*par["NAreas"])
    # neurons separated by networks
    nets=np.split(index,par["NAreas"])
    
    # Array with all excitatory neurons
    # Array with all inhibitory neurons
    tempN=np.split(nets,[int(par["Ne"]*par["N"]),par["N"]],axis=1)[0:2]
    neuronsE=tempN[0].flatten()
    neuronsI=tempN[1].flatten()
    
    # FLN
    conn = loadmat(par["Flnpath"])
    conn=conn['Fln']
    
    # Distance between areas
    dist = loadmat(par["Dpath"])
    dist = dist['distanceMatrix']
    
    # Delay inter-area
    dInter = dist/par['speed'] #[ms]
    
    for areaSource in range(par["NAreas"]):
         
        for neuron in nets[areaSource]:
            
            # Temp for weights
            tempWeights=[] 
            
            # Temp delays
            tempDelays=[]
            
            ######################### Recurrent connections ###################           
            
            netTemp=nets[areaSource][:]
            
            # to avoid autapses
            netTemp=np.delete(netTemp,np.where(nets[areaSource]==neuron))
            
            # compute connections
            tempConn=np.random.rand(len(netTemp))<par["probIntra"]
            
            # save sources
            sources.extend([neuron]*len(netTemp[tempConn]))
            
            # save targets
            targets.extend(netTemp[tempConn])
            
            ##################### Set Weights #################################
            if neuron<int(par["Ne"]*(nets[0][-1]+1)+nets[areaSource][0]): # Excitatory
                
                # array to store excitatory weights for local connections
                wExcTemp=np.ones(len(netTemp[tempConn]))
                
                # E->E
                _,idxE,_=np.intersect1d(netTemp[tempConn],neuronsE,return_indices=True)
                # Gaussian ditribution of weights
                Wtemp=np.random.randn(len(wExcTemp[idxE]))*par["wEEStd"]+par["wEEMean"] 
                # Truncated
                Wtemp[Wtemp<0]=0
                # Set E->E weights
                wExcTemp[idxE]=Wtemp
                
                # E->I 
                _,idxI,_=np.intersect1d(netTemp[tempConn],neuronsI,return_indices=True)
                # Gaussian ditribution of weights
                Wtemp=np.random.randn(len(wExcTemp[idxI]))*par["wEIStd"]+par["wEIMean"]
                # Truncated
                Wtemp[Wtemp<0]=0
                # Set E->I weights
                wExcTemp[idxI]=Wtemp
                
                # local excitatory weights
                wExc.extend(wExcTemp)
                
                # set at zero the inhibitory connections between these neurons     
                wInh.extend([0]*len(netTemp[tempConn]))
                
                # delays for excitatory local connections
                #delays.extend(np.random.randn(len(netTemp[tempConn]))*par["delayStdRR"]+par["delayMeanRR"])
                delays.extend(np.random.uniform(par["delayStdRR"]-par["delayMeanRR"],par["delayStdRR"]+par["delayMeanRR"],len(netTemp[tempConn])))
                
                ############################### Long Range ###################
                
                # delete area source from array of areas
                netsTarget=np.delete(nets,areaSource,0) 
                
                # connections 
                temp=np.random.rand(np.shape(netsTarget)[0],np.shape(netsTarget)[1])<par["probInter"]
                
                # target neurons
                targets.extend(netsTarget[temp])
                
                # source neuron
                sources.extend([neuron]*len(netsTarget[temp]))
                
                # number of target neurons in each area
                nNeuronsTarget=np.sum(temp*1,1)
                
                # Weight for long range connections
                Weights=conn[areaSource,np.delete(np.arange(par["NAreas"]),areaSource,0)] 
                
                # Store weights for long range connections
                [tempWeights.extend(Weights[i].repeat(nNeuronsTarget[i])) for i in range(nNeuronsTarget.size)][0]
                
                # Convert list to array
                tempWeights=np.asarray(tempWeights)
                
                ############## Set Weights according to target ################
                
                # Return idx for connection where the target is excitatory
                _,idxE,_=np.intersect1d(netsTarget[temp],neuronsE,return_indices=True)
                
                # Return idx for connection where the target is inhibitory
                _,idxI,_=np.intersect1d(netsTarget[temp],neuronsI,return_indices=True)
                
                # Adjust long range connections EE
                tempWeights[idxE]=tempWeights[idxE]*par["muE"]
                
                # Adjust long range connections EI
                tempWeights[idxI]=tempWeights[idxI]*par["muI"]
                
                # Store weights
                wExc.extend(tempWeights)
                
                # Set at zero weights for inhibitory connections
                wInh.extend([0]*np.sum(nNeuronsTarget))
                
                ############## Set Delays according to target #################
                
                # Delays for long range connections
                D=dInter[areaSource,np.delete(np.arange(par["NAreas"]),areaSource,0)] 
                
                # Store weights for long range connections
                [tempDelays.extend(D[i].repeat(nNeuronsTarget[i])) for i in range(nNeuronsTarget.size)][0]

                delays.extend(tempDelays)
                
                ################################################################
                    
            else:
                 # array to store inhibitory weights for local connections
                wInhTemp=np.ones(len(netTemp[tempConn]))
                
                # I->E
                _,idxE,_=np.intersect1d(netTemp[tempConn],neuronsE,return_indices=True)
                # Gaussian ditribution of weights
                Wtemp=np.random.randn(len(wInhTemp[idxE]))*par["wIEStd"]+par["wIEMean"] 
                # Truncated
                Wtemp[Wtemp<0]=0
                # Set I->E weights
                wInhTemp[idxE]=Wtemp
               
                # I->I 
                _,idxI,_=np.intersect1d(netTemp[tempConn],neuronsI,return_indices=True)
                 # Gaussian ditribution of weights
                Wtemp=np.random.randn(len(wInhTemp[idxI]))*par["wIIStd"]+par["wIIMean"] 
                # Truncated
                Wtemp[Wtemp<0]=0
                # Set I->E weights
                wInhTemp[idxI]=Wtemp
               
                # set at zero the excitatory connections between these neurons
                wExc.extend([0]*len(netTemp[tempConn]))
                # local inhibitory weights
                wInh.extend(wInhTemp) 
                # delays for inhibitory local connections
                #delays.extend(np.random.randn(len(netTemp[tempConn]))*par["delayStdRR"]+par["delayMeanRR"])
                delays.extend(np.random.uniform(par["delayStdRR"]-par["delayMeanRR"],par["delayStdRR"]+par["delayMeanRR"],len(netTemp[tempConn])))
                
                ###################################################################        
            

    return sources,targets,wExc,wInh,delays


def setNetwork(eqs,par):
    
    #  Set seed
    np.random.seed(par["seed"])
    seed(par["seed"])

    # Set connections
    sources,targets,wExc,wInh,delays=setConnections(par)
    
    # Total number of neurons
    Ntotal=int(par['N']*par['NAreas'])
    
    # Create neuron group
    P = NeuronGroup(Ntotal, model=eqs, threshold='v>-20*mV', refractory='v>-20*mV', \
                  method='exponential_euler',namespace=par)

    # Poisson in all neurons
    Input = PoissonGroup(Ntotal, par['ratePoisson']) 
    Input_P = Synapses(Input, P, model='w : siemens',on_pre='gext+=w')
    Input_P.connect(j='i')
    Input_P.w=(np.random.randn(Ntotal)*par['wInputStd']+par['wInputMean'])*nS
    
    # Excitatory Synapses (RR and LR)
    CE = Synapses(P, P, model='w : siemens',on_pre='ge+=w')
    CE.connect(i=sources, j=targets)
    CE.w=wExc*nS
    CE.delay=delays*ms
    
    # Inhibitory Synapses (RR)
    CI = Synapses(P, P, model='w : siemens',on_pre='gi+=w')
    CI.connect(i=sources, j=targets)
    CI.w=wInh*nS
    CI.delay=delays*ms
            
    # Initialization
    P.v = np.random.randn(Ntotal)*par['VIniStd']+par['VIniMean']
    P.ge = par['gEIni']
    P.gi = par['gIIni']
    P.gext = par['gExtIni']
    P.gext2 = par['gExtIni']
    
    # Capacitance
    P.Cm=setCapacitance(par)        


    return P,Input,CE,CI,Input_P

def setCapacitance(par):
    
    
    np.random.seed(par["seed"])
    seed(par["seed"])
    
    C=[]
    
    for i in range(par['NAreas']):
        C.extend(np.random.rand(int(par['N']*par['Ne']))*par["CmExcStd"]+par["CmExcMean"]) 
        C.extend(np.random.rand(int(par['N']-par['N']*par['Ne']))*par["CmInhStd"]+par["CmInhMean"])
 
    return C

def nRuns_writeFile(P,Input,CE,CI,Input_P,par):
    
    
    np.random.seed(par["seed"])
    seed(par["seed"])
    
    # File names
    filenamesLFP=[par['pathFiles']+'LFP_'+str(i+1)+'.dat' for i in range(0,par['NAreas'])]
    filenamesSpikes=[par['pathFiles']+'Spikes_'+str(i+1)+'.dat' for i in range(0,par['NAreas'])]
    
    # If file exists, Delete it
    for i in range(0,par['NAreas']):
        if os.path.exists(filenamesLFP[i]):
            os.remove(filenamesLFP[i])
        if os.path.exists(filenamesSpikes[i]):
            os.remove(filenamesSpikes[i])
        
    # Open Files
    fileLFP = [open(filename,'a+') for filename in filenamesLFP]
    fileSpikes = [open(filename,'a+') for filename in filenamesSpikes]
    
    defaultclock.dt=par["dt"]
    
    for i in range(par["runs"]):
           
        ##  Monitors
        mon1=StateMonitor(P,['Iampa', 'Igaba', 'Iext', 'Iext2'],record=True,dt=par["dt"])
        mon2 = SpikeMonitor(P)
        
        #run simulation
        run(par["tSim"], report='text')
        
        splitData(par,mon1,mon2,fileLFP,fileSpikes)
        
        del mon1
        del mon2
        
    # Close Files    
    [file.close() for file in fileLFP]
    [file.close() for file in fileSpikes]
    
    convert(filenamesLFP,filenamesSpikes,par)
    
    return 0          


def splitData(par,monCurrents,monSpikes,fLFP,fSpk):
   
    
    np.random.seed(par["seed"])
    seed(par["seed"])
    
    
    # Neurons index from each network
    index=np.arange(par["N"]*par["NAreas"])
    nets=np.asanyarray(np.split(index,par["NAreas"]))
    
    ## LFP
    # Split array
    # It needs to be converted to ampere because it loses the unit during asanyarrat)
    tempIampa=np.asanyarray(np.split(monCurrents.Iampa,par["NAreas"]))
    tempIgaba=np.asanyarray(np.split(monCurrents.Igaba,par["NAreas"]))
    tempIext=np.asanyarray(np.split(monCurrents.Iext,par["NAreas"]))
    tempIext2=np.asanyarray(np.split(monCurrents.Iext2,par["NAreas"]))
    # Delete information from inhibitory neurons (not necessary to LFP)
    tempIampa=np.delete(tempIampa,np.s_[int(par["N"]*par["Ne"]):int(par["N"])],axis=1)*amp
    tempIgaba=np.delete(tempIgaba,np.s_[int(par["N"]*par["Ne"]):int(par["N"])],axis=1)*amp
    tempIext=np.delete(tempIext,np.s_[int(par["N"]*par["Ne"]):int(par["N"])],axis=1)*amp
    tempIext2=np.delete(tempIext2,np.s_[int(par["N"]*par["Ne"]):int(par["N"])],axis=1)*amp
    # Compute LFP (Mazzoni)
    LFP=1.0*Mohm*(np.mean(np.abs(tempIampa),axis=1)+np.mean(np.abs(tempIgaba),axis=1)+\
        np.mean(np.abs(tempIext),axis=1)+np.mean(np.abs(tempIext2),axis=1))
    
    ## Spikes    
    for i in range(0,par["NAreas"]):
        idx=np.logical_and(monSpikes.i>=nets[i][0], monSpikes.i<=nets[i][-1])
        spikesIdx=monSpikes.i[idx]
        spikesTimes=monSpikes.t[idx]
        
        # Write LFP (in mV) in text file
        np.savetxt(fLFP[i],LFP[i][:]/mvolt)
        # Write Spikes Pop1
        np.savetxt(fSpk[i],np.column_stack([spikesTimes/ms,spikesIdx]))
        
        
    return 0


def convert(fLFP,fSpk,par):
    
    for i in range(0,par['NAreas']):
        if os.path.exists(fLFP[i]):
            lfp=np.loadtxt(fLFP[i])
            np.save(fLFP[i].replace('dat','npy'),lfp) 
            os.remove(fLFP[i])
        if os.path.exists(fSpk[i]):
            spk=np.loadtxt(fSpk[i])
            np.save(fSpk[i].replace('dat','npy'),spk) 
            os.remove(fSpk[i])

    return 0
 


