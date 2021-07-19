from brian2 import *
import mouseModelUtils as mmu
from model import *
import time
import os
import sys


#################################### Parameters ############################## 
par=setParameters()

######################################## Seed  ################################

par["seed"]=int(sys.argv[1])       

start=time.time()

################################## Create folder to store data ################
par["pathFiles"]=par["pathFiles"]+"/SimulationData/Seed"+str(par["seed"])+"/"

if not os.path.exists(par["pathFiles"]):
        os.makedirs(par["pathFiles"])

################################### Model Equations ###########################
eqs=setEquations()
#################################### Define Network ###########################
P,Input,CE,CI,Input_P=mmu.setNetwork(eqs,par)
################################## Run large scale model ######################
mmu.nRuns_writeFile(P,Input,CE,CI,Input_P,par)
################################### Save Simulation time ######################
tempo=time.time()-start
f1=open(par['pathFiles']+'Tempo.txt','w+')
f1.write(str(tempo))
f1.close()

f2=open(par['pathFiles']+'Parameters.txt','w+')
f2.write(str(par))
f2.close()
  
