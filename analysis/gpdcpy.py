#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:51:28 2019

@author: ronaldo
"""
import numpy as np
import control
import scipy


def pdc_function(A,SIG,N,Fs):  
    
    pmax=np.shape(A)[2] #
    
    Am=np.reshape(A,(np.shape(A)[0],np.shape(A)[0]*pmax),order="F") #
    
    M=np.shape(Am)[0] #
    p=int(np.shape(Am)[1]/M) #
    
    tmp3=np.zeros(M,dtype=complex) #
    tmp4=np.zeros(M,dtype=complex) #
    
    Cd=np.diag(np.diag(SIG)) #
    invCd=np.array(np.linalg.inv(Cd),dtype=complex) #
    
    f=np.arange(0,N)*(Fs/(2*N)) #
    z=np.complex(0,1)*2*np.pi/Fs #
    
    A=np.c_[np.identity(M), -Am] #
    
    # Convert to complex
    A=np.array(A,dtype=complex)
    
    
    PDC=np.zeros((M,M,N),dtype=complex) #
    GPDC=np.zeros((M,M,N),dtype=complex) #
    
    for n in range(0,N):
        
        As=np.zeros((M,M),dtype=complex) #
        
        # Aqui ta errado
        for k in range(1,p+2):
           As=As+A[:,(k*M-1)+(1-np.arange(M,0,-1))]*np.exp(-z*(k-1)*f[n]) #
    
        # Aqui ta errado
        for m in range(0,M):
            tmpp1=np.squeeze(As[:,m]) #
            tmp3[m]=np.sqrt(np.transpose(np.conj(tmpp1))@tmpp1) #
            tmp4[m]=np.sqrt(np.transpose(np.conj(tmpp1))@invCd@tmpp1) #   
    
        # Calcula o PDC
        PDC[:,:,n]=As/np.tile(tmp3,(M,1)) 
        # Calcula o GPDC
        GPDC[:,:,n]=(scipy.linalg.sqrtm(invCd)@As)/np.tile(tmp4,(M,1))
    
    return np.absolute(PDC)**2,np.absolute(GPDC)**2,f
    
    
def bootstrap(u,A,ef,nfreq,Fs,nsamps,flag): # Implementar
    
    n,m=np.shape(u)
    
    # Create matrices to receive values
    EB=np.zeros((n,m))
    
    
    pdcBootstrap=np.zeros((n,n,nsamps,nfreq))
        
    
    for i in range(n):
        for j in range(n):
            
            # It makes Var coefficients ij equal to zero for all orders
            Ah0=np.copy(A)
            Ah0[i,j,:]=0
            
            if(j!=i):
                
                for s in range(nsamps):
                    
                    # Resample residuals 
                    for k in range(n):
                        EB[k,:]=np.random.choice(ef[k,:],m)
                        
                    # Generate bootstrap time-series (H0 j->i)
                    Utemp=genvar(u,Ah0,EB)
                
                    # Compute Var 
                    _,pftemp,Atemp,_,_,_,_,_,_ = mvar(Utemp,np.shape(A)[2],2,5)
                    
                    
                    # Compute PDC and GPDC in null hypothesis 
                    
                    if flag=='gpdc':
                       _, GpdcTemp, _=pdc_function(Atemp,pftemp,nfreq,Fs)
                       pdcBootstrap[i,j,s,:]=GpdcTemp[i,j,:]
                    elif flag=='pdc':
                       pdcTemp,_, _=pdc_function(Atemp,pftemp,nfreq,Fs)
                       pdcBootstrap[i,j,s,:]=pdcTemp[i,j,:]
                    
             
    return pdcBootstrap



def genvar(serie,A,e):
    
    n,m=np.shape(serie)
    _,_,p=np.shape(A)
    
    # Only initial conditions
    serie[:,p:]=np.copy(e[:,p:])
    
    for t in range(p,m):
        for k in range(p):
            serie[:,t]=serie[:,t]+A[:,:,k]@serie[:,t-k]
    
    return serie


def significance(pdcBootstrap,alpha,nChannels,nfreq):
    
    significancePDC=np.zeros((nChannels,nChannels,nfreq))
    
    for i in range(nChannels):
        for j in range(nChannels):
            significancePDC[i,j,:]=np.quantile(pdcBootstrap[i,j,:,:],1-alpha,0)

    return significancePDC
 


def mvar(u,minIP,maxIP,alg,criterion):  # Pronta
    
    nSegLength,nChannels=np.transpose(u).shape
    
    if criterion==5:
    
        IP=maxIP
        if alg==1: # Nuttal-Strand
            pf,A,pb,B,ef,eb=mcarns(u,IP)
            pf=pf/nSegLength
        elif alg==2: # OLS
            pf,A,ef=cmlsm(u,IP)
            B=[]
            eb=[]
            pb=[]
            pf=pf/nSegLength
        elif alg==3: # Vieira-Morf
            pf,A,pb,B,ef,eb=mcarvm(u,IP)
            pf=pf/nSegLength
       
        vaic=np.shape(u)[1]*np.log(np.linalg.det(pf))+2*(nChannels*nChannels)*IP
        Vaicv=np.copy(vaic)
    
        return IP,pf,A,pb,B,ef,eb,vaic,Vaicv

    vaicv=0
    maxOrder=np.copy(maxIP)
    UpperbounderOrder=np.copy(maxIP)
    IP=minIP
    Vaicv=np.zeros(maxOrder+1)
    
    while IP<=UpperbounderOrder:
        
        if alg==1: # Nuttal-Strand
            npf,na,npb,nb,nef,neb=mcarns(u,IP)  
        elif alg==2: # OLS
            npf,na,nef=cmlsm(u,IP)
        elif alg==3: # Vieira-Morf    
            npf,na,npb,nb,nef,neb=mcarvm(u,IP)
        
        if criterion==1: # AIC 
            vaic=np.shape(u)[1]*np.log(np.linalg.det(npf))\
            +2*(nChannels*nChannels)*IP
        elif criterion==2: # Hannan-Quinn
            vaic=np.shape(u)[1]*np.log(np.linalg.det(npf))\
            +2*np.log(np.log(np.shape(u)[1]))*(nChannels*nChannels)*IP
        elif criterion==3: # Schwartz
            vaic=np.shape(u)[1]*np.log(np.linalg.det(npf))\
            +np.log(np.shape(u)[1])*(nChannels*nChannels)*IP    
        elif criterion==4: # FPE
            vaic=np.log(np.linalg.det(npf)*((np.shape(u)[1]+nChannels*IP+1)/\
            np.shape(u)[1]-nChannels*IP-1)**nChannels)
            
            
        Vaicv[IP]=vaic
        print("IP=%f,vaic=%f\n" % (IP,vaic))
        if (vaic>vaicv) and (IP!=minIP):
            
            vaic=np.copy(vaicv)
            
            break   
        
        vaicv=np.copy(vaic)
        pf=np.copy(npf)
        A=np.copy(na)
        ef=np.copy(nef)

        if alg==1:     
            B=np.copy(nb)
            eb=np.copy(neb)
            pb=np.copy(npb)
        else:
            B=[]
            eb=[]
            pb=[]
        
        IP=IP+1
        
    IP=IP-1
    vaic=np.copy(vaicv)
    Vaicv=np.copy(Vaicv[1:(IP+1)])
    

    pf=pf/nSegLength        
            
    return IP,pf,A,pb,B,ef,eb,vaic,Vaicv       
            
            
def cmlsm(u,IP): # Pronta
    
    m,n=u.shape
    bvar,SU,e=mlsmx(u,IP)
    
    na=np.reshape(bvar,(m,m,IP),order="F")
    npf=SU*n
    nef=e
    
    return npf,na,nef

def mlsmx(Y,p): # Pronta
    
    K,T=Y.shape
    Z=zmatrm(Y,p)
    Gamma=Z@np.transpose(Z)
    U1=np.linalg.inv(Gamma)@Z
    SU=Y@np.transpose(Y)-Y@np.transpose(Z)@U1@np.transpose(Y)
    SU=SU/(T-K*p-1);
    bvar=np.kron(U1,np.identity(K))@np.reshape(Y,K*T,order="F") #b#
    etemp=np.reshape(Y,K*T,order="F")-np.kron(np.transpose(Z),np.identity(K))@bvar
    e=np.reshape(etemp,(K,T),order="F")
    
    return bvar,SU,e

def zmatrm(Y,p): # Pronta
    
    K,T=Y.shape
    y1=np.r_[np.zeros(K*p), np.reshape(np.flipud(Y),K*T,order="F")]
    Z=np.zeros((K*p,T))
    
    for i in range(0,T):
        Z[:,i]=np.flipud(y1[K*i:K*i+K*p]) 
    
    return Z


def mcarns(u,IP): # Pronta
    
    lx,cx=u.shape
    
    NUMCHS=lx         
    N=max(u.shape)    

    
    ef=np.copy(u)       
    eb=np.copy(u)       
    pf=u@np.transpose(u);     
    pb=np.copy(pf)       

    M=-1
    flag=0
    
    
    A=np.zeros((NUMCHS,NUMCHS,IP))
    B=np.zeros((NUMCHS,NUMCHS,IP))
    
    while flag==0:
        
        pfhat=ef[:,M+2:N]@np.transpose(ef[:,M+2:N])      
        pbhat=eb[:,M+1:N-1]@np.transpose(eb[:,M+1:N-1])  
        pfbhat=ef[:,M+2:N]@np.transpose(eb[:,M+1:N-1])   
        
        M=M+1
        
        RHO=control.lyap(pfhat@np.linalg.inv(pf),np.linalg.inv(pb)@pbhat,-2*pfbhat) 
        
        AM=-RHO@np.linalg.inv(pb)                    
        BM=np.transpose(-RHO)@np.linalg.inv(pf)  
        
        A[:,:,M]=AM  
        B[:,:,M]=BM  
        
        pf=pf-AM@BM@pf 
        pb=pb-BM@AM@pb 
        
        
        if M != 0:
            for K in range(0,M):
                    temp1=np.copy(A[:,:,K])
                    A[:,:,K]=A[:,:,K]+AM@B[:,:,(M-1)-K]   
                    B[:,:,(M-1)-K]=B[:,:,(M-1)-K]+BM@temp1 
        
        
        efTemp=np.fliplr(ef[:,M+1:N]) 
        ebTemp=np.fliplr(eb[:,M:N-1]) 
        
        matrix1Temp=efTemp+AM@ebTemp 
        matrix2Temp=ebTemp+BM@efTemp 
        
        ef[:,M+1:N]=np.fliplr(matrix1Temp)   
        eb[:,M+1:N]=np.fliplr(matrix2Temp)  
       
       
        if M==(IP-1):
            A=-A
            B=-B
            flag=1
      
   
    return pf,A,pb,B,ef,eb
    
def mcarvm(u,IP): # Pronta
    
    lx,cx=u.shape
    
    NUMCHS=lx         
    N=max(u.shape)    

    
    ef=np.copy(u)       
    eb=np.copy(u)       
    pf=u@np.transpose(u);     
    pb=np.copy(pf)       

    M=-1
    flag=0
    
    
    A=np.zeros((NUMCHS,NUMCHS,IP))
    B=np.zeros((NUMCHS,NUMCHS,IP))
    
    while flag==0:
        
        pfhat=ef[:,M+2:N]@np.transpose(ef[:,M+2:N])      
        pbhat=eb[:,M+1:N-1]@np.transpose(eb[:,M+1:N-1])  
        pfbhat=ef[:,M+2:N]@np.transpose(eb[:,M+1:N-1])   
        
        M=M+1
        
        Spfhat=scipy.linalg.sqrtm(pfhat)
        Spbhat=scipy.linalg.sqrtm(pbhat)
        ISpfhat=np.linalg.inv(Spfhat)
        ISpbhat=np.linalg.inv(Spbhat)
        RHO=ISpfhat@pfbhat@np.transpose(ISpbhat)
        
        
        AM=-Spfhat@RHO@ISpbhat
        BM=-Spbhat@np.transpose(RHO)@ISpfhat
        
        A[:,:,M]=AM  
        B[:,:,M]=BM  
        
        pf=pf-AM@BM@pf 
        pb=pb-BM@AM@pb 
        
        
        if M != 0:
            for K in range(0,M):
                    temp1=np.copy(A[:,:,K])
                    A[:,:,K]=A[:,:,K]+AM@B[:,:,(M-1)-K]   
                    B[:,:,(M-1)-K]=B[:,:,(M-1)-K]+BM@temp1 
        
        
        efTemp=np.fliplr(ef[:,M+1:N]) 
        ebTemp=np.fliplr(eb[:,M:N-1]) 
        
        matrix1Temp=efTemp+AM@ebTemp 
        matrix2Temp=ebTemp+BM@efTemp 
        
        ef[:,M+1:N]=np.fliplr(matrix1Temp)   
        eb[:,M+1:N]=np.fliplr(matrix2Temp)  
       
       
        if M==(IP-1):
            A=-A
            B=-B
            flag=1
      
   
    return pf,A,pb,B,ef,eb   
    
    
    
