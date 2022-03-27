# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 2017
Modified on Fri Aug 31 2018
@author: gveraofe
This script includes the generation of the wave spectrum for running long 
crested and short crested wave simulations.
PM = Pierson-Moskovitz
JS = JONSWAP
DSP = DIRECTIONAL SPREADING METHOD.The directional spreading method calculates 
a random wave angle for each frequency based on the wave spreading factor s1.
DPM = TIM CODE FOR RUNNING PM DIRECTIONAL WAVES. IT DOES NOT WORK AT THE MOMENT
DJS = TIM CODE FOR RUNNING JS DIRECTIONAL WAVES. IT DOES NOT WORK AT THE MOMENT
"""

import numpy as np
import random
import scipy.special as spe
import random
#-----------------------------------------------------------------------------
#                              1.PM - WAVE SPECTRUM
#-----------------------------------------------------------------------------
#
# A Pierson - Moskovitz Spectra will be used
def PM(NEM_ini):
    NEM_ini['Spm'] = []
    NEM_ini['H'] = []
    NEM_ini['amp'] = []
    #----------------------------FREQUENCY SPECTRUM--------------------------
    for fp,Hs,Tp,f in zip(NEM_ini['fp'],NEM_ini['Hs'],NEM_ini['Tp'],NEM_ini['f']):#Zip iterates over Hs and Tp as a pair of Tupples
        df = (np.max(f)-np.min(f))/(np.shape(f)[0]-1)
        B = 5.0/16.0 * Hs**2.0/Tp**4
        C = 5.0/4.0 * 1.0/Tp**4
        Spm=([B / fii**5 *np.exp( - C/ fii**4) for fii in f])
        Hi = [2.0 * np.sqrt(2.0*Spm*df) for Spm in Spm]
        amp = [Hi/2.0 for Hi in Hi]
        NEM_ini['Spm'].append(Spm)
        NEM_ini['H'].append(Hi)
        NEM_ini['amp'].append(amp)
               
    NEM_ini['Hs_m01'] = []
    for ii in range(len(NEM_ini['fp'])):
        m01 = 0.0

        for ff in range(len(NEM_ini['f'][ii])-1):
            m01 = m01 + (NEM_ini['f'][ii][ff]-NEM_ini['f'][ii][ff+1] )*NEM_ini['Spm'][ii][ff]
        
        Hs_m01 = 4*np.sqrt(m01)    
        NEM_ini['Hs_m01'].append(Hs_m01)
    
# A JONSWAP Spectra will be used
def JS(NEM_ini):
    NEM_ini['Spm'] = []
    NEM_ini['H'] = []
    NEM_ini['amp'] = []
    #----------------------------FREQUENCY SPECTRUM--------------------------
    for fp,Hs,Tp,f in zip(NEM_ini['fp'],NEM_ini['Hs'],NEM_ini['Tp'],NEM_ini['f']):
        df = (np.max(f)-np.min(f))/(np.shape(f)[0]-1)
        B=Hs**2.0/Tp**4# the 5/16 term is included already on the alpha calculation
        C = 5.0/4.0 * 1.0/Tp**4
        gamma = 3.3
        alpha = 0.0624/(0.230 + 0.0336*gamma - (0.185/(1.9+gamma)))
        beta = [np.exp(-(fii-fp)**2/(2*0.07**2*fp**2)) if (fii <= fp) else np.exp(-(fii-fp)**2/(2*0.09**2*fp**2)) for fii in f]
        Spm=([(alpha* B / fii**5 *np.exp( - C/ fii**4) * gamma **beta) for fii,beta in zip(f,beta)])
        Hi = [2.0 * np.sqrt(2.0*Spm*df) for Spm in Spm]
        amp = [Hi/2.0 for Hi in Hi]
        NEM_ini['Spm'].append(Spm)
        NEM_ini['H'].append(Hi)
        NEM_ini['amp'].append(amp)
        
    NEM_ini['Hs_m01'] = []
    for ii in range(len(NEM_ini['fp'])):
        m01 = 0.0

        for ff in range(len(NEM_ini['f'][ii])-1):
            m01 = m01 + (NEM_ini['f'][ii][ff]-NEM_ini['f'][ii][ff+1] )*NEM_ini['Spm'][ii][ff]
        
        Hs_m01 = 4*np.sqrt(m01)    
        NEM_ini['Hs_m01'].append(Hs_m01)

    return(NEM_ini)

# DSM DIRECTIONAL SPREDING METHOD
def DSM(NEM_ini):
    #THIS SCRIPT IS COPIED FROM THE MW VERSION OF PVI TO CREATE DIRECTIONAL WAVES TO REPLICATE THE EXACT WAVE DIRECTIONS FOR NEMOH INPUT
    sector = []    
    for ii in range(len(NEM_ini['s1'])):
        if NEM_ini['s1'][ii] <= 3.0:
            sector_aux = np.pi
        elif NEM_ini['s1'][ii] <= 15.0:
            sector_aux = np.pi / 2.0
        elif NEM_ini['s1'][ii] <= 50.0:
            sector_aux = np.pi / 4.0
        else:
            sector_aux = np.pi / 8.0
        sector.append(sector_aux)
        
    index = list(range(NEM_ini['Nf']))
    thetam = []
    for sec in sector:
        thetam.append([sec*(index/float(NEM_ini['Nf']-1)-0.5) for index in index])

   
    CDF = []
    DSF = []
    cos2s1 = []
    for theta_m,s1 in zip(thetam,NEM_ini['s1']):
        CDF_ini = 0
        CDF_aux =np.zeros(NEM_ini['Nf'])
        cos2s1=[np.cos(theta) for theta in theta_m]
        DSF = [cos **(2*s1) if cos >0 else 0 for cos in cos2s1]
        for ii in range(len(DSF)):
            CDF_aux[ii] =CDF_ini +DSF[ii]
            CDF_ini = CDF_aux[ii]
        CDF.append(CDF_aux)

    CDF = [CDF/CDF[NEM_ini['Nf']-1] for CDF in CDF]

    thetan = []
    for theta,CDF_,deg in zip(thetam,CDF,NEM_ini['deg_main']):
        thetan_aux = np.zeros(NEM_ini['Nf'])
        for jj in range (NEM_ini['Nf']):
            rval = random.uniform(0,1)
            for ii in range (1,NEM_ini['Nf'],1):
                if CDF_[ii-1] < rval < CDF_[ii]:
                    thetan_aux[jj] = (theta [ii-1] + (theta[ii]-theta[ii-1])*((rval - CDF_[ii-1])/(CDF_[ii]-CDF_[ii-1])))
        thetan.append(thetan_aux) 
    return(thetan)

#DIRECTIONAL SPECTRUMS FROM TIM CODE AT THE MOMENT THEY HAVE NO USE
'''
# A DIRECTIONAL JONSWAP Spectra will be used
# G = D(theta)
def JSDIR(Hs,Tp,Nf,fini,ffin,Ndir,dini,dfin,sig):

    # Discretisation
    freqs = np.linspace(fini,ffin,Nf)/Tp
    T = 1.0/freqs

    #freqs =np.array( [1/Tp])
    sdir = np.linspace(dini,dfin,Ndir)
    
    #JONSWAP SPECTRUM
    #def makeSpec(Hs,Tp,sdir,sig,freqs):#define as funtion
    # Determine necessary parameters    
    fp = 1.0/Tp
    gamma = 3.3
    alpha = 0.0624/(0.230 + 0.0336*gamma - (0.185/(1.9+gamma))) 
    freq_n = len(freqs)
    beta = np.zeros(freq_n)
    specSS = np.zeros(freq_n)
    # Calculate spectrum
    for iS in range(0,freq_n):
        if freqs[iS] < fp:
            sigma = 0.07
        else:
            sigma = 0.09
        beta[iS] = np.exp(-(freqs[iS]-fp)**2/(2*sigma**2*fp**2))
        specSS[iS] = alpha*Hs**2*fp**4*freqs[iS]**(-5)*gamma**beta[iS]*np.exp(-5/4.0*(fp/freqs[iS])**4)
    # Directionality - Longuet
    #C = np.sqrt(np.pi)/(2.0*np.pi)*spe.gamma(sig+1)/spe.gamma(sig+1/2.0)
    #G = C*np.cos((sdir-0.0)/2.0)**(2.0*sig)
    G = np.zeros((len(sdir),len(sig)))
    for ii in range(len(sig)):
        for kk in range(len(sdir)):
            G[kk,ii] = 2.0**(2.0*sig[ii]-1)/np.pi*(spe.gamma(sig[ii]+1)**2.0)/(spe.gamma(2.0*sig[ii]+1))*np.cos((sdir[kk]-0.0)/2.0)**(2.0*sig[ii])    
    
    SG = np.zeros((len(sdir),len(freqs)))
    amp = np.zeros((len(sdir),len(freqs)))
    
    df = (np.max(freqs)-np.min(freqs))/(np.shape(freqs)[0]-1)
    ddir = (np.max(sdir)-np.min(sdir))/(np.shape(sdir)[0]-1)
    
    for iDir in range(len(sdir)):
        SG[iDir,:] = specSS*G[iDir]
        amp[iDir,:] = np.sqrt(2.0*specSS*G[iDir]*df*ddir)
        
        
    S = SG
    
    ##CALCULTE Hs for a single frequency in all directions 
    
    m01 = 0.0
    for ii in range(0,np.shape(S)[1]-1):
        m01 = m01 + (freqs[ii+1]-freqs[ii] )*S[19][ii+1]
            
    Hs_m01 = 4*np.sqrt(m01)
    Hi = 2.0 * amp
    
    return(T,freqs,sdir,Hi,amp,Hs_m01,SG,G)

#P-M Spectrum
# G = D(theta)
def PMDIR(Hs,Tp,Nf,fini,ffin,Ndir,dini,dfin,sig):
    
    freqs = np.linspace(fini,ffin,Nf)/Tp
    T = 1.0/freqs
    sdir = np.linspace(dini,dfin,Ndir)
    Spm = np.zeros(Nf)
    #w = np.linspace(0.1,10.0,100)*wp
    Spm = np.zeros((np.shape(freqs)[0]))
    amp = np.zeros((len(sdir),len(freqs)))
    df = (np.max(freqs)-np.min(freqs))/(np.shape(freqs)[0]-1)
    ddir = (np.max(sdir)-np.min(sdir))/(np.shape(sdir)[0]-1)

    #----------------------------FREQUENCY SPECTRUM--------------------------
    B = 5.0/16.0 * Hs**2.0/Tp**4
    C = 5.0/4.0 * 1.0/Tp**4
    for ii in range(np.shape(freqs)[0]):
        Spm[ii] = B / freqs[ii]**5 *np.exp( - C/ freqs[ii]**4)
    
    G = np.zeros((len(sdir),len(sig)))
    for ii in range(len(sig)):
        for kk in range(len(sdir)):
            G[kk,ii] = 2.0**(2.0*sig[ii]-1)/np.pi*(spe.gamma(sig[ii]+1)**2.0)/(spe.gamma(2.0*sig[ii]+1))*np.cos((sdir[kk]-0.0)/2.0)**(2.0*sig[ii])
    
    SG = np.zeros((len(sdir),len(freqs)))
    
    for iDir in range(len(sdir)):
        SG[iDir,:] = Spm*G[iDir]
        amp[iDir,:] = np.sqrt(2.0*Spm*G[iDir]*df*ddir)

    S = SG
        
    m01 = 0.0
    for ii in range(0,np.shape(S)[1]-1):
        m01 = m01 + (freqs[ii+1]-freqs[ii] )*S[19][ii+1]
            
    Hs_m01 = 4*np.sqrt(m01)
    Hi = 2.0 * amp

    return(T,freqs,sdir,Hi,amp,Hs_m01,SG,G)        
'''