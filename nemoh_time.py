# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:11:42 2019

@author: gveraofe
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import random

def irr_wave(NEM_ini,NEM_GRID,NEM_OUT,time,dt):
    time_kd = time
    time = np.linspace(0,time,time/dt)
    kd_time = []
    #TODO KILL ZIP
    for jj in range(len(NEM_ini['Tp'])):
        Tjj = NEM_ini['T'][jj]
        Hjj = NEM_ini['H'][jj]
        ajj = [H/2 for H in Hjj]
        wjj = [2*np.pi*1/T for T in Tjj]       
        Nvar_aux = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']))
        r_phase = np.zeros(len(Tjj))
        for ii in range(len(Tjj)):
            r_phase[ii] = random.uniform(0,2)*np.pi
        for tt in range(np.shape(time)[0]):
            nm_aux = np.zeros((NEM_GRID['Ny'],NEM_GRID['Nx']))
            for ii in range(len(Tjj)):
                nm_tw_P = ajj[ii]*NEM_OUT['tot_wave'][jj][ii] 
                nm_tw_P = np.real(np.abs(nm_tw_P) * np.exp(-1j*(wjj[ii]*time[tt]-r_phase[ii]+np.angle(nm_tw_P))))
                nm_aux += (nm_tw_P)
            Nvar_aux += ((nm_aux)**2)
    
        Nvar = Nvar_aux /time_kd *dt
        kd = 4*np.sqrt(Nvar) / NEM_ini['Hs'][jj]
        kd[kd>=3]=np.nan
        kd_time.append(kd)
    return(kd_time)
    
dirname = "E:/A_Ch_4_Simulations_E/GRID_SIZE/Tp_08.00/dx_L_20"
MW = np.load(os.path.join(dirname,'Irregular_20_MW.npz'))
MW_OUT_2 = MW['MW_OUT']

    
dirname = "C:/Users/gveraofe/Documents/A_Ch_4_Simulations/GRID_SIZE/T_8.00/dx_L_25"
MW = np.load(os.path.join(dirname,'Irregular_20_MW.npz'))
MW_OUT_3 = MW['MW_OUT']

NM = np.load(os.path.join(dirname,'nem_Irregular_20.npz'))
NEM_INI_2 = NM['NEM_ini']
NEM_GRID_2 = NM['NEM_GRID']
NEM_OUT_2 = NM['NEM_OUT']

time = 3999.6 - 30.0

dt = 0.6

kd_time = irr_wave(NEM_INI_2[0],NEM_GRID_2[()],NEM_OUT_2[0],time,dt)
kd_time[0] = kd_time[0] * 1.98 /2.0
plt.figure()
plt.plot(NEM_OUT_2[0]['xvect'][int(201/2.0),:],kd_time[0][int(201/2),:],label = 'NEMOH')
plt.plot(MW_OUT_2[0]['x_vector'][int(169/2),:],MW_OUT_2[0]['Kd_TW'][0][0][int(169/2),:]/MW_OUT_2[0]['Kd_EB'][0][0][int(169/2),:],label = 'MW L/20')
plt.plot(MW_OUT_3[0]['x_vector'][int(203/2),:],MW_OUT_3[0]['Kd_TW'][0][0][int(203/2),:]/MW_OUT_3[0]['Kd_EB'][0][0][int(203/2),:],label = 'MW L/25')
plt.xlim(-400,400)
plt.ylim(0.8,1.2)
plt.legend()