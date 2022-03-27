# -*- coding: utf-8 -*-
"""
mw_utilities_post
Created on Thu Feb 28 19:48:35 2019
@author: pbalitsk
those functions that perform post processing of MILDwave-NEMOH coupling
kd
total_wave
"""
import shutil
import numpy as np
import os
from xml.etree import ElementTree as et #Element Tree for modifying xml values

#need to split kd calculation for EB which is one CP which is multiple for each kk copupling loc
def kd_EB(MW_sim,MW_ini,MW_exe):

    # kd_aux = [] GAEL various coupling locations? wjhat aboput EB which is not the sum of multiple cpl locations??
    mwd=MW_ini['mw_dir']   # get the mw dir from the ii loop outside
    # create empty list inside 2-D list to fill with kd values    
#    for kk in range(len(MW_CPL[ii]['ph_LOC_x'])): #ii over various coup.lign lcoations TODO KK LOOP IS INSIDE FUNCTION        
    if MW_sim['regular']:          
        Tkd = MW_ini['T']
    else:
        Tkd = MW_ini['Tp']
    kd = [[] for jj in range(len(Tkd))]
    
    for jj in range(len(Tkd)):
        Tjj = Tkd[jj] 
        # for each mw directory for each period subdir and each run sudir inside (3 dir levels )
        #-----------------------------------------------------------------------------
        #                              3.1 Kd Calculation for EB
        #-----------------------------------------------------------------------------
        #First define and auxiliar value to save the eta loaded in each time step
        #Sort data files
        if MW_sim['regular']:       
            dir_T = os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj))
        else:
            dir_T = os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tjj))
            
        tree = et.parse(os.path.join(dir_T,"MILDwave.xml"))
        #-----------------------------------------------------------------------------
        #                               3.2 OBTAIN 'eta*' file name from first frequency
        #-----------------------------------------------------------------------------
        #obtain data of the wave gauge
        dx = np.float(tree.find(".//dx").text)*MW_ini['dx_bin_step']
        dy = np.float(tree.find(".//dy").text)*MW_ini['dy_bin_step']
        '''
        NWG = 1
        WG_y = [39.2/2.0]#,5.0/dy,5.0/dy,5.0/dy,5.0/dy]
        WG_y = [ '%.0f' % elem for elem in WG_y ]
        WG_x = [49.6/2.0]#,10.0/dx,15.0/dx,20.0/dx,25.0/dx]
        WG_x = [ '%.0f' % elem for elem in WG_x ]
        '''
        time_length = np.float(tree.find(".//etaOutput/end").text) - np.float(tree.find(".//etaOutput/start").text) - 30.0
        dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
        time_length = int(time_length/dt)*dt

        #WG = np.zeros((time_length,NWG))

        kd_size = MW_ini['kd_y'] * MW_ini['kd_x']
        eta = np.load(os.path.join(dir_T,'EB','data','eta.npy'))
        eta = [eta[mm:mm+kd_size] for mm in range(0, len(eta),kd_size)]
        eta = [np.reshape(eta_,( MW_ini['kd_y'], MW_ini['kd_x'])) for eta_ in eta]
        kd_ = [kd1**2 for kd1 in eta]
        kd_ = sum(kd_)
        kd_=kd_/time_length *dt
 
        if MW_sim['regular']: # only one H per T or Hs per Tp
            kd_ = [np.sqrt(8*kd_) /  MW_ini['H'][jj]]
            kd[jj]=kd_ # append kd to this list for each mw dir and each period
        elif MW_sim['irregular']:
            kd_ = [4*np.sqrt(kd_) / MW_ini['Hs'][jj]]
            kd[jj]=kd_#append kd to this list for each mw dir and each peak period
            # GAEL why sum of kds over various dirs  kd_aux.append(kd)
    #kd = [for kd[kd>=3]=np.nan
        #-----------------------------------------------------------------------------
    
        #-----------------------------------------------------------------------------
        #                              5.1 SAVE GENERATED Kd
        #-----------------------------------------------------------------------------
    S = np.shape(kd[0]);#tupple object indexinf for python ? GAEL 
    #     GAEL y and then x??
    Nx_eff = S[2]#the new effective domain corresponds with the domain without sponge layers
         #the domain can further be reduced defining a new value of Nx = target cells in x-direction
    Ny_eff = S[1]#the new effective domain corresponds with the domain without sponge layers
    #the domain can further be reduced defining a new value of Ny = target cells in y-direction
    xwidth=round((Nx_eff-1)*dx,3) #Value of the effective domain (in meters)along the x-axis.
    # Now it will extract all the grid points minus the Sponge Layer
    ylength=round((Ny_eff-1)*dy,3)#value of the effective domain(in meters)along the y-axis
    xvector = np.arange(-xwidth/2,xwidth/2+0.0001,dx)#adding the 0.0001 creates a correct vector for dx = dy = even number
    yvector = np.arange(-ylength/2,ylength/2+0.00001,dy)#adding the 0.0001 creates a correct vector for dx = dy = even number
    [xvector,yvector]=np.meshgrid(xvector,yvector)#Creating vector matrix
    del eta
    return(kd,xvector,yvector)
    
def kd_CP(MW_sim,MW_ini,MW_CPL,MW_exe):
    # kd_aux = [] GAEL various coupling locations? wjhat aboput EB which is not the sum of multiple cpl locations??
    mwd=MW_ini['mw_dir']   # get the mw dir from the ii loop outside
    # create empty list inside 2-D list to fill with kd values    
#    for kk in range(len(MW_CPL[ii]['ph_LOC_x'])): #ii over various coup.lign lcoations TODO KK LOOP IS INSIDE FUNCTION        
    if MW_sim['regular']:          
        Tkd = MW_ini['T']
    else:
        Tkd = MW_ini['Tp']
        #create a nested 2-D array of CPL by T and coupling location
    kd=[[[] for kk in range(len(MW_CPL['ph_LOC_x']))] for jj in range(len(Tkd))]
    for jj in range(len(Tkd)):
        Tjj = Tkd[jj] 
        # for each mw directory for each period subdir and each run sudir inside (3 dir levels )
        #-----------------------------------------------------------------------------
        #                              3.1 Kd Calculation for EB
        #-----------------------------------------------------------------------------
        #First define and auxiliar value to save the eta loaded in each time step
        #Sort data files
        if MW_sim['regular']:       
            dir_T = os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj))
        else:
            dir_T = os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tjj))
            
        tree = et.parse(os.path.join(dir_T,"MILDwave.xml"))
        #-----------------------------------------------------------------------------
        #                               3.2 OBTAIN 'eta*' file name from first frequency
        #-----------------------------------------------------------------------------
        #obtain data of the wave gauge
        dx = np.float(tree.find(".//dx").text)*MW_ini['dx_bin_step']
        dy = np.float(tree.find(".//dy").text)*MW_ini['dy_bin_step']
        '''
        NWG = 1
        WG_y = [39.2/2.0]#,5.0/dy,5.0/dy,5.0/dy,5.0/dy]
        WG_y = [ '%.0f' % elem for elem in WG_y ]
        WG_x = [49.6/2.0]#,10.0/dx,15.0/dx,20.0/dx,25.0/dx]
        WG_x = [ '%.0f' % elem for elem in WG_x ]
        '''
        time_length = np.float(tree.find(".//etaOutput/end").text) - np.float(tree.find(".//etaOutput/start").text) - 30.0
        dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
        time_length = int(time_length/dt)*dt

        #WG = np.zeros((time_length,NWG))

        kd_size = MW_ini['kd_y'] * MW_ini['kd_x']
        for kk in range(len(MW_CPL['ph_LOC_x'])):
            eta = np.load(os.path.join(dir_T,'CP_'+'LOC_'+str(kk),'data','eta.npy'))
            eta = [eta[mm:mm+kd_size] for mm in range(0, len(eta),kd_size)]
            eta = [np.reshape(eta_,( MW_ini['kd_y'], MW_ini['kd_x'])) for eta_ in eta]
            kd_ = [kd1**2 for kd1 in eta]
            kd_ = sum(kd_)
            kd_=kd_/time_length *dt

            if MW_sim['regular']: # only one H per T or Hs per Tp
                kd_ = [np.sqrt(8*kd_) /  MW_ini['H'][jj]]
                kd[jj][kk]=kd_ # append kd to this list for each and each period and coupling location
            elif MW_sim['irregular']:
                kd_ = [4*np.sqrt(kd_) / MW_ini['Hs'][jj]]
                kd[jj][kk]=kd_#append kd to this list for each and each period and coupling location
                # GAEL why sum of kds over various dirs  kd_aux.append(kd)
        #kd = [for kd[kd>=3]=np.nan
            #-----------------------------------------------------------------------------

            #-----------------------------------------------------------------------------
            #                              5.1 SAVE GENERATED Kd
            #-----------------------------------------------------------------------------
    S = np.shape(kd[0][0]);#tupple object indexinf for python ? GAEL 
    #     GAEL y and then x??
    Nx_eff = S[2]#the new effective domain corresponds with the domain without sponge layers
         #the domain can further be reduced defining a new value of Nx = target cells in x-direction
    Ny_eff = S[1]#the new effective domain corresponds with the domain without sponge layers
    #the domain can further be reduced defining a new value of Ny = target cells in y-direction
    xwidth=round((Nx_eff-1)*dx,3) #Value of the effective domain (in meters)along the x-axis.
    # Now it will extract all the grid points minus the Sponge Layer
    ylength=round((Ny_eff-1)*dy,3)#value of the effective domain(in meters)along the y-axis
    xvector = np.arange(-xwidth/2,xwidth/2+0.0001,dx)#adding the 0.0001 creates a correct vector for dx = dy = even number
    yvector = np.arange(-ylength/2,ylength/2+0.00001,dy)#adding the 0.0001 creates a correct vector for dx = dy = even number
    [xvector,yvector]=np.meshgrid(xvector,yvector)#Creating vector matrix
    del eta
    return(kd,xvector,yvector)
    #TODO aadd location k loop
    
# This is  how to initilzie a 2-D Array
# kd = [[[]] * len(MW_ini['T']) for ii in range(len(MW_exe['mw_dir']))]
def total_wave(MW_sim,MW_ini,MW_exe,MW_CPL): # GAEL postpro?
    #-----------------------------------------------------------------------------
    #                              3.1 OBTAIN 'eta*' file name from first frequency
    #-----------------------------------------------------------------------------
    #data data_folder inside the MW run data_folder
    # kd_tw = [] various coupling locations?
    # create empty list inside 2-D list to fill with kd values            
    mwd=MW_ini['mw_dir']   # get the mw dir from the ii loop outside       
    if MW_sim['regular']:          
        Tkd = MW_ini['T']
        kd = [[] for jj in range(len(MW_ini['T']))]
    else: # irreg
        Tkd = MW_ini['Tp']
        kd = [[] for jj in range(len(MW_ini['Tp']))]
        
    for jj in range(len(Tkd)):
        
        if MW_sim['regular']: 
            eta_cp=[[] for kk in range(len(MW_CPL['ph_LOC_x']))] # empty list of coupled :Eta of length kk
            Tjj = MW_ini['T'][jj]
            dir_T = os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj))
            eta = np.load(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),'EB','data','eta.npy'))   
            #First define and auxiliar value to save the eta loaded in each time step
            tree = et.parse(os.path.join(dir_T,"MILDwave.xml"))
            #obtain data of the wave gauge
            dx = np.float(tree.find(".//dx").text)*MW_ini['dx_bin_step']
            dy = np.float(tree.find(".//dy").text)*MW_ini['dy_bin_step']
            '''
            NWG = 1
            WG_y = [39.2/2.0]#,5.0/dy,5.0/dy,5.0/dy,5.0/dy]
            WG_y = [ '%.0f' % elem for elem in WG_y ]
            WG_x = [49.6/2.0]#,10.0/dx,15.0/dx,20.0/dx,25.0/dx]
            WG_x = [ '%.0f' % elem for elem in WG_x ]
            '''
            kd_size = MW_ini['kd_y'] * MW_ini['kd_x']
            time_length = np.float(tree.find(".//etaOutput/end").text) - np.float(tree.find(".//etaOutput/start").text) - 30.0
            dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
            #WG = np.zeros((time_length,NWG))
            #we load incident wave and save it in the total wave variable
            for kk in range(len(MW_CPL['ph_LOC_x'])):
                tw_subdir = 'TW_'+'LOC_'+str(kk)
                CP_subdir = 'CP_'+'LOC_'+str(kk)
                if os.path.exists(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir,'data')):
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir,'data'),os.path.join(mwd,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
    
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir),os.path.join(mwd,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
    
                    os.mkdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir))
                    os.mkdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir,'data'))
                    #If folder name does not exist it creates the folder name and data file
                elif os.path.isdir(os.path.join(mwd,tw_subdir)):
    
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(mwd,tw_subdir),os.path.join(mwd,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
    
                    os.mkdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir))
                    os.mkdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir,'data'))
    
                else:
                    if not(os.path.isdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir))): #in case wherePL direcory already exists but data subdir does not
                            os.mkdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir))
                    os.mkdir(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir,'data'))

                eta_cp[kk] = np.load(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),CP_subdir,'data','eta.npy'))#we load perturbed wave
            #-----------------------------------------------------------------------------
            #                              3.2 Kd Calculation Reg
            #-----------------------------------------------------------------------------    
            eta_p=[] #empty list for perturbed wave            
            eta_p = eta + sum(eta_cp)#we sum the incident and the perturbed wave
            del eta_cp#we delete the perturbed wave to save memory
            np.save(os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj),tw_subdir,'data','eta'),eta_p)
            eta_p = [eta_p[mm:mm+kd_size] for mm in range(0, len(eta_p),kd_size)]
            eta_p = [np.reshape(eta_,( MW_ini['kd_y'], MW_ini['kd_x'])) for eta_ in eta_p]
            kd_p = [kd1**2 for kd1 in eta_p] #udnerscore is temp variables
            kd_p = sum(kd_p)
            kd_p= kd_p/time_length *dt
            
        else: #irregular    
            #we load incident wave only once for each total wave run and save it in the total wave variable
            eta_cp=[[] for kk in range(len(MW_CPL['ph_LOC_x']))]# empty list of coupled :Eta of length kk
            Tpjj = MW_ini['Tp'][jj]
            eta = np.load(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),'EB','data','eta.npy'))   
            dir_T = os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj))       
            tree = et.parse(os.path.join(dir_T,"MILDwave.xml"))
            #obtain data of the wave gauge
            dx = np.float(tree.find(".//dx").text)*MW_ini['dx_bin_step']
            dy = np.float(tree.find(".//dy").text)*MW_ini['dy_bin_step']
            kd_size = MW_ini['kd_y'] * MW_ini['kd_x']
            time_length = np.float(tree.find(".//etaOutput/end").text) - np.float(tree.find(".//etaOutput/start").text) - 30.0
            dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
            #WG = np.zeros((time_length,NWG))
            #loop over the various coupling locations to get the :eta and sum them up with the EB :eta
            for kk in range(len(MW_CPL['ph_LOC_x'])):
                tw_subdir = 'TW_'+'LOC_'+str(kk)
                CP_subdir = 'CP_'+'LOC_'+str(kk)
                if os.path.exists(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir,'data')):
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir,'data'),os.path.join(mwd,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
    
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir),os.path.join(mwd,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
    
                    os.mkdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir))
                    os.mkdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir,'data'))
                    #If folder name does not exist it creates the folder name and data file
                    
                elif os.path.isdir(os.path.join(mwd,tw_subdir)):
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(mwd,tw_subdir),os.path.join(mwd,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
    
                    os.mkdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir))
                    os.mkdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir,'data'))
    
                else:
                    if not(os.path.isdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir))): #in case wherePL direcory already exists but data subdir does not
                            os.mkdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir))
                    os.mkdir(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir,'data'))
            #First define and auxiliar value to save the eta loaded in each time step
                eta_cp[kk] = np.load(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),CP_subdir,'data','eta.npy'))#we load perturbed wave
            #-----------------------------------------------------------------------------
            #                              3.2 Kd Calculation Irreg
            #-----------------------------------------------------------------------------
            eta_p=[] #empty list for perturbed wave
            eta_p = eta + sum(eta_cp)#we sum the incident and the perturbed wave for each location kk
            del eta_cp#we delete the perturbed wave to save memory
            np.save(os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj),tw_subdir,'data','eta'),eta_p)
            eta_p = [eta_p[mm:mm+kd_size] for mm in range(0, len(eta_p),kd_size)]
            eta_p = [np.reshape(eta_,( MW_ini['kd_y'], MW_ini['kd_x'])) for eta_ in eta_p]
            kd_p = [kd1**2 for kd1 in eta_p] #udnerscore is temp variables
            kd_p = sum(kd_p)
            kd_p = kd_p/time_length *dt

        if MW_sim['regular']:  # only one H per T or Hs per Tp
            kd_p = [np.sqrt(8*kd_p) / MW_ini['H'][jj]]
            kd[jj]=kd_p #for each coupling locatuion kd to the sumthis list for each mw dir and each period
        elif  MW_sim['irregular']: #irreg    #NB chenge .append to = so that kd is not nested   
            kd_p = [4*np.sqrt(kd_p) / MW_ini['Hs'][jj]]
            kd[jj]=kd_p#append kd to this list for each mw dir and each peak period
            #for kd[kd>=3]=np.nan TODO REMOVE HIGH KD VALUES
        #kd_aux.append(kd)
            #-----------------------------------------------------------------------------
            #                              5.1 SAVE GENERATED Kd
            #-----------------------------------------------------------------------------
    S = np.shape(eta_p)#tupple object indexinf for python
    Nx_eff = S[2]#the new effective domain corresponds with the domain without sponge layers
             #the domain can further be reduced defining a new value of Nx = target cells in x-direction
    Ny_eff = S[1]#the new effective domain corresponds with the domain without sponge layers
             #the domain can further be reduced defining a new value of Ny = target cells in y-direction
    xwidth=round(Nx_eff*dx,3) #Value of the effective domain (in meters)along the x-axis. Now it will extract all the grid points minus the Sponge Layer
    ylength=round(Ny_eff*dy,3)#value of the effective domain(in meters)along the y-axis
    xvector = np.arange(-xwidth/2,xwidth/2,dx)
    yvector = np.arange(-ylength/2,ylength/2,dy)
    [xvector,yvector]=np.meshgrid(xvector,yvector)#Creating vector matrixe
    
    return(kd)



# function which calls MW executable# -*- coding: utf-8 -*-