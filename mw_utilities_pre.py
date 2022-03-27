
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 14:12:31 2018
@author: gveraofe
mod Phil Dog  Feb 28 18:20:25 2019
mw_utilitiy functions used for preparing for the MILDwave runs i.e. prepro
create_input
reg_init
irr_irr_init
and helper Functions
"""
import numpy as np
import os
import spectral_functions as sf
from xml.etree import ElementTree as et #Element Tree for modifying xml values
import math as math
import shutil
#Creates all the input variables needed for running MILDwave once the type of case and input wave conditions are defined
#save oin MW_ini dict
def create_input (MW_sim,MW_ini): #NB variables that do not depend on dd will be the same over all the simulation mw_dirs
    regCheck = MW_sim['regular']
    irrCheck = MW_sim['irregular']
    dirCheck = MW_sim['directional']
    if regCheck:
        MW_ini['f']=[]
        MW_ini['w']=[]
        MW_ini['amp']=[]
        for jj in range(len(MW_ini['T'])):
#                    MW_ini['f'].append([1/T for T in Ti])
            MW_ini['f'].append(1/MW_ini['T'][jj])
        for ii in range(len( MW_ini['f'])):
            MW_ini['w'].append(2*np.pi*MW_ini['f'][jj])
        for ii in range(len(MW_ini['H'])):
            MW_ini['amp'].append(MW_ini['H'][0]/2)

    elif irrCheck:
        MW_ini['f']=[[] for jj in range(len(MW_ini['Tp']))]
        MW_ini['T']=[[] for jj in range(len(MW_ini['Tp']))]
        MW_ini['w']=[[] for jj in range(len(MW_ini['Tp']))]
        #check for list errors here too ^
        MW_ini['fp'] = [1/Tp for Tp in MW_ini['Tp']]# Peak frequency in Hz
        MW_ini['wp'] = [2*np.pi*fp for fp in MW_ini['fp']] 
        for jj in range(len(MW_ini['Tp'])):
#            fp = MW_ini['fp'][jj]           
            MW_ini['f'][jj]=(list(reversed(np.round(np.linspace(MW_ini['fini'][jj],MW_ini['fend'][jj],MW_ini['Nf'],endpoint = True),3)))) #Linespacing the frequency maybe is better to do it with the period?? 
#                MW_ini['f'][jj]=(list(reversed((np.round(np.linspace(MW_ini['fini'][jj],MW_ini['fend'][jj],MW_ini['Nf'],endpoint = True),3)*fp)))) #Linespacing the frequency maybe is better to do it with the period?? 
                          #Limits should be adjusted depending on the peak period
#            for fp in MW_ini['f']:
            MW_ini['T'][jj]=([round((1/f),3) for f in MW_ini['f'][jj]  ])
            MW_ini['w'][jj]=([2*np.pi*f for f in MW_ini['f'][jj]  ])
            
        if MW_ini['spectra'] == 'JS':
            sf.JS(MW_ini)
        elif MW_ini['spectra'] =='PM':
            sf.PM(MW_ini)
            
        if dirCheck:
            #filename = 'dir_nemoh_input.txt'#Directiosn used in the MILDwave SC simulation
            #for ii in MW_ini['dirfile']:
                #deg_MW = np.loadtxt(os.path.join(MW_ini['dirfile'][ii],filename),skiprows=0)
                #deg_MW = deg_MW[:,1]
                #MW_ini['deg'] =list(reversed(deg_MW))   #TODO IMPLEMENT FOR DIFFERENT DIRECTIONAL WAVES CASES
           deg_aux = sf.DSM(MW_ini)
        #CORRECTION OF THE ANGLES ACORDING TO MILDWAVE CORRECTION IN THE WAVE GENERATION LINE IN RELATION WITH GRID SIZE AND WAVE GENERATION LINE SIZE
           for ii in range(len(MW_ini['deg_main'])):
           #TODO zips
               (L_T) = [wavlen(H,T,MW_ini['depth'][ii]) for H,T in zip(MW_ini['H'][ii],MW_ini['T'][ii]) ]#L of each T
               k_T = [2*np.pi/L for L in L_T]#angular wave number
               tree = et.parse(os.path.join(MW_ini['mw_dir'],"MILDwave.xml"))
               jm1 = int(tree.find(".//Ny").text)
               #dx = 0.08#float(tree.find(".//dx").text)
               dy = float(tree.find(".//dy").text)
               sinRw = [math.sin(Theta) for Theta in deg_aux[ii]]
               rval = [np.round(k*sin*(jm1 - 1)*dy/(2*np.pi))*2*np.pi/((jm1-1)*dy*k) for sin,k in zip(sinRw,k_T)]
               deg = np.zeros(MW_ini['Nf'])
               #TODO INPUT ERROR MESSAGE IF THE PERIODIC LENGHT  IS NOT ENOUGH THE PROGRAM GIVES THE MAIN DIRECTION
               #     PERIODIC LENGHT HAS TO BE 3L AT LEAST
               for jj in range(len(rval)):
                  if (rval[jj] > 1.0):
                      deg[jj] = math.asin(rval[jj] - 2*np.pi/((jm1-1)*dy*k_T[jj]))*180/np.pi + MW_ini['deg_main'][ii] #INPUT in MW AND NEMOH IS DEGREES and WE HAVE TO SHIFT FOR THE MAIN DIRECTION
                  elif (rval[jj] < -1.0):
                      deg[jj] = math.asin(rval[jj]+2*np.pi/((jm1-1)*dy*k_T[jj]))*180/np.pi + MW_ini['deg_main'][ii]   #INPUT in MW AND NEMOH IS DEGREES and WE HAVE TO SHIFT FOR THE MAIN DIRECTION
                  else:
                      deg[jj] = math.asin(rval[jj])*180/np.pi + MW_ini['deg_main'][ii]     #INPUT in MW AND NEMOH IS DEGREES and WE HAVE TO SHIFT FOR THE MAIN DIRECTION
               MW_ini['deg'].append(deg)
        else:
            MW_ini['deg'] = [deg*MW_ini['Nf'] for deg in MW_ini['deg']]
            #            TODO this needs to loop over list of dir_MW
            # extract values from xml file for coupling params            
    return(MW_ini)

# intitiliatize regular wave run
# overwrite MILDwave.xml file
def reg_init(MW_ini,MW_exe):
    #create MILDwave.xml file in each MILDwave run directory
    # for mw_dir,deg,depth in zip(MW_exe['mw_dir'],MW_ini['deg'],MW_ini['depth']):
    mwd=MW_ini['mw_dir']
    # T is an list of lists!
    deg,depth = MW_ini['deg'],MW_ini['depth'][0]
    # for TT in MW_ini['T']:
    # dir_T = [os.path.join(mw_dir,'EB','T_'+'{:06.3f}'.format(Tjj)) for Tjj in T ]
    # TODO implement for various periods
    # set up for future use with command line ofptions for HH in MW_ini['H']
            # for deg in NEM_ini['deg']:
                # for depth in NEM_ini['depth']:
        #for now will have to deal with single values
        # change so that the .xml file is in the mw_dir and not the EB subduir   
    # these variables do not depend on T and only loaded onze
    tree = et.parse(os.path.join(mwd,"MILDwave.xml"))# MILDwave.ini input file?
#        T =  float(tree.find(".//wavePeriod").text)
    dx = float(tree.find(".//dx").text)
    dt_ini = float(tree.find(".//delt").text)
    eta_Increment_ini = int(float(tree.find(".//etaOutput/increment").text))
    dt_eta = dt_ini*eta_Increment_ini
    width = (int(tree.find(".//Nx").text)-2*int(tree.find(".//ixsL").text))*float(tree.find(".//dx").text)
    
    # create SUBdirectory for each period in T in MW_ini - disregard the T in the mw.xml   
    for jj in range(len(MW_ini['T'])):
        Tjj = MW_ini['T'][jj]
        Hjj = MW_ini['H'][jj]
        # create the subdirectory for each peak period isnot already there
        dir_T = os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj))
        if not(os.path.isdir(os.path.join(mwd,dir_T))):
            os.mkdir(os.path.join(mwd,dir_T))
        #copy existing MILDwave.xml file to subfodler for each period
        shutil.copyfile(os.path.join(mwd,"MILDwave.xml"),os.path.join(mwd,dir_T,"MILDwave.xml"))
        shutil.copyfile(os.path.join(mwd,"depth.dat"),os.path.join(mwd,dir_T,"depth.dat"))
        shutil.copyfile(os.path.join(mwd,"obstacle.dat"),os.path.join(mwd,dir_T,"obstacle.dat"))
                        
        L_T = 9.81*Tjj**2/(2*np.pi)#new dt for calculating the different periods based on the Courant-Friedrichs-Lewy criterion
        k_T = 2*np.pi/L_T#wave number for each T
        w_T = 2*np.pi/Tjj#angular frequency for each T
        c_w = w_T/k_T #maximum phase velocity
        dt_C = float(truncate(dx/c_w,2))#Time step suggested in base of the Courant-Firedrichs-Lewy Criterion use as reference to calcualte the new dt for the irregular wave calculation
        dt_C = dt_C - 0.05#The idea of substraction 0.05 is that it is never a value close to the dt_C which can cause MILDwave to rash
        etaIncrement = (first(range(1,10), lambda ii: dt_eta/ii < dt_C))
        dt = dt_eta/etaIncrement #new dt for calculating the different periods based on the Courant-Friedrichs-Lewy criterion
    
        c = L_T/Tjj#Celerity for each Tjj
        cg_T = 0.5 * (1+4*np.pi*depth/L_T/np.sinh(4*np.pi*depth/L_T))*c #Group Celerity of each T
        wave_time_start = width/cg_T+MW_exe['time_ramp'] #commence calculating Kd when the first fully develop wave arrives
        wave_time_end = wave_time_start+MW_exe['time_Kd']
        extra_time = 10 #extra time for the simulation not computed in seconds
    
        etaPhi_Start = np.round(np.round(wave_time_start/Tjj,0)*Tjj/dt,0)*dt
        #Round up the initial time for surface elevation calculation so it is a multiple of the dt
        etaPhi_End = etaPhi_Start + MW_exe['time_Kd']
        #Parse the MILDwave.xml that has ben copied in the subfolder for each period
        tree = et.parse(os.path.join(mwd,dir_T,"MILDwave.xml"))# MILDwave.ini input file?

        tree.find(".//twfin").text = '%.0f'%(wave_time_end+extra_time)
        tree.find(".//delt").text = '%.3f'%(dt)
        tree.find(".//stopWGEN").text = '%.0f'%(wave_time_end+extra_time)
        tree.find(".//etaOutput/start").text = '%.2f'%etaPhi_Start
        tree.find(".//etaOutput/increment").text = '%.2f'%etaIncrement
        tree.find(".//etaOutput/end").text = '%.2f'%etaPhi_End
        tree.find(".//phiOutput/start").text = '%.2f'%etaPhi_Start
        tree.find(".//phiOutput/increment").text = '%.2f'%etaIncrement
        tree.find(".//phiOutput/end").text = '%.2f'%(etaPhi_Start+2*dt)
        tree.find(".//wavePeriod").text = '%.6f'%Tjj
        tree.find(".//waveHeight").text = '%.6f'%Hjj
        
        # save overwritten xml file to root dirctory
        tree.write(os.path.join(mwd,dir_T,"MILDwave.xml"))
        #create EB doler if not exist already
        if not(os.path.isdir(os.path.join(mwd,dir_T,'EB'))):
            os.mkdir(os.path.join(mwd,dir_T,'EB'))
        # and copy to EB subdirectory
        shutil.copyfile(os.path.join(mwd,dir_T,"MILDwave.xml"),os.path.join(mwd,dir_T,'EB',"MILDwave.xml"))
        shutil.copyfile(os.path.join(mwd,dir_T,"depth.dat"),os.path.join(mwd,dir_T,'EB',"depth.dat"))
        shutil.copyfile(os.path.join(mwd,dir_T,"obstacle.dat"),os.path.join(mwd,dir_T,'EB',"obstacle.dat"))
            
# intitiliatize irregular wave run
# overwrite MILDwave.xml file

def irr_init(MW_ini,MW_exe):#creates all the MW_exe files for running a irregular empty simulation in MILDwave from the Tp MILDwave.xml file
 #TODO add EB directory creation and movement of "DUMMY" Xml file and depth and obstacle .dat to subdirectory
    # for mw_dir,T,H,deg,depth in zip(MW_exe['mw_dir'],MW_ini['T'],MW_ini['H'],MW_ini['deg'],MW_ini['depth']):
        mwd,depth = MW_ini['mw_dir'],MW_ini['depth'][0]
        
        tree = et.parse(os.path.join(mwd,"MILDwave.xml"))# MILDwave.ini input file?
        #        T =  float(tree.find(".//wavePeriod").text)
        dx = float(tree.find(".//dx").text)
        dt_ini = float(tree.find(".//delt").text)
        eta_Increment_ini = int(float(tree.find(".//etaOutput/increment").text))
        dt_eta = dt_ini*eta_Increment_ini
        width = (int(tree.find(".//Nx").text)-2*int(tree.find(".//ixsL").text))*float(tree.find(".//dx").text)
        
        for jj in range(len(MW_ini['Tp'])):
            Tpjj = MW_ini['Tp'][jj]
            Tn = MW_ini['T'][jj]
            H = MW_ini['H'][jj]
            deg = MW_ini['deg'][jj]
            # for HH in MW_ini['H']
            # for deg in NEM_ini['deg']:
            # for depth in NEM_ini['depth']:
            # works for one direction and one depth only
            # create the subdirectory for each peak period isnot already there
            dir_Tp = os.path.join(mwd,'Tp_'+'{:05.2f}'.format(Tpjj))
            if not(os.path.isdir(os.path.join(mwd,dir_Tp))):
                os.mkdir(os.path.join(mwd,dir_Tp))
            #copy existing MILDwave.xml file to subfodler for each period
            shutil.copyfile(os.path.join(mwd,"MILDwave.xml"),os.path.join(mwd,dir_Tp,"MILDwave.xml"))
            shutil.copyfile(os.path.join(mwd,"depth.dat"),os.path.join(mwd,dir_Tp,"depth.dat"))
            shutil.copyfile(os.path.join(mwd,"obstacle.dat"),os.path.join(mwd,dir_Tp,"obstacle.dat"))
            
            L_Tp = 9.81*Tpjj**2/(2*np.pi)#new dt for calculating the different periods based on the Courant-Friedrichs-Lewy criterion
            k_Tp = 2*np.pi/L_Tp#wave number for each T
            w_Tp = 2*np.pi/Tpjj#angular frequency for each T
            c_wp = w_Tp/k_Tp #maximum phase velocity
            dt_C = float(truncate(dx/c_wp,2))#Time step suggested in base of the Courant-Firedrichs-Lewy Criterion use as reference to calcualte the new dt for the irregular wave calculation
            dt_C = dt_C - 0.05#The idea of substraction 0.05 is that it is never a value close to the dt_C which can cause MILDwave to rash
            etaIncrement = (first(range(1,10), lambda ii: dt_eta/ii < dt_C))
            dt = dt_eta/etaIncrement #new dt for calculating the different periods based on the Courant-Friedrichs-Lewy criterion
        
            cp = L_Tp/Tpjj#Celerity for each Tpjj
            cg_Tp = 0.5 * (1+4*np.pi*depth/L_Tp/np.sinh(4*np.pi*depth/L_Tp))*cp #Group Celerity of each T
            wave_time_start = width/cg_Tp+MW_exe['time_ramp'] #commence calculating Kd when the first fully develop wave arrives
            wave_time_end = wave_time_start+MW_exe['time_Kd']
            extra_time = 10 #extra time for the simulation not computed in seconds
        
            etaPhi_Start = np.round(np.round(wave_time_start/Tpjj,0)*Tpjj/dt,0)*dt
            #Round up the initial time for surface elevation calculation so it is a multiple of the dt
            etaPhi_End = etaPhi_Start + MW_exe['time_Kd']
            #Parse the MILDwave.xml that has ben copied in the subfolder for each period
            tree1 = et.parse(os.path.join(mwd,dir_Tp,"MILDwave.xml"))# MILDwave.ini input file?
    
            tree1.find(".//twfin").text = '%.0f'%(wave_time_end+extra_time)
            tree1.find(".//delt").text = '%.3f'%(dt)
            tree1.find(".//stopWGEN").text = '%.0f'%(wave_time_end+extra_time)
            tree1.find(".//etaOutput/start").text = '%.2f'%etaPhi_Start
            tree1.find(".//etaOutput/increment").text = '%.2f'%etaIncrement
            tree1.find(".//etaOutput/end").text = '%.2f'%etaPhi_End
            tree1.find(".//phiOutput/start").text = '%.2f'%etaPhi_Start
            tree1.find(".//phiOutput/increment").text = '%.2f'%etaIncrement
            tree1.find(".//phiOutput/end").text = '%.2f'%(etaPhi_Start+2*dt)
        
            # save overwritten xml file to root dirctory
            tree1.write(os.path.join(mwd,dir_Tp,"MILDwave.xml"))
            #create EB doler if not exist already
            if not(os.path.isdir(os.path.join(mwd,dir_Tp,'EB'))):
                os.mkdir(os.path.join(mwd,dir_Tp,'EB'))
            
            tree1 = et.parse(os.path.join(dir_Tp,"MILDwave.xml"))#DUMMY MILDwave.ini input file not used in final calculation
            
            #create subdirectory for each irr wave period in the simulation
            dir_Tn = [os.path.join(dir_Tp,'EB','T_'+'{:06.3f}'.format(Tnn)) for Tnn in Tn]
            
            L_T = [9.81*Tnn**2/(2*np.pi) for Tnn in Tn]#L of each T in deep water as it is in the MILDwave source code SettingsForm1.ccp
            k_T = [2*np.pi/Lii for Lii in L_T]#wave number for each T
            w_T = [2*np.pi/Tnn for Tnn in Tn]#angular frequency for each T
            c_w = [w/k for w,k in zip(w_T,k_T)]#maximum phase velocity
            dt_C = [float(truncate(dx/c,2)) for c in c_w]#Time step suggested in base of the Courant-Firedrichs-Lewy Criterion use as reference to calcualte the new dt for the irregular wave calculation
            dt_C = [dt - 0.05 for dt in dt_C]#The idea of substraction 0.05 is that it is never a value close to the dt_C which can cause MILDwave to rash
            eta_Increment = [(first(range(1,10), lambda ii: dt_eta/ii < dt_Ci)) for dt_Ci in dt_C]
            #GAEL at what point isd the dteta seta and how is criteria determined 
            dt_irr = [dt_eta/eta_Inc for eta_Inc in eta_Increment]#new dt for calculating the different periods based on the Courant-Friedrichs-Lewy criterion

            c = [ L/Tnn for L,Tnn in zip(L_T,Tn)]#Celerity for each T
            cg_T = [0.5 * (1+4*np.pi*depth/L/np.sinh(4*np.pi*depth/L))*cii for L,cii in zip(L_T,c)]#Group Celerity of each T
            wave_time_start = [width/cg+MW_exe['time_ramp'] for cg in cg_T ] #strat calculating Kd when the first fully develop wave arrives
            wave_time_end = [time_value+MW_exe['time_Kd'] for time_value in wave_time_start]
            extra_time = 10#extra time for the simulation not computed in seconds
            # loop over the subperiods T in the irr wave run
            for dir_Tnn,Tnn,Hnn,beta_0,etaPhiOutStart,etaPhiOutEnd,etaIncrement,dt in zip(dir_Tn,Tn,H,deg,wave_time_start,wave_time_end,eta_Increment,dt_irr):
               
                tree_sub = et.parse(os.path.join(dir_Tp,"MILDwave.xml"))
                        
                if os.path.exists(os.path.join(dir_Tnn,'data')):
                
                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(dir_Tnn,'data'),os.path.join(dir_Tnn,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(dir_Tnn,del_name),ignore_errors=True)

                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(dir_Tnn),os.path.join(dir_Tp,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(dir_Tp,del_name),ignore_errors=True)

                    os.mkdir(os.path.join(dir_Tnn))
                    
                elif os.path.exists(os.path.join(dir_Tnn)):

                    del_name = 'id1'#Auxiliar name for deleting existing folder
                    os.rename(os.path.join(dir_Tnn),os.path.join(dir_Tp,del_name))#renames existing folder to delete and avoid system permision issues
                    shutil.rmtree(os.path.join(dir_Tp,del_name),ignore_errors=True)
                    os.mkdir(os.path.join(dir_Tnn))
                else:
                    os.mkdir(os.path.join(dir_Tnn))

                tree_sub.find(".//wavePeriod").text = '%.6f'%Tnn
                tree_sub.find(".//waveHeight").text = '%.6f'%Hnn
                tree_sub.find(".//beta_0").text = '%.6f'%beta_0

                etaPhi_Start = np.round(np.round(etaPhiOutStart/Tnn,0)*Tnn/dt,0)*dt
                #Round up the initial time for surface elevation calculation so it is a multiple of the dt
                etaPhi_End = etaPhi_Start + MW_exe['time_Kd']

                tree_sub.find(".//twfin").text = '%.0f'%(etaPhiOutEnd+extra_time)
                tree_sub.find(".//delt").text = '%.3f'%(dt)
                tree_sub.find(".//stopWGEN").text = '%.0f'%(etaPhiOutEnd+extra_time)
                tree_sub.find(".//etaOutput/start").text = '%.2f'%etaPhi_Start
                tree_sub.find(".//etaOutput/increment").text = '%.2f'%etaIncrement
                tree_sub.find(".//etaOutput/end").text = '%.2f'%etaPhi_End
                tree_sub.find(".//phiOutput/start").text = '%.2f'%etaPhi_Start
                tree_sub.find(".//phiOutput/increment").text = '%.2f'%etaIncrement
                tree_sub.find(".//phiOutput/end").text = '%.2f'%(etaPhi_Start+2*dt)
                tree_sub.write(os.path.join(dir_Tnn,"MILDwave.xml"))
                # TODO check for EB here
                shutil.copyfile(os.path.join(dir_Tp,"depth.dat"),os.path.join(dir_Tnn,"depth.dat"))
                shutil.copyfile(os.path.join(dir_Tp,"obstacle.dat"),os.path.join(dir_Tnn,"obstacle.dat"))
                #if WGCheck:#TO DO IMPLEMENT WAVE GUAGES
                    #shutil.copyfile(os.path.join(dir_T,"Pos_WG.xml"),os.path.join(dir_T,"Pos_WG.xml"))


# least common criterion condition
            
def first(the_iterable, condition = lambda x: True):
    for ii in the_iterable:
        if condition(ii):
            return ii

def wavlen(H,T,d):
    g = 9.81
    dpi = 2*np.pi

    L0 = g*T**2/dpi
    L1 = g*T**2/dpi*np.tanh(dpi*d/L0)

    while (np.abs(L1-L0) > 0.001):
        L0 = L1
        L1 = g*T**2/dpi*np.tanh(dpi*d/L0)
    return L1

def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])