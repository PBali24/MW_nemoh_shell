# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 19:48:35 2019

@author: pbalitsk
functions Necessary for runnign MILDwave and calling the MW_CLI
run_mw set up the MILDwave runs and shuffle around the outputs
execute_mw -call subprocess to run mw_cli.exe
"""
import numpy as np
import os
from xml.etree import ElementTree as et #Element Tree for modifying xml values
import subprocess
import shutil
import time
from subprocess import Popen

def run_mw(MW_sim,MW_ini,mw_dir,Tjj,MW_exe,run_type,Tnn):

    dir_exe = MW_exe['mw_cli']
    dx_bin_step=MW_ini['dx_bin_step'] # command linie input of the :eta.bin output step (~file size)
    if MW_sim['regular']:
        dir_data = 'data'
        T_dir= os.path.join(mw_dir,'T_'+'{:05.2f}'.format(Tjj),run_type)

        #call to the functions which runs the MILDwave CLI
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF            
        execute_mw(dir_exe,T_dir,MW_sim['regular'],dx_bin_step)
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF            
        tree = et.parse(os.path.join(T_dir,"MILDwave.xml"))
        time_length = float(truncate(np.float(tree.find(".//etaOutput/end").text),3)) - float(truncate(np.float(tree.find(".//etaOutput/start").text),3)) - 30.0#remove thirty seconds in case the number of etas is not the same for each simulation
        dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
        time_length = int(time_length/dt)*dt
        npoints = int((time_length/dt) * MW_ini['kd_x']  * MW_ini['kd_y'] )
        with open(os.path.join(T_dir,dir_data,'eta.bin')) as feta:
            eta = np.fromfile(feta,dtype = np.float32,count = npoints,sep='')
        os.remove(os.path.join(T_dir,dir_data,'eta.bin'))
        np.save(os.path.join(T_dir,dir_data,'eta'),eta)
        del eta
    else:  # Irregular
        n = MW_exe['n_cores']
        dir_data = 'data'
        Tp_dir = os.path.join(mw_dir,'Tp_'+'{:05.2f}'.format(Tjj))
        if os.path.exists(os.path.join(Tp_dir,run_type,dir_data)):
            del_name = 'id1'#Auxiliar name for deleting existing folder
            os.rename(os.path.join(Tp_dir,run_type,dir_data),os.path.join(Tp_dir,run_type,del_name))#renames existing folder to delete and avoid system permision issues
            shutil.rmtree(os.path.join(Tp_dir,run_type,del_name))
            os.mkdir(os.path.join(Tp_dir,run_type,dir_data))
            #If folder name does not exist it creates the folder name and data file
        else:
            os.mkdir(os.path.join(Tp_dir,run_type,dir_data))
            #TODO check jj and ii
            
        dir_Tnn = [os.path.join(Tp_dir,run_type,'T_'+'{:06.3f}'.format(Tn)) for Tn in Tnn ]

        #f includes n sublist with the total number of simulations to run each time in each core
        f_runs = [dir_Tnn[i * n:(i + 1) * n] for i in range((len(dir_Tnn) + n - 1) // n )]

        tree = et.parse(os.path.join(Tp_dir, "MILDwave.xml"))
        time_length = float(truncate(np.float(tree.find(".//etaOutput/end").text),3)) - float(truncate(np.float(tree.find(".//etaOutput/start").text),3)) - 30.0#remove thirty seconds in case the number of etas is not the same for each simulation
        dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
        time_length = int(time_length/dt)*dt
        npoints = int((time_length/dt) * MW_ini['kd_x'] * MW_ini['kd_y'])

        for subdir_mw in f_runs:
            print(f'MW CLI calculation started: '+str(time.asctime()))
            # calls executing function in chunks
            # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF            

            execute_mw(dir_exe,subdir_mw,MW_sim['regular'],dx_bin_step)#runs n numbers of frequencies and waits for all of them to finish
            # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF            

            print(f'MW CLI calculation finished: '+str(time.asctime()))

            print(f'MW loading eta.bin started: '+str(time.asctime()))

            if os.path.exists(os.path.join(Tp_dir,run_type,dir_data,'eta.npy')):
                eta_irr = np.load(os.path.join(Tp_dir,run_type,dir_data,'eta.npy'))
            else:
                eta_irr = 0

            for dir_name in subdir_mw:
                with open(os.path.join(dir_name,dir_data,'eta.bin')) as feta:
                    eta_irr += np.fromfile(feta,dtype = np.float32,count = npoints,sep='')
                os.remove(os.path.join(dir_name,dir_data,'eta.bin'))
            np.save(os.path.join(Tp_dir,run_type,dir_data,'eta'),eta_irr)
            del eta_irr
            print(f'MW loading eta.bin sfinished: '+str(time.asctime()))
    return None
    
# funcion to execute MILDwave CLI.exe bia subprocess     
def execute_mw(dir_exe,T_dir,regular,dx_bin_step):
    dir_data = 'data'
    if regular:
        MW_settings_file = os.path.join(T_dir,'MILDwave.xml')
        if os.path.exists(os.path.join(T_dir,dir_data)):
            del_name = 'id1'#Auxiliar name for deleting existing folder
            os.rename(os.path.join(T_dir,dir_data),os.path.join(T_dir,del_name))#renames existing folder to delete and avoid system permision issues
            shutil.rmtree(os.path.join(T_dir,del_name))
            os.mkdir(os.path.join(T_dir,dir_data))
            #If folder name does not exist it creates the folder name and data file
        else:
            os.mkdir(os.path.join(T_dir,dir_data))
        subprocess.check_call([dir_exe,'-X', str(dx_bin_step), MW_settings_file])
    else:
        for folder in T_dir:
            if os.path.exists(os.path.join(folder,dir_data)):
                del_name = 'id1'#Auxiliar name for deleting existing folder
                os.rename(os.path.join(folder,dir_data),os.path.join(folder,del_name))#renames existing folder to delete and avoid system permision issues
                shutil.rmtree(os.path.join(folder,del_name))
                os.mkdir(os.path.join(folder,dir_data))
                #If folder name does not exist it creates the folder name and data file
            else:
                os.mkdir(os.path.join(folder,dir_data))
        cmds_list = [[dir_exe,'-X', str(dx_bin_step), os.path.join(MW_settings_file,'MILDwave.xml')] for MW_settings_file in T_dir]
        procs_list = [Popen(cmd) for cmd in cmds_list]
        #procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE) for cmd in cmds_list] Use this code line changing PIPE for the simulation folder to write a log file including any errors.
        #Waits untill all the processes are finished to start a new process.
        for proc in procs_list:
            proc.wait()
            
            
def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])