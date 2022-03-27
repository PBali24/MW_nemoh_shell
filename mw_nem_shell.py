"""
mw_nem_shell.py
Created on Wed Oct 24 08:22:50 2018
@author: gveraofe
SCRIPT TO RUN A COUPLED SIMULATION
1. RUN NEMOH OR SELECT NEMOH RESULTS
2. RUN EB MW OR SELECT EXISTING EMPTY BASIN
3. RUN COUPLING FOLDERS OR SELECT A COUPLED FOLDER
4. INSTRUCTUIONS MUST BE PLACED HERE AND NOT IN THE FUNCTIONS
"""
#FUNCTIONS FOR RUNNING MW
import os
import numpy as np
import sys
# the mock-0.3.1 dir contains testcase.py, testutils.py & mock.py
# por si acaso assure that the direcotry which you are running is the nm_nemoh_shell directory with the files for running MILDwave
NEMOHDir = "E://NEM//nemoh_py_complete//NEMOH"  #DIRECTORY WHERE THE NEMOH PROGRAM IS LOCATED
CLI_dir = "C://MW_RUN//CLI_DIRECTION//MW_CLI.exe" #Location of the CLI exe for MILDwave
#FUNCIONS FOR RUNNING NEMOH IMPORT FROM MW NEMOH SHARED
sys.path.append(os.path.abspath(os.path.join(os.path.dirname("..//MW_NEM_shared"),".." )) )
from MW_NEM_shared import nem_shell_utilities as nsu
from MW_NEM_shared import nem_utilities as ne
from MW_NEM_shared import meshTypes as mt
#----------------------------------------------------------------------------------------
import mw_utilities_pre as pre
import mw_utilities_run as run
import mw_utilities_post as post
import coupling as cpl
import matplotlib.pyplot as plt
from matplotlib import cm
from xml.etree import ElementTree as et #Element Tree for modifying xml values
import time
import logging #writes log file withe the simulation times
#-----------------------------------------------------------------------------
#                                1.0 DEFINE SHELL PARAMETERS
#-----------------------------------------------------------------------------
print(f'Simulation Started: '+ str(time.asctime()))
MW_sim = {}#data structure containing the Simulation Parameters for MW-NEMOH Coupled MODEL
MW_sim['run_name']='b5_reg'  #this is the name common to all the MW dir runs and the NEMOH runs and is the name under which the main dicsts will be saved
#MW_sim['run_name']='cli_test'  #this is the name common to all the MW dir runs and the NEMOH runs and is the name under which the main dicsts will be saved
#MW_sim['mw_dir'] = ["E://phd_runz//ch4//mw//b1/","E://phd_runz//ch4//mw//b2/","E://phd_runz//ch4//mw//b3/"] # a list of directories for separate top level MILDwave runs i.e. cases [e.g. number of bodies, WEC types etc]
#NB to get names for plotting etc from the dir names make sure that the folder is not followed by a slash or two slashed ie '../../some_name'
#MW_sim['mw_dir'] = ["E://phd_runz//ch4//mw//b4_irr","E://phd_runz//ch4//mw//b5_irr","E://phd_runz//ch4//mw//b6_irr"]
MW_sim['mw_dir'] = ["E://phd_runz//ch7//mw//b5_reg//"]
MW_sim['create_ini'] = True
MW_sim['run_NEMOH'] = False#RUN A NEW NEMOH SIMULATION OR USED AND EXISTING ONE
MW_sim['run_MW_EB'] = False#RUN AN INCIDENT WAVE IN MILDwave OR USE AND EXISTING ONE
MW_sim['run_MW_CPL'] = True#RUN A PERTURBED WAVE IN MILDwave OR USE AND EXISTING ONE
MW_sim['Coupling'] = True#CREATE COUPLING FILES OR USE EXISTING COUPLING FILES IN THE FOLDER
MW_sim['regular'] = True #REGULAR WAVE
MW_sim['irregular'] = False#IRREGULAR WAVE
MW_sim['directional'] = False#DIRECTIONAL IRREGULAR WAVE
MW_sim['interpolation'] = False#TODO REMOVE INTERPOLATION

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(filename=os.path.join(os.path.abspath(os.path.join(MW_sim['mw_dir'][0],"..")),MW_sim['run_name'] +"_log.txt"),level=logging.DEBUG)

#-----------------------------------------------------------------------------
#                                1.1 DEFINE MW INPUT WAVE CONDITIONS
#----------------------------------------------------------------------------- 
MW_ini = {}
MW_ini = [{} for _ in range(len(MW_sim['mw_dir']))] #do not ask me how this is possible
if MW_sim['create_ini']:
    for ii in range(len(MW_ini)):
        MW_ini[ii]['mw_dir'] = MW_sim['mw_dir'][ii]#put the directories in the ini so they are readable From the outside level
        MW_ini[ii]['depth'] = [10]#depth has to be a float
        if MW_sim['regular']:
            MW_ini[ii]['T'] = [10.]  #NB have to take values out of lists
            MW_ini[ii]['H'] = [2. ]
            MW_ini[ii]['deg'] = [0.0]
            MW_ini[ii]['Nf'] = 2 #TODO search for use
            #TODO to do various frequencies at the same time from NEMOH
        elif MW_sim['irregular']:
            MW_ini[ii]['Tp'] = [6., 8.]
            MW_ini[ii]['Hs'] = [2. , 2.]      
            MW_ini[ii]['spectra'] = 'PM' #Define the type of Spectrum you want to use 'JS' for JONSWAP or 'PM' for Pierson-Moskovitz
            #TODO GAEL maybe eaiser entering range of perdios rather than freq
            MW_ini[ii]['fini'] = [2*(1/Tp) for Tp in MW_ini[ii]['Tp']]
            MW_ini[ii]['fend'] = [1/2*(1/Tp) for Tp in MW_ini[ii]['Tp']]
            MW_ini[ii]['Nf'] = 20
            MW_ini[ii]['deg'] = [[0.0],[0.0]]
            if MW_sim['directional']:
                MW_ini[ii]['deg'] = []
                MW_ini[ii]['deg_main'] = [0.0]
                MW_ini[ii]['s1'] = [18.37]    
        # call to createInput function

        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF            
        pre.create_input(MW_sim,MW_ini[ii])
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        tree = et.parse(os.path.join(MW_sim['mw_dir'][ii] ,"MILDwave.xml"))# MILDwave.xml input file 

        MW_ini[ii]['Nx']=int(tree.find(".//Nx").text) #TODO GAEL check
        MW_ini[ii]['Ny']=int(tree.find(".//Ny").text)
        #  left and right sponge layers
        MW_ini[ii]['Lsponge']=int(tree.find(".//ixsL").text)
        MW_ini[ii]['Rsponge']=int(tree.find(".//ixsR").text)
    
        MW_ini[ii]['dx_bin_step'] = 2#Increase of dx in the binary output files ->only the x value sets both for now in the mildwave grid 
        MW_ini[ii]['dy_bin_step'] = 2#Increase of dy in the binary otuput files
        
        MW_ini[ii]['kd_x'] = int(len(np.arange( MW_ini[ii]['Lsponge'], MW_ini[ii]['Nx'] - MW_ini[ii]['Rsponge'],MW_ini[ii]['dx_bin_step'] )))
        MW_ini[ii]['kd_y'] = int(len(np.arange(MW_ini[ii]['Lsponge'], MW_ini[ii]['Ny'] - MW_ini[ii]['Rsponge'],MW_ini[ii]['dy_bin_step'] )))
        # for complete domain == MILDwave domain size  -2 
        
    # snaphot if rest of simulations fails save initial conditions 
    MW = {'MW_sim' : MW_sim, 'MW_ini' : MW_ini}
    np.savez(os.path.join(os.path.abspath(os.path.join(MW_sim['mw_dir'][0],"..")),MW_sim['run_name'] +"_snap"), **MW)
    del MW
    
else: #TODO change to list of dicts
    for ii in range(len(MW_ini)):
        MW = np.load(os.path.join(os.path.abspath(os.path.join(MW_sim['mw_dir'][0],"..")),MW_sim['run_name'] +"_snap.npz"))
        MW_ini[ii] = MW['MW_ini'][ii] #have to do separately for it to be picakble (vis in variuable window)
    del MW
#-----------------------------------------------------------------------------
#                                1.2 DEFINE COUPLING SIMULATION PARAMETERS 
#-----------------------------------------------------------------------------
#COUPLING  all parameterss in cells exceopt where indicated
#Information for coupling
MW_CPL = [{} for _ in range(len(MW_sim['mw_dir']))] #do not ask me how this is possible
for ii in range(len(MW_sim['mw_dir'])):
    MW_CPL[ii]['circ'] = False #generates a circular coupling zone
    MW_CPL[ii]['rect'] = True#generates a rectangular coupling zone
    MW_CPL[ii]['output_check'] = False# GAEL ACTIVATE FOR PLOTTING THE COUPLING FIGURES TO CHECK THAT THE COUPLING IS DONE CORRECTLY
    MW_CPL[ii]['frequency_check'] = 0#
    # TODO GAEL why not in original MILDwave.xml
   #GAEL TODO list of lists for variaous locations eg [[185, 245] let's make it list of ints not list of lists 
	# incl soponge layers or no?
	#relative to complete domain (sans sponge layers)
    MW_CPL[ii]['ph_LOC_x'] = [int((MW_ini[ii]['Nx'] - MW_ini[ii]['Lsponge'] - MW_ini[ii]['Rsponge'])/2)-100,int((MW_ini[ii]['Nx'] - MW_ini[ii]['Lsponge'] - MW_ini[ii]['Rsponge'])/2)+100]
       #[[185]] Location where you want to couple the nem wave field with respect to the MW coordinates [0,0{top left}]
    MW_CPL[ii]['ph_LOC_y'] = [int((MW_ini[ii]['Ny'] - MW_ini[ii]['Lsponge'] - MW_ini[ii]['Lsponge'])/2),int((MW_ini[ii]['Ny'] - MW_ini[ii]['Lsponge'] - MW_ini[ii]['Lsponge'])/2)]#
	#[[121]] #Y CO-ORDINATES OF THE DIFFERENT COUPLING LOCATIONS
    if MW_sim['irregular']:
        MW_CPL[ii]['FFT_BLOCK'] = 256 #256 normally gives a good result.Other used values 128, 256, 512 or 1024 . The higher Nf the higher FFT_BLOCK
        MW_CPL[ii]['Overlap']   = 0.6 #Between 0 < Overlap < 1  normally a higher % in the overlap gives a better agreement between the target WAVE and the MW WAVE
        
    if MW_CPL[ii]['circ']:
        MW_CPL[ii]['radius'] = 50.0#coupling radius in (m) (radius in meters) radiusM = radius_coupling/dx = X.5
        #Example: nemoh(Ny,Nx)(201,201) centercircle_x=99 & centercircle_y=99 is
        #the central point of the domain AS PYTHON INDEXES FROM 0
    

        #TODO multiple locoations in a list??
    if MW_CPL[ii]['rect']:
        MW_CPL[ii]['width'] =  50 # it has to be defined in m
        MW_CPL[ii]['length'] = 50 # it has to be defined in m

        #X CO-ORDINATES OF THE DIFFERENT COUPLING LOCATIONS [[136-3 ,111]] 2 COUPLING LOCATIONS FOR THE SAME INCIDENT WAVE AND NEMOH RUN [[136-3] ,[111]] TWO DIFFERENT INCIDENT MILDWAVE WAVES AND NEMOH RUNS ARE REQUIRED
#-----------------------------------------------------------------------------
      
#-----------------------------------------------------------------------------
#                                1.3 DEFINE MW SIMULATION EXECUTABLE PARAMETERS 
#-----------------------------------------------------------------------------
MW_exe ={} 
# change name
MW_exe['mw_cli'] = CLI_dir #local machine directory with CLI executable
MW_exe['mw_dir'] = MW_sim['mw_dir']#each regular wave or irregular wave calculated requires its on directory
MW_exe['time_ramp'] = 200#expectd time for the ramping function to create the regular wave(for the longest time)
MW_exe['time_Kd'] = 802 #605 time frame to calculate the Kd (try +/- 1)
MW_exe['n_cores'] = 10#number of cores to run the irregular wave simulation
MW_OUT = {}#TO SAVE ALL MW OUTPUTS must be list of dictonaries of size of mw_dir (cases)
MW_OUT = [{} for _ in range(len(MW_sim['mw_dir']))] #do not ask me how this is possibleyy
#-----------------------------------------------------------------------------
#                                1.4 DEFINE NEMOH SIMULATION PARAMETERS and RUN NEMOH 
#-----------------------------------------------------------------------------

#NEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEM  
#    NEMOH PARAMETERS AND RUN START 
#NEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEM  
if MW_sim['run_NEMOH']:
    #data structure containing the Simulation Parameters for NEMOH
    NEM_ini = [{} for _ in range(len(MW_sim['mw_dir']))] #do not ask me how this is possible
        # different mw dirs must have the ability to have different nemoh runs
    for ii in range(len(MW_exe['mw_dir'])):        
        print(f'Simulation Started NEMOH: '+str(time.asctime()))
        logging.info('Simulation Started NEMOH: '+str(time.asctime()))
        #lets' set these and the wave conditions depth etc. equal to the mw_dir so reg=reg
        NEM_ini[ii]['case_dir'] = MW_sim['mw_dir'] #the name equivalent of NEMOH for the mw dirs had to be equal to the size of mw_dir list
        NEM_ini[ii]['runNEMOH'] = MW_sim['run_NEMOH']
        NEM_ini[ii]['regular'] = MW_sim['regular']
        NEM_ini[ii]['irregular'] = MW_sim['irregular']
        NEM_ini[ii]['directional'] = False#MW_sim['directional']
        NEM_ini[ii]['interpolation'] = MW_sim['interpolation']
        #-----------------------------------------------------------------------------
        #                                1.4.1 WAVE CONDITIONS
        #-----------------------------------------------------------------------------
        NEM_ini[ii]['depth'] = MW_ini[ii]['depth']  #TODO Should be set different for different mw_dirs but leave same for each for now
        #copy info from MW dict to NEM dict so the simulations parameters for NEMOH exist separately
        if NEM_ini[ii]['regular']:
            # NB this is set from the first case of the mw_dirs beware the periods and waveheights and degrees must be the same for all mw_dir cases
            NEM_ini[ii]['T'] = MW_ini[ii]['T']# [[1.26]]#define each T as list
            NEM_ini[ii]['H'] = MW_ini[ii]['H']#[[0.104]]#define each H as list
            NEM_ini[ii]['deg'] = MW_ini[ii]['deg']#[[0.0]]
            NEM_ini[ii]['Nf'] = 1
        elif NEM_ini[ii]['irregular']:
            #NB since we are NOT running NEMOH in multifreq Tp has got to be T for the loops in NEMOH shell to run correctly
            #Tp is only for completion sake and for plotting
            NEM_ini[ii]['T'] = MW_ini[ii]['T']#[1.26]
            NEM_ini[ii]['H'] = MW_ini[ii]['H']#[1.26]
            NEM_ini[ii]['Tp'] = MW_ini[ii]['Tp']#[1.26]
            NEM_ini[ii]['Hs'] = MW_ini[ii]['Hs']# [0.104]      
            NEM_ini[ii]['spectra'] =  MW_ini[ii]['spectra']#'JS' #Define the type of Spectrum you want to use 'JS' for JONSWAP or 'PM' for Pierson-Moskovitz
            NEM_ini[ii]['fini'] =  MW_ini[ii]['fini']#0.7#0.705
            NEM_ini[ii]['fend'] = MW_ini[ii]['fend']#2.1#2.15
            NEM_ini[ii]['Nf'] = MW_ini[ii]['Nf']#20
            NEM_ini[ii]['deg'] = MW_ini[ii]['deg']# [[0.0]]
            if MW_sim['directional']:
                NEM_ini[ii]['deg'] = []
                NEM_ini[ii]['deg_main'] = MW_ini[ii]['deg_main']
                NEM_ini[ii]['s1'] = MW_ini[ii]['s1']
    
    if MW_sim['directional']:
        NEM_ini[ii]['deg'] = MW_ini[ii]['deg']#WE HAVE TO USE THE SAME DIRECTIONS FOR MILDwave and NEMOH for the SHORT CRESTED WAVES
        #check for lists-<floats TODO
    #-----------------------------------------------------------------------------
    #                              1.4.2 GRID
    #-----------------------------------------------------------------------------  
    NEM_GRID = {}  #NEMOH input in meteres!!!
    NEM_GRID['Lg'] = 200 #lenght and width of the basin  (m)
    NEM_GRID['Wg'] = 200 #width of the basin (m)
    NEM_GRID['dx'] = 5#introduce the target deltax for NEMOH simulation
    NEM_GRID['dy'] = 5#introduce the target deltay for NEMOH simulation
    NEM_GRID['Nx'] = int(np.round(NEM_GRID['Lg']/NEM_GRID['dx']) + 1) #redifine number of gridpoints as Nxg = np.round(Lg/dx) + 1 Maximum number of GRID CELLS IN NEMOH 200
    NEM_GRID['Ny'] = int(np.round(NEM_GRID['Wg']/NEM_GRID['dy']) + 1) #redifine number of gridpoints as Nyg = np.round(Wg/dy) + 1
    if NEM_ini[ii]['interpolation']:
        NEM_GRID['dxN'] = 1.
        NEM_GRID['dxN'] = 1.
    #TODO !!! for now only change the body settings over the diofferent mw_dirs (cases) not the domains becase the .sxml imputs for domains are the same
    NEM_OUT=[{} for _ in range(len(MW_exe['mw_dir']))] #magic empty list of dicstionary creation
    for ii in range(len(MW_exe['mw_dir'])):
         #-----------------------------------------------------------------------------
        #                              1.4.3 NEM BODY
        #-----------------------------------------------------------------------------  
        NEM_BODY = [{} for _ in range(len(MW_exe['mw_dir']))] #do not ask me how this is possible
        NEM_BODY[ii]['dof'] = [0,0,1,0,0,0]#number of degrees of freedom for each simulation
        NEM_BODY[ii]['ndof'] = np.sum(NEM_BODY[ii]['dof'])
#        
#        NEM_BODY[0]['nbody'] = 5 #number of bodies
#        NEM_BODY[0]['xBody']= [30,-30,30,-30,30]#x-coordinates of each body (0,0) center of the domain
#        NEM_BODY[0]['yBody'] =  [60,30,0,-30,-60]  #y-coordinates of each body (0,0) center of the domain
    #   nbody = 5 #number of bodies
        NEM_BODY[0]['nbody'] = 1 #number of bodies
        NEM_BODY[0]['xBody'] = [0.] #x-coordinates of each body (0,0) center of the domain
        NEM_BODY[0]['yBody'] = [0]#y-coordinates of each body (0,0) center of the domain
        
        # NEM_BODY[2]['nbody'] = 3 #number of bodies
        # NEM_BODY[2]['xBody']= [25,25. , -25.] #x-coordinates of each body (0,0) center of the domain
        # NEM_BODY[2]['yBody'] = [ 25 ,-25., 0.]#y-coordinates of each body (0,0) center of the domain
        
        cylmesh = [0]*NEM_BODY[ii]['nbody']
    #    ob1 = [0]*NEM_BODY['nbody']
    #    ob2 = [0]*NEM_BODY['nbody']
        NEM_BODY[ii]['cG'] = -10.   #z-coordinate of the gravity center of each body 
        nPanels = 200   #number of panels for the NEMOH mesh
        nsym = 0#MESH THE BODY WITH SIMETRY AXIS
        NEM_BODY[ii]['rho'] = 1025.0
        NEM_BODY[ii]['PTOtype'] = 'Bpto_L'#Theory_LINEAR = Bpto_L WECSIM_LINEAR = Bpto_wsL WECSIM_HYDRA = Bpto_wsH
        NEM_BODY[ii]['Bpto'] = 121 * 10**6
        NEM_BODY[ii]['ptoProp'] = [0.0,NEM_BODY[ii]['Bpto'],0.0]#[Mextra,Bpto,Kextra]
        #-----------------------------------------------------------------------------
        #                              1.4.4 MESHING
        #----------------------------------------------------------------------------+-

        os.chdir(os.path.join(NEMOHDir,'Calculation'))#In the directory Calculation all the NEMOH runs for each frequency will be calculated
                               #Regular wave calculations folder is saved with corresponding H and T
                               #The results folder for irregular wave is changed for a folder with name Tp and Hs
        if not(os.path.isdir('mesh')) :    
            os.mkdir('mesh')
        # and f'inresults forlder                
        if not(os.path.isdir('results')) :    
            os.mkdir('results')
        for iB in range(NEM_BODY[ii]['nbody']):
            #ob1[iB] = mt.cylinder(0.315,0.1575+0.0082,[NEM_BODY[ii]['xBody'][iB],NEM_BODY[ii]['yBody'][iB],0.0])#diameter and draft
            #ob2[iB] = mt.hemisphere (0.315,[NEM_BODY[ii]['xBody'][iB],NEM_BODY[ii]['yBody'][iB],-0.1575-0.0082])#diameter 
            #cylmesh[iB] = mt.Mesh()
            #cylmesh[iB].combineMesh(ob1[iB],ob2[iB])
            #TODO this also has to be customaizable for different cases
            #HPA
#            cylmesh[iB] = mt.cylinder(10.0,2.0,[NEM_BODY[ii]['xBody'][iB],NEM_BODY[ii]['yBody'][iB],0.0])#diameter and draft
            #OSWEC
            cylmesh[iB] = mt.box(1.0,20.0, 12.0,[NEM_BODY[ii]['xBody'][iB],NEM_BODY[ii]['yBody'][iB],-4])# thickness width height

            if NEM_BODY[ii]['nbody'] > 1:
                mt.writeMesh(cylmesh[iB],'./mesh/axisym{0:d}'.format(iB+1))
            else:
                mt.writeMesh(cylmesh[iB],'./mesh/axisym') 
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        ne.createMeshOpt(NEM_BODY[ii]['cG'],nPanels,nsym,rho=NEM_BODY[ii]['rho'],g=9.81,nbody=NEM_BODY[ii]['nbody'],
                         xG=NEM_BODY[ii]['xBody'],yG=NEM_BODY[ii]['yBody'])#CALLS THE MESH FUNCION INCLUDED IN NEMOH
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        (NEM_BODY[ii]['Mass'],NEM_BODY[ii]['Kh']) = ne.calcM(rho=1026.0)
    	
    # =============================================================================
    #     ne.plotmesh(NEM_BODY[ii]['nbody'],NEMOHDir,MW_CPL[ii]['radius'],NEM_BODY[ii]['xBody'],NEM_BODY[ii]['yBody'])
    #     
    #     print('press 1 if the NEMOH mesh is correct')
    #     mesh=int(input())
    #     if mesh==1:
    #         print('Continue')
    #         plt.close('all')
    #     else:
    #         raise ValueError('Mesh is not correct')
    # #	del ob1,ob2,iB,nPanels,nsym,cylmesh  for 3 bodies	
    #     del iB,nPanels,nsym,cylmesh #for 1 body
    #     #TODO PLOT THE MESH GEOMETRY BEFORE RUNNING NEMOH
    # =============================================================================
        #-----------------------------------------------------------------------------
        #                              1.4.5. NEMOH SIMULATION PARAMETERS (ADVANCED OPTIONS)
        #-----------------------------------------------------------------------------
        # Basic Options (RAO calculation)
        nrFreq = 1         #each NEMOH simulation is run with one frequency TODO run NEMOH for several frequencies
        # Advanced  TODO keep same for differenct cases for now
        NEM_advOps = {}
        NEM_advOps['rhoW'] = 1025.0      #water density
        NEM_advOps['dirCheck'] =  True#Activate to change wave direction
        NEM_advOps['dirStep'] = 1
        NEM_advOps['irfCheck']  = False#Activate for IRF calculations
        NEM_advOps['irfDur'] = 40.0
        NEM_advOps['irfStep'] = 0.01
        NEM_advOps['kochCheck'] = False#Activate to calculate Kochin Function
        NEM_advOps['kochStart'] = 0.0
        NEM_advOps['kochStop'] = 360.0
        NEM_advOps['kochStep'] = 24
        NEM_advOps['fsCheck'] = True#Activate to calculate free surface elevation
        NEM_advOps['fsNx'] = NEM_GRID['Nx']
        NEM_advOps['fsNy'] = NEM_GRID['Ny']
        NEM_advOps['fsLengthX'] = NEM_GRID['Lg']
        NEM_advOps['fsLengthY'] = NEM_GRID['Wg']  
        NEM_advOps['Show_Console'] = False #Toggle ON or OFF the WINDOWS CONSOLE
        #-----------------------------------------------------------------------------
        #                              1.4.6. RUNNING NEMOH
        #-----------------------------------------------------------------------------
#        # call runNEM for each of The nem dirs (ii)
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        (NEM_OUT[ii]) = nsu.runNEM(nrFreq,NEM_ini[ii],NEM_BODY[ii],NEM_advOps,NEM_GRID,NEMOHDir)
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

#-----------------------------------------------------------------------------
#                              1.4.7. GENERATING COUPLING OUTPUTS FROM NEMOH
#-----------------------------------------------------------------------------

    NEM = {'NEM_ini' : NEM_ini, 'NEM_GRID' :NEM_GRID, 'NEM_BODY' : NEM_BODY, 'NEM_advOps' : NEM_advOps, 'NEM_OUT' : NEM_OUT}
    np.savez(os.path.join(os.path.abspath(os.path.join(MW_sim['mw_dir'][0],'..')),'nem_'+MW_sim['run_name']),**NEM)
    print(f'Simulation Finished NEMOH: '+str(time.asctime()))
    logging.info('Simulation Finished NEMOH: '+str(time.asctime()))
    
# NB makje sure that the NEMOGH oaded has  been run over the same number of mw_dirs as the present run
else:#we load and existing NEMOH run from the root dir with the NEM dictionary which is of dimension of mw_dir (ii) (cases)
            #TODO figure out a better root dir name which will apply to both NEM and MW Dicts npz 
        nemoh_dir = os.path.abspath(os.path.join(MW_sim['mw_dir'][0], '..')) #we get the parent dir of the 1st mw-dir which must always exist
        NEM = np.load(os.path.join(nemoh_dir,'nem_'+ MW_sim['run_name']+'.npz')) # TODO check to make sure this actually works
        NEM_ini = NEM['NEM_ini']
        NEM_GRID = NEM['NEM_GRID'].item()
        NEM_BODY = NEM['NEM_BODY']
        NEM_advOps = NEM['NEM_advOps'].item()
        NEM_OUT = NEM['NEM_OUT']
        
#NEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEM          
#  NEMOH PARAMETERS AND RUN END        
#NEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEMNEM    
#-----------------------------------------------------------------------------
#                                1.5.1 RUN MW EMPTY BASIN
#-----------------------------------------------------------------------------   
if MW_sim['run_MW_EB']:
    # regular wave run

        # initialization functions that ovewrite MILDwave.xml files
    for ii in range(len(MW_sim['mw_dir'])):  
        if MW_sim['regular']:
            #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
            pre.reg_init(MW_ini[ii],MW_exe)#creates the dummy MILDwave.xml file for a irregular wave simulation or de desired MILDwave.xml file for a regular simulation
            # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
            mwd = MW_exe['mw_dir'][ii]
            for Tjj in MW_ini[ii]['T']:
                run_type ='EB'
                print(f'Simulation Incident Wave MW Started: '+ str(time.asctime()) +' for case '+mwd +'for T = '+str(Tjj))
                logging.info('Simulation Incident Wave MW Started: '+str(time.asctime())+' for case '+mwd +'for T = '+str(Tjj))
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                run.run_mw(MW_sim,MW_ini[ii],mwd,Tjj,MW_exe,run_type,Tnn="")
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                print(f'Simulation Incident Wave MW Finished: '+ str(time.asctime()) +' for case '+mwd +'for T = '+str(Tjj))
                logging.info('Simulation Incident Wave MW Finished: '+str(time.asctime())+' for case '+mwd +'for T = '+str(Tjj))
                
        elif MW_sim['irregular']: #TODO not checked!!!!
            pre.irr_init(MW_ini[ii],MW_exe)#creates all the MW_exe files for running a irregular empty simulation in MILDwave from the Tp MILDwave.xml file
        #    only running over various mw_dirs and various periods T (or Tp) for now    
            mwd = MW_exe['mw_dir'][ii]
            for jj in range(len(MW_ini[ii]['Tp'])):
                Tjj = MW_ini[ii]['Tp'][jj]
                Tnn = MW_ini[ii]['T'][jj]
                run_type ='EB'
                print(f'Simulation Incident Wave MW Started: '+ str(time.asctime()) +' for case '+mwd +'for T = '+str(Tjj))
                logging.info('Simulation Incident Wave MW Started: '+str(time.asctime())+' for case '+mwd +'for T = '+str(Tjj))
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                run.run_mw(MW_sim,MW_ini[ii],mwd,Tjj,MW_exe,run_type,Tnn)
                # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
                print(f'Simulation Incident Wave MW Finished: '+ str(time.asctime()) +' for case '+mwd +'for T = '+str(Tjj))
                logging.info('Simulation Incidennt Wave MW Finished: '+str(time.asctime())+' for case '+mwd +'for T = '+str(Tjj))

#----------------------------------------------------------------------[ii]-------
#                                1.5.2 CALCULATE MW INCIDENT WAVE KD
#-----------------------------------------------------------------------------   
# POSTPRO - run the dir and T loops inside function
run_type = 'EB'
# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
for ii in range(len(MW_exe['mw_dir'])): # run over the directories outside but ove the periods (and ccoupling locations indside the fcuntion)
#because the size of the sdimualtion domain may change over the dirs but not over the T or Tp's
    (MW_OUT[ii]['Kd_EB'],MW_OUT[ii]['xvect'],MW_OUT[ii]['yvect']) = post.kd_EB(MW_sim,MW_ini[ii],MW_exe)

# FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF   
    
 #=============================================================================
#-----------------------------------------------------------------------------
#                                1.5.3 RUN MW PERTURBED WAVE SIMULATION
#-----------------------------------------------------------------------------   
#=============================================================================
if MW_sim['run_MW_CPL']:
#   for jj,mw_dir,T in zip(range(len(MW_exe['mw_dir'])),MW_exe['mw_dir'],MW_ini['T']):
    for ii in range(len(MW_sim['mw_dir'])):
        if MW_sim['Coupling']:
            # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
            cpl.coupling(MW_sim,MW_ini[ii],MW_exe,MW_OUT[ii],MW_CPL[ii],NEM_ini,NEM_GRID,NEM_OUT[ii])
            # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        # we use the run function with The outside jj (period) loop
        mwd = MW_exe['mw_dir'][ii]
        if MW_sim['regular'] :
             for jj in range(len(MW_ini[ii]['T'])):
                 Tjj=MW_ini[ii]['T'][jj]
                 print(f'Simulation Coupled Wave MW Started: '+ str(time.asctime()) +' for case '+mwd +'for T = '+str(Tjj))
                 logging.info('Simulation Coupled Wave MW Started: '+str(time.asctime())+' for case '+mwd +'for T = '+str(Tjj))
                 for kk in range(len(MW_CPL[ii]['ph_LOC_x'])): # kk is various coupling lcoations
                     run_type = 'CP_'+'LOC_'+str(kk)
                     run.run_mw(MW_sim,MW_ini[ii],mwd,Tjj,MW_exe,run_type,Tnn="") 
                     #TODO add for location
                 print(f'Simulation Coupled Wave MW Finished: '+str(time.asctime())+' for case '+ mwd +'for T =',str(Tjj))
                 logging.info('Simulation Coupled Wave MW Finished: '+str(time.asctime())+' for case ' +mwd +'for T =',str(Tjj))
        elif MW_sim['irregular'] :
            for jj in range(len(MW_ini[ii]['Tp'])):
                Tjj = MW_ini[ii]['Tp'][jj]
                Tnn = MW_ini[ii]['T'][jj]
                for kk in range(len(MW_CPL[ii]['ph_LOC_x'])): # kk is various coupling lcoations
                    run_type = 'CP_'+'LOC_'+str(kk)
                    run.run_mw(MW_sim,MW_ini[ii],mwd,Tjj,MW_exe,run_type,Tnn)
                print(f'Simulation Coupled Wave MW Finished: '+str(time.asctime())+' for case '+ mwd +'for T =',str(Tjj))
                logging.info('Simulation Coupled Wave MW Finished: '+str(time.asctime())+' for case ' +mwd +'for T =',str(Tjj))
#=============================================================================
#-----------------------------------------------------------------------------
#                                1.5.4 CALCULATE KD PERTURBED WAVE SIMULATION
#-----------------------------------------------------------------------------   
#=============================================================================            

for ii in range(len(MW_sim['mw_dir'])):
    #the periods and the cpoupling lcoations are handled inside
        (MW_OUT[ii]['Kd_CP'],MW_OUT[ii]['xvect'],MW_OUT[ii]['yvect']) = post.kd_CP(MW_sim,MW_ini[ii],MW_CPL[ii],MW_exe)
   
#-----------------------------------------------------------------------------
#                                1.5.5 CALCULATE MW TOTAL WAVE
#-----------------------------------------------------------------------------  
# POSTPRO run loops over mwdirs and Ts inside funxion
t1_tot=time.time()
 # run over the directories outside but ove the periods (and ccoupling locations indside the fcuntion)
#because the size of the sdimualtion domain may change over the dirs but not over the T or Tp's
for ii in range(len(MW_exe['mw_dir'])):
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        MW_OUT[ii]['Kd_TW'] = post.total_wave(MW_sim,MW_ini[ii],MW_exe,MW_CPL[ii])
        # FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
t2_tot=time.time()
#create dictionary of dictionaries for saving all data
MW = {'MW_sim' : MW_sim, 'MW_ini' : MW_ini, 'MW_exe' : MW_exe, 'MW_CPL' : MW_CPL,'MW_OUT' : MW_OUT}
#for now save to change to superdirectory i.e. root_dir of the list of mw_dirs 
np.savez(os.path.join(os.path.abspath(os.path.join(MW_sim['mw_dir'][0],"..")),MW_sim['run_name']+'_mw_full_run') ,**MW)
logging.info('Simulation Finished: '+str(time.asctime()))

#-----------------------------------------------------------------------------
#                                1.6. PLOTTING RESULTS
#-----------------------------------------------------------------------------  
#mw_plot.MW_CPL[ii]ontourf(MW_OUT['Kd_EB'][0][0],MW_OUT['x_vector'],MW_OUT['y_vector'],0,elev_min = 0.8,elev_max = 1.2,mid_val = 1.0,title = 'Kd')
#mw_plot.MW_CPL[ii]ontourf(MW_OUT['Kd_CP'][0][0],MW_OUT['x_vector'],MW_OUT['y_vector'],0,elev_min = 0,elev_max = 0.35,mid_val = 0.17,title = 'Kd')
#mw_plot.MW_CPL[ii]ontourf(MW_OUT['Kd_TW'][0][0],MW_OUT['x_vector'],MW_OUT['y_vector'],MW_CPL[ii]['radius'],
#   elev_min = 0.8,elev_max = 1.2,mid_val = 1.0,title = 'Kd')
#plot for each mw_dior and each perdio

cmap = cm.plasma       
for ii in range(len(MW_exe['mw_dir'])):
    if MW_sim['regular']:          
            Tpl = MW_ini[ii]['T']
    else:
            Tpl = MW_ini[ii]['Tp']  
    for jj in range(len(Tpl)):
        plt.figure()    
        plt.title(MW_exe['mw_dir'][ii]+str(Tpl[jj])+'EB')
        plt.contourf(MW_OUT[ii]['xvect'],MW_OUT[ii]['yvect'],MW_OUT[ii]['Kd_EB'][jj][0])    
        plt.colorbar()
        ##     total wave        
        plt.figure()    
        plt.title(MW_exe['mw_dir'][ii]+str(Tpl[jj])+'_TW')
        plt.contourf(MW_OUT[ii]['xvect'],MW_OUT[ii]['yvect'],MW_OUT[ii]['Kd_TW'][jj][0],cmap=cm.get_cmap(cmap))    
        plt.colorbar()
        for kk in range(len(MW_CPL[ii]['ph_LOC_x'])): #TODO for when various coupling locations are simulated 
            # CP for each coupling lcoation
            plt.figure()                
            plt.title(MW_exe['mw_dir'][ii]+str(Tpl[jj])+'_CP_location'+str(MW_CPL[ii]['ph_LOC_x']))
            plt.contourf(MW_OUT[ii]['xvect'],MW_OUT[ii]['yvect'],MW_OUT[ii]['Kd_CP'][jj][kk][0],cmap=cm.get_cmap(cmap))    
            plt.colorbar()
        
        
'''    
plt.figure()
plt.plot(MW_OUT['x_vector'][137,:]+0.32,MW_OUT['Kd_TW'][0][0][0][34,:])
plt.plot(NEM_OUT['xvect'][49,:],NEM_OUT['Kd_time'][0][99,:])

plt.ylim(0.8,1.2)
'''