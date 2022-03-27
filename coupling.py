"""
coupling.py the functions for coupling NEMOH to MILDwave 
@author: gveraofe
mod on Thu Feb 28 18:20:25 2019 PhilDog
incldues 
coupling
coupling_input_irregular
"""
import numpy as np
import os
from xml.etree import ElementTree as et #Element Tree for modifying xml values
import shutil
from scipy import signal
import fnmatch #GALE was ist dies?
def coupling(MW_sim,MW_ini,MW_exe,MW_OUT,MW_CPL,NEM_ini,NEM_GRID,NEM_OUT):
    """
    GVF 24/01/2017
    This function calculates the coupling input files for MILDwave
    nemoh_data.txt contains the pert_wave generation circle
    It does not consider interaction between the DOF of the buoy(RAO
    calculation equation)
    NEMOH domain needs to have an odd number of cells in both directions
    MILDwave domain needs to have an odd number of cells in both directions
    The pert_wave generation circle is obtained around a single grid cell
    The pert_wave generation rectangle is obtained from the top left corner cell
    GRID SIZE for NEMOH and MILDwave has to be the same
    NEMOH results can be interpolated to a finner grid using post_nemoh.py
    COMMON INPUT
    phase data at the coupling locations
    wave amplitude at the coupling locations
    phase_dimensions = Nf * Locations
    Nf = number of frequencies
    Locations = locations
    for regular waves dimensions is amp and phase 1xLocations
    for irregular wave dimensions for amp and phase is 1xLocations
    TODO INCLUDE R_AMP CREATION FOR DIFFERENT LOCATIONS EACH LOCATION WILL HAVE A DIFFERENT ESPECTRUM BASED ON Hs
    Inputs circular coupling
    radius (m) Total radius of the pert_wave generation circle (COUPLING RADIUS)
    centercircle_x = x-coordinate of the center of the pert_wave generation circle
    in the NEMOH domain in grid cells
    centercircle_y = y-coordinate of the center of the pert_wave generation circle
    in the NEMOH domain in grid cells
    tlc_x_mw = x-coordinates of the center of the coupling location in the MILDwave domain in m for circular coupling
    tlc_y_mw = y-coordinates of the center of the coupling location in the MILDwave domain in m for circular coupling
    Inputs rectangular coupling
    tlc_x_mw = x-coordinates of the top left corner of the coupling location in the MIDLwave doamin in cells for rectangular coupling
    tlc_y_mw = y-coordinates of the top left corner of the coupling location in the MIDLwave doamin in cells for rectangular coupling
    width = width in number of cells of the rectangular coupling region
    length = length in number of cells of the rectangular coupling region
    Nx_mw = x-grid size of mildwave domain
    Ny_mw = y-grid size of the mildwave domain
    """
#-----------------------------------------------------------------------------
#                                  2. CREATING PERTURBED WAVE FOR COUPLNG
#-----------------------------------------------------------------------------
#    for mw_dir,kd_eb,a,T,f,loc_x,loc_y,pert_wave,xoff,yoff in zip(MW_exe['mw_dir'],MW_OUT['Kd_EB'],MW_ini['amp'],MW_ini['T'],MW_ini['f'],MW_CPL['LOC_x'],MW_CPL['LOC_y'],NEM_OUT['pert_wave'],MW_CPL['xoff'],MW_CPL['yoff']):
#    TODO replace zip
    mwd=MW_ini['mw_dir'] # MW ini is an instance for each directory (ii) outside
                       
    for jj in range(len(MW_ini['T'])):
        Tjj = MW_ini['T'][jj] 
        # intilaize variables instead of zip
        kd_eb,wave_amp,loc_x,loc_y,pert_wave =  MW_OUT['Kd_EB'][jj],MW_ini['amp'][0], MW_CPL['ph_LOC_x'],MW_CPL['ph_LOC_y'],NEM_OUT['pert_wave'][jj]
        #name ofsubfolder for each period (already created in reg_init or irreg_init so need to check pour l'existence
        #offset in number of cells from the center of the MILDwave domain to the coupling LOCATION
        xoff = [locx - int((MW_ini['Nx']-MW_ini['Lsponge'] - MW_ini['Rsponge'])/2) for locx in MW_CPL['ph_LOC_x']]#-int(1) to account for Python 0 indexing
        #[[0]] #offset in number of cells from the center of the MILDwave domain to the coupling LOCATION 
        yoff = [locy - int((MW_ini['Ny']-MW_ini['Lsponge'] - MW_ini['Rsponge'])/2) for locy in MW_CPL['ph_LOC_y']]#-int(1) to account for Python 0 indexing
        # GAEL ???? wave_ampl = ?
        if MW_sim['regular']:
            wave_ampw = [wave_amp]*len(loc_x)
            f= 1/Tjj #TODFO check that this corresponds to the MWW_ini['f'][jj]
            dir_T = os.path.join(mwd,'T_'+'{:05.2f}'.format(Tjj))
        else:#calculate FFT`1`
            f = [1/Tj for Tj in Tjj] #TODFO check that this corresponds to the MWW_ini['f'][jj]
            dir_T = os.path.join(mwd,'Tp_'+'{:05.2f}'.format(MW_ini['Tp'][jj]))
            wave_ampw = coupling_input_irregular(MW_ini['Nf'],f,dir_T,MW_exe,MW_ini,MW_CPL)
            
        for kk in range(len(loc_x)):#TODO SAVE MORE COUPLING REGIONS IN THE SAME NEMOH txt file
            dir_cp = 'CP_'+'LOC_'+str(kk)
            #and still copy existing MILDwave.xml file to subfodler for each period to account for hald run simulations
            if os.path.exists(os.path.join(dir_T,dir_cp,'data')):
                del_name = 'id1'#Auxiliar name for deleting existing folder
                os.rename(os.path.join(dir_T,dir_cp,'data'),os.path.join(dir_T,dir_cp,del_name))#renames existing folder to delete and avoid system permision issues
                shutil.rmtree(os.path.join(dir_T,dir_cp,del_name),ignore_errors=True)
                del_name = 'id1'#Auxiliar name for deleting existing folder
                os.rename(os.path.join(dir_T,dir_cp),os.path.join(dir_T,del_name))#renames existing folder to delete and avoid system permision issues
                shutil.rmtree(os.path.join(mwd,del_name),ignore_errors=True)
                os.mkdir(os.path.join(dir_T,dir_cp))
                os.mkdir(os.path.join(dir_T,dir_cp,'data'))
            elif os.path.exists(os.path.join(dir_T,dir_cp)):
                del_name = 'id1'#Auxiliar name for deleting existing folder
                os.rename(os.path.join(dir_T,dir_cp),os.path.join(dir_T,del_name))#renames existing folder to delete and avoid system permision issues
                shutil.rmtree(os.path.join(dir_T,del_name,dir_cp),ignore_errors=True)
                os.mkdir(os.path.join(dir_T,dir_cp))
                os.mkdir(os.path.join(dir_T,dir_cp,'data'))
            else:
                if not(os.path.isdir(os.path.join(dir_T,dir_cp))):
                    os.mkdir(os.path.join(dir_T,dir_cp)) # GAEL this fails
                os.mkdir(os.path.join(dir_T,dir_cp,'data'))
            #copy the three horesemen of the apocalyse to the created coupling directory for each location
            #shutil.copyfile(os.path.join(dir_T,"MILDwave.xml"),os.path.join(dir_T,dir_cp,"MILDwave.xml")),
            #shutil.copyfile(os.path.join(dir_T,"depth.dat"),os.path.join(dir_T,dir_cp,"depth.dat"))
            #shutil.copyfile(os.path.join(dir_T,"obstacle.dat"),os.path.join(dir_T,dir_cp,"obstacle.dat"))
                
            if MW_sim['regular']:
                
                dir_EB = os.path.join(dir_T,'EB')
                dir_cp_full = os.path.join(dir_T,dir_cp)
                phi_name = fnmatch.filter(os.listdir(os.path.join(dir_EB,'data')), 'phi*')
                phi_name.sort()
                eta_name = fnmatch.filter(os.listdir(os.path.join(dir_EB,'data')), 'ephi*')
                eta_name.sort()
                with open(os.path.join(dir_EB,'data',phi_name[0])) as fphi_aux:
                    phi_aux = np.array([line.strip().split() for line in fphi_aux],float)
                with open(os.path.join(dir_EB,'data',eta_name[0])) as feta_phi:
                    eta_phi = np.array([line.strip().split() for line in feta_phi],float)
                eta_instant = np.zeros((np.shape(eta_phi)[0],np.shape(eta_phi)[1]),dtype=np.complex_)
                tree = et.parse(os.path.join(dir_EB,"MILDwave.xml"))
                Tw = np.float(tree.find(".//wavePeriod").text)
                w = 2*np.pi / Tw
                eta_instant = eta_phi+1j*phi_aux*w/9.81
                phase = np.angle(eta_instant)
                pert_wave_aux=pert_wave*np.exp(1j*phase[[loc_y[kk]],[loc_x[kk]]]) 
                pert_wave_coup = wave_amp* kd_eb[0][ loc_y[kk]//MW_ini['dy_bin_step'],loc_x[kk]//MW_ini['dx_bin_step']]* pert_wave_aux

                tree.find(".//waveHeight").text = '%.6f'% (wave_amp* kd_eb[0][ loc_y[kk]//MW_ini['dy_bin_step'],loc_x[kk]//MW_ini['dx_bin_step']]*2.0)
                tree.find(".//generationType").text = '%.0f'%1 #generationType 1 in MW is coupling
                #tree.find(".//ixsL").text = '%0f'%ixsL_new
                #tree.find(".//ixsR").text = '%0f'%ixsR_new
                jysT =  int(tree.find(".//ixsL").text)
                tree.find(".//jysT").text = '%.0f'%jysT
                tree.find(".//jysB").text = '%.0f'%jysT
                tree.write(os.path.join(dir_cp_full,"MILDwave.xml"))    
                shutil.copyfile(os.path.join(dir_T,"depth.dat"),os.path.join(dir_cp_full,"depth.dat"))
                shutil.copyfile(os.path.join(dir_T,"obstacle.dat"),os.path.join(dir_cp_full,"obstacle.dat"))
                #if WGCheck:
                #copyfile(os.path.join(dirname_eb,dir_name_Tp,"Pos_WG.xml"),os.path.join(dir_MW_CPLp,dir_name,"Pos_WG.xml"))
                Lgmw = (int(tree.find(".//Nx").text)-1)*float(tree.find(".//dx").text)
                Wgmw = (int(tree.find(".//Ny").text)-1)*float(tree.find(".//dy").text)

                xFS = np.arange(-NEM_GRID['Lg']/2.0, NEM_GRID['Lg']/2.0+NEM_GRID['dx'],NEM_GRID['dy'])
                xFS = xFS[:NEM_GRID['Nx']]
                yFS = (-1)*np.arange(-NEM_GRID['Wg']/2.0, NEM_GRID['Wg']/2.0+NEM_GRID['dy'],NEM_GRID['dy'])
                yFS = yFS[:NEM_GRID['Ny']]
                #I multiply yFS times (-1) to get the correct values for the coordiante matrix in m
                X,Y = np.meshgrid(xFS,yFS) #TODO figure out how to do for a non-circular grid
                xFS_mw = np.arange(-Lgmw/2.0, Lgmw/2.0+NEM_GRID['dx'],NEM_GRID['dx'])
                #xFS_mw = xFS_mw[:Nx_mw]
                yFS_mw = (-1)*np.arange(-Wgmw/2.0, Wgmw/2.0+NEM_GRID['dy'],NEM_GRID['dy'])#there is no need to multiply here in order to obtain teh right input position for MILDwave
                #yFS_mw = yFS_mw[:Ny_mw]                
                Xmw,Ymw = np.meshgrid(xFS_mw,yFS_mw)
                
                if MW_CPL['circ']:
                    center_x = 0
                    center_y = 0
                    center_xmw = 0
                    center_ymw = 0
                    
                    mask1 = np.sqrt((X-center_x)**2+(Y-center_y)**2) < MW_CPL['radius'] + NEM_GRID['dx'] #it is necessary to add the extra  NEM_GRID['dx'] to obtain the maxium radius
                    mask2 = np.sqrt((Xmw-center_xmw)**2+(Ymw-center_ymw)**2) < MW_CPL['radius'] + NEM_GRID['dy']
    
                    yynm, xxnm = np.where(mask1)
                    yymw, xxmw = np.where(mask2)
    
                    
                elif MW_CPL['rect']:  
                    xxh = int(NEM_GRID['Nx']/2) - int(MW_CPL['length']/(NEM_GRID['dx'] *2))
                    yyh = int(NEM_GRID['Ny']/2) - int(MW_CPL['width']/(NEM_GRID['dy']*2))
                    xxh_mw =  int(MW_ini['Nx']/2) - int(MW_CPL['length']/(NEM_GRID['dx']*2))
                    yyh_mw =  int(MW_ini['Ny']/2) - int(MW_CPL['width']/(NEM_GRID['dy']*2))
                    mask1 = np.zeros((np.shape(Y)[0],np.shape(X)[1]), dtype=bool)
                    mask1[yyh:yyh+int(MW_CPL['width']/(NEM_GRID['dy'])) ,xxh:xxh+ int(MW_CPL['length']/(NEM_GRID['dx']))] = True
                    mask2 = np.zeros((np.shape(Ymw)[0],np.shape(Xmw)[1]), dtype=bool)
                    mask2[yyh_mw:yyh_mw+int(MW_CPL['width']/(NEM_GRID['dy'])) ,xxh_mw:xxh_mw+int(MW_CPL['length']/(NEM_GRID['dx']))] = True
                    yynm, xxnm = np.where(mask1)
                    yymw, xxmw = np.where(mask2)
                    
                amp = np.zeros((np.shape(xxnm)[0]))
                phase = np.zeros((np.shape(xxnm)[0]))
                amp[:] = np.abs(pert_wave_coup[yynm,xxnm])
                phase[:] = np.angle(pert_wave_coup[yynm,xxnm])
    
            #-----------------------------------------------------------------------------
            #                                  5. CREATING NEMOH_DATA.TXT
            #-----------------------------------------------------------------------------
                fid = open(os.path.join(dir_cp_full,'nemoh_data.txt'),'w')  #one frequency
                for mm in range(np.shape(amp)[0]):
                    fid.write('{:.0f}  {:.0f}  {:.8f}  {:.8f}\n'.format(xxmw[mm]+xoff[kk],yymw[mm]+yoff[kk],amp[mm],phase[mm]))
                fid.close()
                
            else: # irregular 
                #check where the xml depth and obst files are copied from               
                dir_EB = [os.path.join(mwd,dir_T,'EB','T_'+'{:06.3f}'.format(Tj)) for Tj in Tjj ]
                dir_cp_full = [os.path.join(mwd,dir_T,dir_cp,'T_'+'{:06.3f}'.format(Tj)) for Tj in Tjj ]
            #general for loop for coupling
            # dir Tp is taken care of by the T loop in ll 59-60 dirCP  loop kk in l79                 
                for nn in range(len(dir_cp_full)): 
                    
                    if os.path.exists(dir_cp_full[nn]):
                        del_name = 'id1'#Auxiliar name for deleting existing folder
                        os.rename(os.path.join(dir_cp_full[nn]),os.path.join(dir_T,dir_cp,del_name))#renames existing folder to delete and avoid system permision issues
                        shutil.rmtree(os.path.join(dir_T,dir_cp,del_name),ignore_errors=True)
                        os.mkdir(dir_cp_full[nn])
                        #If folder name does not exist it creates the folder name and data file
                    else:
                        os.mkdir(dir_cp_full[nn])
                        
                    phi_name = fnmatch.filter(os.listdir(os.path.join(dir_EB[nn],'data')), 'phi*')
                    phi_name.sort()
                    eta_name = fnmatch.filter(os.listdir(os.path.join(dir_EB[nn],'data')), 'ephi*')
                    eta_name.sort()
                    with open(os.path.join(dir_EB[nn],'data',phi_name[0])) as fphi_aux:
                        phi_aux = np.array([line.strip().split() for line in fphi_aux],float)
                    with open(os.path.join(dir_EB[nn],'data',eta_name[0])) as feta_phi:
                        eta_phi = np.array([line.strip().split() for line in feta_phi],float)
                    eta_instant = np.zeros((np.shape(eta_phi)[0],np.shape(eta_phi)[1]),dtype=np.complex_)
                    tree = et.parse(os.path.join(dir_EB[nn],"MILDwave.xml"))
                    Tw = np.float(tree.find(".//wavePeriod").text)
                    w = 2*np.pi / Tw
                    eta_instant = eta_phi+1j*phi_aux*w/9.81
                    phase = np.angle(eta_instant)
                    amp_aux = wave_ampw[nn]
                    pert_wave_aux = pert_wave[nn]
                    phase_aux = 1j*phase[[loc_y[kk]],[loc_x[kk]]]
                    pert_wave_coup = amp_aux*pert_wave_aux*np.exp(phase_aux)
                    tree.find(".//waveHeight").text = '%.6f'%(wave_ampw[nn]*2.0)
                    tree.find(".//generationType").text = '%.0f'%1 #generationType 1 in MW is coupling
                    #tree.find(".//ixsL").text = '%0f'%ixsL_new
                    #tree.find(".//ixsR").text = '%0f'%ixsR_new
                    jysT =  int(tree.find(".//ixsL").text)
                    tree.find(".//jysT").text = '%.0f'%jysT
                    tree.find(".//jysB").text = '%.0f'%jysT
                    tree.write(os.path.join(dir_cp_full[nn],"MILDwave.xml"))    
                    shutil.copyfile(os.path.join(dir_T,"depth.dat"),os.path.join(dir_cp_full[nn],"depth.dat"))
                    shutil.copyfile(os.path.join(dir_T,"obstacle.dat"),os.path.join(dir_cp_full[nn],"obstacle.dat"))

                    Lgmw = (int(tree.find(".//Nx").text)-1)*float(tree.find(".//dx").text)
                    Wgmw = (int(tree.find(".//Ny").text)-1)*float(tree.find(".//dy").text)
    
                    xFS = np.arange(-NEM_GRID['Lg']/2.0, NEM_GRID['Lg']/2.0+NEM_GRID['dx'],NEM_GRID['dy'])
                    xFS = xFS[:NEM_GRID['Nx']]
                    yFS = (-1)*np.arange(-NEM_GRID['Wg']/2.0, NEM_GRID['Wg']/2.0+NEM_GRID['dy'],NEM_GRID['dy'])
                    yFS = yFS[:NEM_GRID['Ny']]
                    #I multiply yFS times (-1) to get the correct values for the coordiante matrix in m
                    X,Y = np.meshgrid(xFS,yFS) #TODO figure out how to do for a non-circular grid
                    xFS_mw = np.arange(-Lgmw/2.0, Lgmw/2.0+NEM_GRID['dx'],NEM_GRID['dx'])
                    #xFS_mw = xFS_mw[:Nx_mw]
                    yFS_mw = (-1)*np.arange(-Wgmw/2.0, Wgmw/2.0+NEM_GRID['dy'],NEM_GRID['dy'])#there is no need to multiply here in order to obtain teh right input position for MILDwave
                    #yFS_mw = yFS_mw[:Ny_mw]                
                    Xmw,Ymw = np.meshgrid(xFS_mw,yFS_mw)
                    
                    if MW_CPL['circ']:
                        center_x = 0
                        center_y = 0
                        center_xmw = 0
                        center_ymw = 0
                        
                        mask1 = np.sqrt((X-center_x)**2+(Y-center_y)**2) < MW_CPL['radius'] + NEM_GRID['dx'] #it is necessary to add the extra  NEM_GRID['dx'] to obtain the maxium radius
                        mask2 = np.sqrt((Xmw-center_xmw)**2+(Ymw-center_ymw)**2) < MW_CPL['radius'] + NEM_GRID['dy']
        
                        yynm, xxnm = np.where(mask1)
                        yymw, xxmw = np.where(mask2)
        
                        
                    elif MW_CPL['rect']:  
                        xxh = int(NEM_GRID['Nx']/2) - int(MW_CPL['length']/(NEM_GRID['dx'] *2))
                        yyh = int(NEM_GRID['Ny']/2) - int(MW_CPL['width']/(NEM_GRID['dy']*2))
                        xxh_mw =  int(MW_ini['Nx']/2) - int(MW_CPL['length']/(NEM_GRID['dx']*2))
                        yyh_mw =  int(MW_ini['Ny']/2) - int(MW_CPL['width']/(NEM_GRID['dy']*2))
                        mask1 = np.zeros((np.shape(Y)[0],np.shape(X)[1]), dtype=bool)
                        mask1[yyh:yyh+int(MW_CPL['width']/(NEM_GRID['dy'])) ,xxh:xxh+ int(MW_CPL['length']/(NEM_GRID['dx']))] = True
                        mask2 = np.zeros((np.shape(Ymw)[0],np.shape(Xmw)[1]), dtype=bool)
                        mask2[yyh_mw:yyh_mw+int(MW_CPL['width']/(NEM_GRID['dy'])) ,xxh_mw:xxh_mw+int(MW_CPL['length']/(NEM_GRID['dx']))] = True
                        yynm, xxnm = np.where(mask1)
                        yymw, xxmw = np.where(mask2)
                        
                    amp = np.zeros((np.shape(xxnm)[0]))
                    phase = np.zeros((np.shape(xxnm)[0]))
                    amp[:] = np.abs(pert_wave_coup[yynm,xxnm])
                    phase[:] = np.angle(pert_wave_coup[yynm,xxnm])

        
                #-----------------------------------------------------------------------------
                #                                  5. CREATING NEMOH_DATA.TXT
                #-----------------------------------------------------------------------------
                    fid = open(os.path.join(dir_cp_full[nn],'nemoh_data.txt'),'w')  #one frequency
                    for mm in range(np.shape(amp)[0]):
                        fid.write('{:.0f}  {:.0f}  {:.8f}  {:.8f}\n'.format(xxmw[mm]+xoff[kk],yymw[mm]+yoff[kk],amp[mm],phase[mm]))
                    fid.close()        
        
def coupling_input_irregular(Nf,f,dir_Tp,MW_exe,MW_ini,MW_CPL):
    #load the dummy MILDwave.xml for the irregular wave case
    tree = et.parse(os.path.join(dir_Tp,"MILDwave.xml"))
    #obtain data of the wave gauge
    dx = np.float(tree.find(".//dx").text)*MW_ini['dx_bin_step']
    dy = np.float(tree.find(".//dy").text)*MW_ini['dy_bin_step']
    dx_ini = np.float(tree.find(".//dx").text)
    dy_ini = np.float(tree.find(".//dy").text)
    Tp = np.float(tree.find(".//wavePeriod").text)
    time_length = float(truncate(np.float(tree.find(".//etaOutput/end").text),3)) - float(truncate(np.float(tree.find(".//etaOutput/start").text),3)) - 30.0#remove ten seconds in case the number of etas is not the same for each simulation
    dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
    time_length = int(time_length/dt)*dt
    npoints = int((time_length/dt) * MW_ini['kd_x']  * MW_ini['kd_x'])
    neta = int(npoints/(MW_ini['kd_x']  *MW_ini['kd_y'] ))
    WG = []
    WG_aux = []
    [WG_aux.append(np.zeros(neta)) for locx in MW_CPL['ph_LOC_x']]
    WG.append(WG_aux)
    with open(os.path.join(dir_Tp,'EB','data','eta.npy')) as feta:
        eta_irr = np.fromfile(feta,dtype = np.float32,count = npoints,sep='')
    eta_size = MW_ini['kd_y']  * MW_ini['kd_x'] 
    eta_irr = [eta_irr[ii:ii+eta_size] for ii in range(0, len(eta_irr),eta_size)]
    eta_irr = [np.reshape(eta,(MW_ini['kd_y'] ,MW_ini['kd_x'] )) for eta in eta_irr]
    for tt in range(neta):
        for jj,locx,locy in zip(range(len(WG_aux)),MW_CPL['ph_LOC_x'],MW_CPL['ph_LOC_y']):
        
            WG[0][jj][tt] = eta_irr[tt][np.int(locy*dy_ini/dy),np.int(locx*dx_ini/dx)]
    del eta_irr
    fs = 1/dt#Sampling Frequency
    #FFT_BLOCK = 256 #256 normally gives a good result.Other used values 128, 256, 512 or 1024 . The higher Nf the higher FFT_BLOCK
    #Overlap   = 0.99 #Between 0 < Overlap < 1  normally a higher % in the overlap gives a better agreement between the target WAVE and the MW WAVE
    #https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.welch.html
    #nperseg should be N^2
    df = (np.max(f)-np.min(f))/(np.shape(f)[0]-1)
    S_coupling_loc = []
    Hi_coupling_loc = []
    amp_coupling_loc = []
    for ww in range(len(WG[0])):
        S_coupling_aux =  np.zeros(Nf)
        Hi_coupling_aux =  np.zeros(Nf)
        amp_coupling_aux =  np.zeros(Nf)
        fsd, S = signal.welch(WG[0][ww],fs,nperseg = MW_CPL['FFT_BLOCK'],noverlap = MW_CPL['FFT_BLOCK'] * MW_CPL['Overlap'],scaling = 'density')# density (returns the power spectral density) spectrum (returns the power spectrum)
        for ii in range(len(f)):
            for kk in range(1,np.shape(S)[0]):
                if fsd[kk-1] <= f[ii] <= fsd[kk]:
                    if f[ii]<(1/Tp):
                        S_coupling_aux[ii] = S[kk-1]
                        Hi_coupling_aux[ii] =  2.0 * np.sqrt(2.0*S_coupling_aux[ii]*df)
                        amp_coupling_aux[ii] = Hi_coupling_aux[ii]/2.0
                    else:
                        S_coupling_aux[ii] = S[kk]
                        Hi_coupling_aux[ii] =  2.0 * np.sqrt(2.0*S_coupling_aux[ii]*df)
                        amp_coupling_aux[ii] = Hi_coupling_aux[ii]/2.0
        S_coupling_loc.append(S_coupling_aux)
        Hi_coupling_loc.append(Hi_coupling_aux)
        amp_coupling_loc.append(amp_coupling_aux)
    return(amp_coupling_loc[0])# -*- coding: utf-8 -*-
def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])