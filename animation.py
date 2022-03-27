# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 15:38:01 2019

@author: gveraofe
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from xml.etree import ElementTree as et#Element Tree for modifying xml values
from matplotlib import cm
from matplotlib.patches import Circle

dir_MW = "C:/Users/gveraofe/Documents/RUN_MW/5_FLAP/Tp_08.00_dx_2.0_dout_4.0_nf_50_s1_18.37_t_4000_depth_10.0_FLAP_TEST_1_Lx_Ly_2000_VIDEO"
#dir_MW = 'E:/REAL_BATHYMETRY/Tp_08.00_dx_4.00_dout_8.00_nf_20_s1_0.0_t_600_depth_REAL_TEST_2_COMPLETE_NO_PERIODIC'
tree = et.parse(os.path.join(dir_MW,'EB',"MILDwave.xml"))
#TODO change the dict names to explicit
MW = np.load(os.path.join(dir_MW,'mw.npz'))
MW_sim = MW['MW_sim'].item()
MW_ini = MW['MW_ini'].item()
MW_exe = MW['MW_exe'].item()
MW_CPL = MW['MW_CPL'].item()
MW_OUT = MW['MW_OUT'].item()
   
kd_x = 501
kd_y = 501
 
#obtain data of the wave gauge
dx = np.float(tree.find(".//dx").text)*MW_exe['dxstep']
dy = np.float(tree.find(".//dy").text)*MW_exe['dystep']
dx_ini = np.float(tree.find(".//dx").text)
dy_ini = np.float(tree.find(".//dy").text)
Tp = np.float(tree.find(".//wavePeriod").text)
time_length = float(truncate(np.float(tree.find(".//etaOutput/end").text),3)) - float(truncate(np.float(tree.find(".//etaOutput/start").text),3)) - 30.0#remove ten seconds in case the number of etas is not the same for each simulation
dt = np.float(tree.find(".//etaOutput/increment").text)*np.float(tree.find(".//delt").text)
time_length = int(time_length/dt)*dt
npoints = int((time_length/dt) * kd_x * kd_y)
neta = int(npoints/(kd_x *kd_y))

eta_size = kd_y * kd_x

eta_irr_cp = np.load(os.path.join(dir_MW,'EB','data','eta.npy'))
eta_irr_cp = [eta_irr_cp[ii:ii+eta_size] for ii in range(0, len(eta_irr_cp),eta_size)]
eta_irr_cp = [np.reshape(eta,(kd_y,kd_x)) for eta in eta_irr_cp]

import pylab as plt
import numpy
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = r"C:/some_path/ffmpeg.exe"   # if necessary

# Generate data for plotting
Lx = Ly = 3
Nx = Ny = 11
Nt = 116
x = MW_OUT['x_vector']
y = MW_OUT['y_vector']
Xmin = 0
Xmax = 501
Ymin = 0
Ymax = 501

elev_min= 0.85
elev_max= 1.1501
mid_val=  1.0

levels = np.arange(elev_min, elev_max, 0.025)
cmap = cm.ocean

fig = plt.figure()
ax = plt.axes(xlim=(Xmin, Xmax), ylim=(Ymin, Ymax), xlabel='Length of the Basin', ylabel='Width of the Basin')

cvals = numpy.linspace(-2.0,2.0,20+1)      # set contour values 
cont = plt.contourf(eta_irr_cp[0], cvals,cmap=cm.get_cmap(cmap, len(cvals) - 1))    # first image on screen
#ax.add_patch(Circle((0,0),20.0, facecolor="white"))#ls = '--'
plt.colorbar()

# animation function
def animate(i):
    global cont
    z = eta_irr_cp[i]
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf( z, cvals,cmap=cm.get_cmap(cmap, len(cvals) - 1))
    #ax.add_patch(Circle((0,0),20.0, facecolor="white"))#ls = '--'
    plt.title('Simulation Time = %i s' % (i*dt))
    return cont

anim = animation.FuncAnimation(fig, animate, frames=Nt,interval = 10, repeat=False)

#anim.save('animation.mp4', writer=animation.FFMpegWriter())

def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])