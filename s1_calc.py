# CODE TO CALCULATE S1 FROM EXISTING SIGMA THETA OR S
"""
Created on Tue Nov 20 11:56:11 2018
@author: gveraofe
"""

import numpy as np
import scipy.special as spe
from scipy.optimize import fsolve

s = 75

sigmatheta = np.sqrt(2.0/(s+1.0))

func = lambda s1 : sigmatheta - np.sqrt( 2 - (2* (spe.gamma(s1+1))**2)/(spe.gamma(s1+0.5) * spe.gamma(s1 +1.5))) 

s1_ini = 15.8

s1_sol = fsolve(func,s1_ini)