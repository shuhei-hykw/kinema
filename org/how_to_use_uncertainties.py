# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 13:08:54 2016

@author: a
>pip install uncertainties
"""

from uncertainties import ufloat
from uncertainties.umath import *

# Calculate number of photons per pulse
#      P * tau
# N = ---------
#         E
#
# Where P -> power, tau -> pulse duration, E -> photon energy
 
tau = ufloat(3.2e-6, 0.2e-6)
# in seconds
wavelambda = ufloat(780, 0.25)
# in nanometers (note, "lambda" is a special word in python)
P = ufloat(3.2e-9, 0.5e-9)
# in watts
 
# photon energy is hc/lambda, hc = 1240 eV*nm or 1.986e-16 J*nm
 
N = (P * tau)/(1.986e-16/wavelambda)
 
print N




