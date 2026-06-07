'''
CFD @ FIE
A.Y. 2025-2026
Professor: Edgardo A. Serafin
Student: Andrea Grimaldi

Partial Exam I
'''

import numpy as np
import matplotlib.pyplot as plt

########################################################################################################################
########################################################################################################################
# Starting parameters
xLen = 0.1 # meters
yLen = 0.1 # meters
xLenA = 0.04 # meters
yLenA = 0.04 # meters
height = 0.01 # meters
kA = 100 # W/m/K
kB = 550 # W/m/K
qA = 5 * 10**5 # W/m3
TS = 30 # Celsius
TWE = 15 # Celsius
hb = 200 # Wm2/K

# Derived parameters
nx = 20 # number of volumes in the x direction
ny = 20 # number of volumes in the y direction
assert nx % 5 == 0, "The number of volumes in x must be a multiple of 5"
assert ny % 5 == 0, "The number of volumes in y must be a multiple of 5"
dx = xLen / nx # length of a single volume in x
dy = yLen / ny # length of a single volume in y

# Gauss-Seidel parameters
iterMax = 30000
toler = 10**(-3)
########################################################################################################################
########################################################################################################################
# Finite volume discretization

# Conservation equations

# Padding

########################################################################################################################
########################################################################################################################
# GS setup

# Gauss-Seidel

########################################################################################################################
########################################################################################################################
# Post processing
