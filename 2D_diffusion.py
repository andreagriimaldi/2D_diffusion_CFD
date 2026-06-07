'''
CFD @ FIE
A.Y. 2025-2026
Professor: Edgardo A. Serafin
Student: Andrea Grimaldi

Partial Exam I
'''

import numpy as np
import matplotlib.pyplot as plt
from time import time


t0 = time()
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
nVol = dx * dy

# Helpers
def inSectorA(i, j):
    nAx = round(xLenA/dx)
    nAy = round(yLenA/dy)
    return (i < nAx) and (j >= ny - nAy)

def cellK(i, j):
    return kA if inSectorA(i, j) else kB

def faceK(i, j, ni, nj):
    kP = cellK(i, j)
    kNB = cellK(ni, nj)
    return (2 * kP * kNB)/(kP + kNB)


# Gauss-Seidel parameters
iterMax = 30000
toler = 10**(-3)

########################################################################################################################
########################################################################################################################
# Finite volume discretization
xp = np.linspace(dx/2, xLen-dx/2, nx)  # x coordinates for centroids
yp = np.linspace(dy/2, yLen-dy/2, ny)  # y coordinates for centroids

# Conservation equations
a_w = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        if i != 0:
            a_w[j, i] = (faceK(i, j, i - 1, j) * dy * height)/dx


a_e = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        if i != nx - 1:
            a_e[j, i] = (faceK(i, j, i + 1, j) * dy * height)/dx


a_n = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        if j != ny - 1:
            a_n[j, i] = (faceK(i, j, i, j + 1) * dx * height)/dy


a_s = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        if j != 0:
            a_s[j, i] = (faceK(i, j, i, j - 1) * dx * height)/dy
# Padding


########################################################################################################################
########################################################################################################################
# GS setup

# Gauss-Seidel


t1 = time()
time = t1 - t0

########################################################################################################################
########################################################################################################################
# Post processing
