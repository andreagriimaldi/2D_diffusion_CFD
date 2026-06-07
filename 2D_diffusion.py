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
hb = 200 # W/m2/K

# Derived parameters
nx = 20 # number of volumes in the x direction
ny = 20 # number of volumes in the y direction
assert nx % 5 == 0, "The number of volumes in x must be a multiple of 5"
assert ny % 5 == 0, "The number of volumes in y must be a multiple of 5"
dx = xLen / nx # length of a single volume in x
dy = yLen / ny # length of a single volume in y
nVol = nx * ny # number of volumes

# Helpers
config = 1   # 1 for case 1 (top-left), 2 for case 2 (bottom right)
def inSectorA(i, j):
    nAx = round(xLenA/dx)
    nAy = round(yLenA/dy)
    if config == 1:
        return (i < nAx) and (j >= ny - nAy)
    else:
        return (i >= nx - nAx) and (j < nAy)

def cellK(i, j):
    return kA if inSectorA(i, j) else kB

def faceK(i, j, ni, nj):
    kP = cellK(i, j)
    kNB = cellK(ni, nj)
    return (2 * kP * kNB)/(kP + kNB)


# Gauss-Seidel parameters
iterMax = 30000
toll = 10**(-6)

########################################################################################################################
########################################################################################################################
# Finite volume discretization
xp = np.linspace(dx/2, xLen-dx/2, nx)  # x coordinates for centroids
yp = np.linspace(dy/2, yLen-dy/2, ny)  # y coordinates for centroids

# Building coefficients
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


b = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        if inSectorA(i, j):
            b[j, i] += qA * dx * dy * height


a_p = np.zeros((ny, nx))
for j in range(ny):
    for i in range(nx):
        a_p[j, i] = a_w[j, i] + a_e[j, i] + a_n[j, i] + a_s[j, i]


for i in range(nx):
    a_b = cellK(i, 0) * dx * height / (dy/2)  # Dirichlet boundary condition
    a_p[0, i] += a_b
    b[0, i] += a_b * TS


for j in range(ny):
    a_rob = (dy * height)/((dx/(2 * cellK(0, j))) + (1/hb))  # Robin boundary condition (West face)
    a_p[j, 0] += a_rob
    b[j, 0] += a_rob * TWE


for j in range(ny):
    a_rob = (dy * height)/((dx/(2 * cellK(nx - 1, j))) + (1/hb))  # Robin boundary condition (East face)
    a_p[j, nx - 1] += a_rob
    b[j, nx - 1] += a_rob * TWE


# Scarborough criterion
neighbors = a_w + a_e + a_n + a_s
assert np.all(a_p > 0)
assert np.all(a_p >= neighbors - 1e-9)
assert np.any(a_p > neighbors + 1e-9)

# Padding
a_w = np.pad(a_w, 1)
a_e = np.pad(a_e, 1)
a_n = np.pad(a_n, 1)
a_s = np.pad(a_s, 1)
b = np.pad(b, 1)
a_p = np.pad(a_p, 1)

########################################################################################################################
########################################################################################################################
# GS setup
T = np.zeros((ny + 2, nx + 2))
iter = 0

# Gauss-Seidel
for k in range(0, iterMax):
    Told = np.copy(T)
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            RHS = a_w[j, i] * T[j, i - 1] + a_e[j, i] * T[j, i + 1] + a_n[j, i] * T[j + 1, i] + a_s[j , i] * T[j - 1, i] + b[j, i]
            T[j, i] = RHS / a_p[j, i]
    iter += 1
    norm = np.linalg.norm(Told - T)
    if norm < toll:
        break


T_final = T[1:ny+1, 1:nx+1].copy()

t1 = time()
time = t1 - t0

########################################################################################################################
########################################################################################################################
# Validating the solution
R = np.zeros((ny+2, nx+2))
for j in range(1, ny+1):
    for i in range(1, nx+1):
        nb = (a_w[j, i] * T[j, i-1] + a_e[j, i] * T[j, i+1] + a_n[j,i] *T [j+1, i] + a_s[j,i] * T[j-1, i] + b[j, i])
        R[j,i] = nb - a_p[j, i] * T[j, i]
assert np.abs(R).max() < 1e-6, "Discrete equation not resolved for all the cells" # The discrete equation must be satisfied for each cell (up to the calculation's precision)


for j in range(ny):
    for i in range(nx-1):
        assert abs(a_e[j + 1, i + 1] - a_w[j + 1, i + 2]) < 1e-9, f"E/W face mismatch at ({i},{j})" # Vertical face mismatch

for j in range(ny-1):
    for i in range(nx):
        assert abs(a_n[j + 1, i + 1] - a_s[j + 2, i + 1]) < 1e-9, f"N/S face mismatch at ({i},{j})" # Horizontal face mismatch

########################################################################################################################
########################################################################################################################
# Energy calculation


########################################################################################################################
########################################################################################################################
# Plotting figure
plt.figure(1)
X,Y = np.meshgrid(xp,yp)
plt.contourf(X,Y,T_final,30)
if config == 1:
    plt.title('Figure 1', fontsize=15)
else:
    plt.title('Figure 2', fontsize=15)
plt.xlabel('longitude x [m]')
plt.ylabel('longitude y [m]')
plt.xlim(0.0, xLen)
plt.ylim(0.0, yLen)
plt.grid(True)
plt.colorbar()
plt.set_cmap('jet')
plt.show()

# Print onscreen
Tmax = T_final.max()
Tmin = T_final.min()
if config == 1:
    print('\n------------------ First part (Figure 1) --------------------')
else:
    print('\n------------------ Second part (Figure 2) --------------------')
print('                   n. of volumes = {} '.format(nVol))
print('                   dx = {:1.4f} m'.format(dx))
print('                   dy = {:1.4f} m'.format(dy))
print('                   T max = {:1.2f} celsius'.format(Tmax))
print('                   T min = {:1.2f} celsius'.format(Tmin))
print('                   n. of iterations = {} '.format(iter))
print('                   Machine time = {:1.1f} sec'.format(time))