"""
 This script estimates the numerical solution of the heat equation in one
 dimension. It does using a three point approximation for the double 
 spatial derivate.

 Numerical inputs:
   L       - The extension of the spatial grid 
   N       - The number of grid points
   Tfinal  - The duration of the simulation
   dt      - The step size in time

 Physical inputs:
   alpha   - The thermal diffusivity
   
 There is also an input for adjusting the pause between each frame displayed
 during the simulation.  
 All inputs, including the initial condition, are hard coded initially.
"""   

# Import libraries  
import numpy as np
from matplotlib import pyplot as plt

# Numerical grid parameters
L = 100
N = 50

# Numerical time parameters
Tfinal = 100
dt = 2

# Pause - used to adjust the speed of the "animation"
PauseVal = 0.02

# Thermal diffusivity
alpha = 1

# Set up grid
x = np.linspace(-L/2, L/2, N+1)
dx = L/N

# Set up initial distribution
uFunk0 = np.sqrt(L/2) - np.sqrt(np.abs(x))
#uFunk0 = 10*np.random.rand(N+1)

# Set up matrix for the double derivative - 3 point formula
# Allocate and declare
M = np.zeros((N+1, N+1), dtype=float)
# Endpoints
M[0, 0:2] = [-2, 1]
M[N, (N-1):(N+1)] = [1, -2]
# Interior points
for n in range (1,N):
  M[n, [n-1, n, n+1]] = [1, -2, 1]
# Correct pre-factors
M = M/dx**2

# Matrices for solving the heat eqaution
Imat = np.eye(N+1)
# Forward Euler
U_FE = Imat + dt*alpha*M
# Crank-Nicolson scheme
U_CN = np.matmul(np.linalg.inv(Imat - dt/2*alpha*M), Imat + dt/2*alpha*M)

# Initiate plots
plt.ion()
fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot()
# Initiate plot of solution
line1, = ax.plot(x, uFunk0, '-', color='blue', label = 'FE')
line2, = ax.plot(x, uFunk0, '--', color='red', label = 'CN')
plt.grid(visible=True)
plt.legend()

# Initiate wave functons and time
# Turn Psi0 into (N+1) \times 1 matrix (column vector)
uFunk_FE = np.matrix(uFunk0)                 
uFunk_FE = uFunk_FE.reshape(N+1,1)
uFunk_CN = np.matrix(uFunk0)                 
uFunk_CN = uFunk_CN.reshape(N+1,1)
t = 0

# Iterate while time t is less than final time Tfinal
while t < Tfinal:
  # Update numerical wave function
  uFunk_FE = np.matmul(U_FE, uFunk_FE)
  uFunk_CN = np.matmul(U_CN, uFunk_CN)
  
  # Update data for plots
  line1.set_ydata(uFunk_FE)
  line2.set_ydata(uFunk_CN)
  
  # Update plots
  fig.canvas.draw()
  #fig.canvas.flush_events()
  plt.pause(PauseVal)
  
  #Update time
  t = t+dt             