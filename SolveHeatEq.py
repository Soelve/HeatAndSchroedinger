"""
 This script estimates the numerical solution of the heat equation in one
 dimension. It does using a three point approximation for the double 
 spatial derivate.

 Numerical inputs:
   L       - The extension of the spatial grid 
   N       - The number of grid 
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
dt = 1

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

# Forward Euler scheme
Imat = np.eye(N+1)
U_FE = Imat + alpha*dt*M

# Initiate plots
plt.ion()
fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot()
# Initiate plot of solution
line1, = ax.plot(x, uFunk0, '-', color='black')
plt.grid(visible=True)


# Initiate wave functons and time
# Turn Psi0 into (N+1) \times 1 matrix (column vector)
uFunk = np.matrix(uFunk0)                 
uFunk = uFunk.reshape(N+1,1)
t = 0

# Iterate while time t is less than final time Tfinal
while t < Tfinal:
  # Update numerical wave function
  uFunk = np.matmul(U_FE, uFunk)
  
  # Update data for plots
  line1.set_ydata(uFunk)
  
  # Update plots
  fig.canvas.draw()
  #fig.canvas.flush_events()
  plt.pause(PauseVal)
  
  #Update time
  t = t+dt             