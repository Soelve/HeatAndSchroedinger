"""
 This script simulates the evolution of a Gaussian wave packet moving
 freely in one dimension. It does so by estimating the numerical solution 
 of the Schrödinger equation using a three point approximation for the 
 kinetic energy operator.

 Numerical inputs:
   L       - The extension of the spatial grid 
   N       - The number of grid points
   Tfinal  - The duration of the simulation
   dt      - The step size in time

 Physical inputs:
   x0      - The mean position of the initial wave packet
   p0      - The mean momentum of the initial wave packet
   sigmaP  - The momentum width of the initial, Gaussian wave packet
 
 There is also an input for adjusting the pause between each frame displayed
 during the simulation.  
 All inputs are hard coded initially.
"""   

# Import libraries  
import numpy as np
from matplotlib import pyplot as plt

# Numerical grid parameters
L = 100
N = 150

# Numerical time parameters
Tfinal = 20
dt = 0.5

# Pause - used to adjust the speed of the "animation"
PauseVal = 0.02

# Inputs for the Gaussian 
x0 = -20
p0 = 3
sigmaP = 0.2


# Set up grid
x = np.linspace(-L/2, L/2, N+1)
dx = L/N

# Set up matrix for the double derivative - 3 point formula
# Allocate and declare
H = np.zeros((N+1, N+1), dtype=complex)
# Endpoints
H[0, 0:2] = [-2, 1]
H[N, (N-1):(N+1)] = [1, -2]
# Interior points
for n in range (1,N):
  H[n, [n-1, n, n+1]] = [1, -2, 1]
# Correct pre-factors
H = H/dx**2
H = -1/2*H

# Crank-Nicolson propagator
Imat = np.eye(N+1)
U_CN = np.matmul(np.linalg.inv(Imat + 1j*dt/2*H), Imat - 1j*dt/2*H)

# Set up initial Gaussian - analytically
InitialNorm = np.power(2/np.pi, 1/4) * np.sqrt(sigmaP)
Psi0 = InitialNorm*np.exp(-sigmaP**2*(x-x0)**2+1j*p0*x)

# Initiate plots
plt.ion()
fig = plt.figure(1)
plt.clf()
ax = fig.add_subplot()
# Initiate plot of wave function
line1, = ax.plot(x, np.abs(Psi0)**2, '-', color='black')
plt.grid(visible=True)


# Initiate wave functons and time
# Turn Psi0 into (N+1) \times 1 matrix (column vector)
Psi = np.matrix(Psi0)                 
Psi = Psi.reshape(N+1,1)
t = 0

# Iterate while time t is less than final time Tfinal
while t < Tfinal:
  # Update numerical wave function
  Psi = np.matmul(U_CN, Psi)
  
  # Update data for plots
  line1.set_ydata(np.power(np.abs(Psi),2))
  
  # Update plots
  fig.canvas.draw()
  #fig.canvas.flush_events()
  plt.pause(PauseVal)
  
  #Update time
  t = t+dt             